"""Calibrate context terms on known families: how often does a term in a family's context
predict that the family itself carries it?

Run after the full pipeline run; it is not a pipeline rule. For each term (amr:x,
defence:Clover ...) in 12_context_and_structure/family_context_terms.tsv:

  benchmark  known families whose own members carry the term: more than half of the
             member proteins have a label giving the term, or, for a system term, have
             an ORF that is a component of that system
  negatives  known families that carry at least one term by that majority rule (they are
             labelled) and of which NO member carries this term
  excluded   families where a minority of members carry the term, and families without
             a SUCCESS row (fewer than 2 lineages), whose conservation is not measured

The predictor is the family's context conservation of the term (0 when it has no row for
it). At each observed conservation level c > 0, precision = benchmark families at or above
c / all evaluated families at or above c. The threshold for a precision level is the lowest
c whose precision reaches the level and stays there at every higher c that holds enough
families to measure it.

MIN_FAMILIES = 10. A precision of 90% can be observed only on at least 10 families: on
fewer, one false positive already gives 8/9 = 89%. So a term needs 10 benchmark families
and 10 negatives (with no negatives every level is 100% precise by construction) to be
CALIBRATED, and a threshold needs 10 families at or above it. Otherwise the term is
UNCALIBRATED and its thresholds stay blank; the curve is written either way.

Outputs: --out one row per term with the thresholds, --curve the precision curve, and
--summary a short text summary.

Usage:
  python tools/calibrate_context.py \\
      --terms results/12_context_and_structure/family_context_terms.tsv \\
      --labels results/08_protein_labels/protein_labels.tsv \\
      --families results/10_clustering/protein_families.tsv \\
      --map results/03_dereplication/protein_map.tsv \\
      --defence results/12_context_and_structure/defence_systems.tsv \\
      --conjugation results/12_context_and_structure/conjugation_systems.tsv \\
      --out context_calibration.tsv --curve context_calibration_curve.tsv \\
      --summary context_calibration.txt
"""
import argparse
import collections
import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from plasmidann.context_terms import SUCCESS, label_term, system_term  # noqa: E402

MIN_FAMILIES = 10
LEVELS = (0.5, 0.9)
csv.field_size_limit(sys.maxsize)


def precision_curve(scores, positives):
    """Precision at each distinct conservation level above 0, ascending.

    `scores` maps evaluated families (benchmark and negatives) to their conservation; a
    family left out scores 0 and is predicted at no level. One pass from the top keeps
    this linear after the sort, since a term can have ~10^5 negatives.
    """
    ranked = sorted(((s, f in positives) for f, s in scores.items() if s > 0),
                    reverse=True)
    curve, n_predicted, n_true = [], 0, 0
    for i, (c, positive) in enumerate(ranked):
        n_predicted += 1
        n_true += positive
        if i + 1 == len(ranked) or ranked[i + 1][0] != c:
            curve.append({"conservation": c, "n_predicted": n_predicted,
                          "n_true": n_true, "precision": n_true / n_predicted})
    return curve[::-1]


def threshold(curve, level, min_families=MIN_FAMILIES):
    """The lowest curve point reaching `level` that no measurable higher point falls below.

    Points holding fewer than `min_families` families are too small to measure precision
    at `level`, so they neither qualify nor disqualify.
    """
    best = None
    for point in sorted(curve, key=lambda p: -p["conservation"]):
        if point["n_predicted"] < min_families:
            continue
        if point["precision"] < level:
            break
        best = point
    return best


def _own_terms(args, family_members):
    """term -> {family_id: fraction of member proteins carrying it}."""
    wanted = {m for members in family_members.values() for m in members}
    terms_of = collections.defaultdict(set)
    with open(args.labels, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["protein_id"] in wanted:
                term = label_term(r["kind"], r["label"], r.get("sub_label", ""))
                if term:
                    terms_of[r["protein_id"]].add(term)
    system_of = collections.defaultdict(set)
    for path, prefix in ((args.defence, "defence"), (args.conjugation, "conj")):
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("orf_id"):
                    system_of[r["orf_id"]].add(system_term(prefix, r["system"]))
    if system_of:
        with open(args.map) as fh:
            for line in fh:
                sid, ids = line.rstrip("\n").split("\t")
                if sid in wanted:
                    for oid in ids.split(","):
                        terms_of[sid] |= system_of.get(oid, set())
    fraction = collections.defaultdict(dict)
    for family_id, members in family_members.items():
        counts = collections.Counter(t for m in members for t in terms_of.get(m, ()))
        for term, n in counts.items():
            fraction[term][family_id] = n / len(members)
    return fraction


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for name in ("terms", "labels", "families", "map", "defence", "conjugation", "out",
                 "curve", "summary"):
        ap.add_argument(f"--{name}", required=True)
    args = ap.parse_args(argv)

    # Known families with a measured conservation, and their context scores per term.
    score = collections.defaultdict(dict)
    evaluable = set()
    with open(args.terms, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["family_set"] != "known" or r["status"] != SUCCESS:
                continue
            evaluable.add(r["family_id"])
            score[r["term"]][r["family_id"]] = float(r["conservation"])

    family_members = {}
    with open(args.families, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["family_id"] in evaluable:
                family_members[r["family_id"]] = r["members"].split(",")
    fraction = _own_terms(args, family_members)
    labelled = {f for per_family in fraction.values()
                for f, x in per_family.items() if x > 0.5}

    rows, curve_rows = [], []
    for term in sorted(set(score) | set(fraction)):
        carried = fraction.get(term, {})
        positives = {f for f, x in carried.items() if x > 0.5}
        negatives = labelled - set(carried)
        scores = {f: s for f, s in score[term].items() if f in positives or f in negatives}
        curve = precision_curve(scores, positives)
        curve_rows += [{"term": term, **p} for p in curve]
        calibrated = len(positives) >= MIN_FAMILIES and len(negatives) >= MIN_FAMILIES
        row = {"term": term, "term_type": term.split(":", 1)[0],
               "n_benchmark": len(positives), "n_negative": len(negatives),
               "status": "CALIBRATED" if calibrated else "UNCALIBRATED"}
        for level in LEVELS:
            tag = int(level * 100)
            point = threshold(curve, level) if calibrated else None
            row[f"threshold_{tag}"] = point["conservation"] if point else ""
            row[f"precision_{tag}"] = round(point["precision"], 4) if point else ""
            row[f"n_at_{tag}"] = point["n_predicted"] if point else ""
        rows.append(row)

    with open(args.out, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=list(rows[0]) if rows else ["term"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(args.curve, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=["term", "conservation", "n_predicted",
                                            "n_true", "precision"], delimiter="\t")
        w.writeheader()
        w.writerows({**p, "precision": round(p["precision"], 4)} for p in curve_rows)

    calibrated = [r for r in rows if r["status"] == "CALIBRATED"]
    lines = [f"context calibration: {len(rows)} terms, {len(calibrated)} CALIBRATED, "
             f"{len(rows) - len(calibrated)} UNCALIBRATED (fewer than {MIN_FAMILIES} "
             f"benchmark or negative families); {len(evaluable)} known families with a "
             "measured conservation",
             "term\tn_benchmark\tn_negative\tthreshold_50\tthreshold_90"]
    lines += [f"{r['term']}\t{r['n_benchmark']}\t{r['n_negative']}\t"
              f"{r['threshold_50'] or '-'}\t{r['threshold_90'] or '-'}" for r in calibrated]
    lines.append("UNCALIBRATED: " + ", ".join(
        r["term"] for r in rows if r["status"] == "UNCALIBRATED"))
    pathlib.Path(args.summary).write_text("\n".join(lines) + "\n")
    print(lines[0])


if __name__ == "__main__":
    main()
