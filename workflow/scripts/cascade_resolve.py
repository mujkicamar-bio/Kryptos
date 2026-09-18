"""S3 final step: turn every tier's hits into one row per protein.

This is where the cascade's evidence becomes a decision. Two principles govern it:

  * The class is a property of the PROTEIN. Coverage is merged across all informative hits
    from all tiers before anything is decided. Classifying per hit put 17% of Pfam-hit
    proteins into DOMAIN_ONLY at >=50% explained - and since DOMAIN_ONLY is
    target-eligible, complete replication initiators and T4SS ATPases were walking into
    the screening pool.

  * min_explained is applied HERE, post hoc, not during the search. The search narrowed on
    narrow_at (0.9), which is permissive, so every protein in the interesting band was seen
    by every tier and this threshold can be swept without re-running anything.

Two columns describe coverage, over disjoint evidence, and they are not redundant:

  explained_fraction   merged INFORMATIVE spans / length. Informative about ANNOTATED
                       proteins. Constant 0 across the dark set by construction.
  dark_covered_fraction merged UNINFORMATIVE spans / length. The discriminating axis for
                       the dark set: 95% covered by 'hypothetical protein' across three
                       databases is a real, conserved, full-length unnamed protein; a
                       single 20-aa fragment hit is not.
"""
import _ctx  # noqa: F401
import collections
import csv

from plasmidann.cascade import (classify, completeness, dark_evidence, is_informative,
                                explained_fraction, uninformative_spans, n_dark_databases,
                                check_thresholds)

cfg = snakemake.params.thresholds
# Fail loudly and early on an incoherent threshold block rather than producing a table
# nobody can interpret. The JSON schema checks ranges; this checks relationships.
check_thresholds(cfg)
tier_order = snakemake.params.tier_order

# ------------------------------------------------------------------------------------
# Gather every hit for every protein, from every tier.
# ------------------------------------------------------------------------------------
by_query = collections.defaultdict(list)
unnamed = collections.defaultdict(list)
named = collections.defaultdict(list)
for f in snakemake.input.hits:
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["coverage"] = float(r["coverage"]) if r["coverage"] else ""
            r["target_coverage"] = float(r["target_coverage"]) if r["target_coverage"] else None
            # Empty for a tier that reports no alignment span (pharokka). Left as the
            # empty string rather than coerced: cascade.classify reads '' as "no span".
            r["start"] = int(r["start"]) if r["start"] else ""
            r["end"] = int(r["end"]) if r["end"] else ""
            by_query[r["query"]].append(r)
            if not is_informative(r["label"]):
                unnamed[r["query"]].append(r)
            else:
                named[r["query"]].append(r)

# The informative explained fraction is taken from the final tier's cumulative spans,
# which were merged as the cascade descended. qlen comes with it.
#
# One file per cascade shard. Each tier writes forward everything it inherited, so the last
# tier's spans hold every protein the cascade ever saw, and a protein appears in exactly one
# shard - it never crosses. Reading only the first file would zero the explained fraction of
# every protein in every other shard, and a zero explained fraction reads as "nothing named
# it", which would push the entire plasmid backbone into the screening pool.
explained, qlen = {}, {}
span_files = snakemake.input.spans
for path in ([span_files] if isinstance(span_files, str) else list(span_files)):
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            explained[r["seq_id"]] = float(r["explained_fraction"])
            qlen[r["seq_id"]] = int(r["qlen"])

seq_ids = [l[1:].split()[0] for l in open(snakemake.input.faa) if l[0] == ">"]

cols = [
    "seq_id",
    # what it is
    "annot_tier", "annot_label", "functional_class", "homology_depth",
    "annot_qcov", "annot_tcov", "annot_evalue", "n_informative_hits",
    # EVERY informative label, not just the winning one. The backbone stop-list matches
    # Pfam family names, and only the Pfam tiers emit those; ranking labels by E-value
    # means an nr hit with free text routinely takes annot_label away from a curated Pfam
    # assignment on the same protein. A complete replication initiator carrying RepA_N at
    # T1 then reached the screening pool with the guard that exists to stop it never
    # firing, because the only column S5 could read no longer held a Pfam name.
    "informative_labels", "informative_tiers",
    # how much of it is accounted for
    "explained_fraction", "annot_completeness", "meets_min_explained",
    # what the dark evidence says
    "dark_covered_fraction", "dark_completeness", "dark_evidence",
    "n_dark_databases", "uninformative_labels", "uninformative_tiers",
    # provenance: every row carries the thresholds that produced it (P4)
    "thr_min_coverage", "thr_min_explained", "thr_narrow_at", "thr_full_at", "thr_partial_at",
]

with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for sid in seq_ids:
        ef = explained.get(sid, 0.0)
        length = qlen.get(sid, 0)
        u = unnamed.get(sid, [])
        n = named.get(sid, [])
        labels = [x["label"] for x in u]

        # Dark coverage: the same merge, over the spans of hits that named nothing.
        dcf = explained_fraction(length, uninformative_spans(
            [x for x in u if x["start"] != "" and x["end"] != ""]))

        classified = classify(by_query.get(sid, []), explained=ef,
                              min_coverage=cfg["min_coverage"], tier_order=tier_order)
        span_measured = classified.pop("span_measured")

        w.writerow({
            "seq_id": sid,
            **classified,
            "explained_fraction": ef,
            # NOT_MEASURED when the class rests on a family-level assignment alone: the
            # explained fraction is then 0 because nothing measured it, and reporting
            # completeness NONE would read as "nothing matched" on a protein whose whole
            # family is known. Spec section 2.9: absence of a measurement is a status.
            "annot_completeness": (
                completeness(ef, full_at=cfg["full_at"], partial_at=cfg["partial_at"])
                if span_measured else "NOT_MEASURED"),
            # min_explained as a reported flag rather than a filter: the protein is in the
            # table either way, and this column can be recomputed at any threshold.
            "meets_min_explained": int(ef >= cfg["min_explained"]),
            "dark_covered_fraction": dcf,
            "dark_completeness": completeness(dcf, full_at=cfg["full_at"],
                                              partial_at=cfg["partial_at"]),
            "dark_evidence": dark_evidence(labels),
            "n_dark_databases": n_dark_databases(u),
            # The labels themselves are kept: someone else has described this protein, and
            # that is evidence FOR screening it, not against.
            "informative_labels": " | ".join(dict.fromkeys(x["label"] for x in n)),
            "informative_tiers": ",".join(dict.fromkeys(x["tier"] for x in n)),
            "uninformative_labels": " | ".join(dict.fromkeys(labels)),
            "uninformative_tiers": ",".join(dict.fromkeys(x["tier"] for x in u)),
            "thr_min_coverage": cfg["min_coverage"],
            "thr_min_explained": cfg["min_explained"],
            "thr_narrow_at": cfg["narrow_at"],
            "thr_full_at": cfg["full_at"],
            "thr_partial_at": cfg["partial_at"],
        })

print(f"resolved {len(seq_ids)} proteins across tiers {tier_order}")
