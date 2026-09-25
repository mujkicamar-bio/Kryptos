"""S5: the quality gate. Halts the run if the cascade has a recall problem.

TWO CONTROLS, ONE OF THEM HALTING

The positive control halts the run; the negative control is measured and reported. That
asymmetry is deliberate. A recall failure means "unannotated" carries no information, so
everything downstream is void and the run should stop. A decoy that slipped through is a
number the reader needs in order to interpret the dark set, not a reason to discard a
collection - and a halting negative gate would stop the pipeline over the hardest
sequences in it.

POSITIVE CONTROL (success criterion SC2)
    Known plasmid biology is run through the same cascade as everything else and must come
    out FUNCTIONAL. ECLIPSE recovered 99.2-100% of 246 virulence, 42 AMR and 75 essential
    genes this way. Anything KNOWN that emerges dark is a recall failure, and a target list
    built on a broken annotation step is worse than no target list - it costs bench time
    and it is not detectably wrong until the assays fail.

    Its absence from v1 was the single point of unanimous reviewer criticism.

WHY THERE IS NO LONGER A BACKBONE STOP-LIST

An earlier version also excluded any protein carrying a curated plasmid-backbone family
name, as a second guard behind the classification. It has been removed deliberately.

The pipeline is two steps: annotate everything that can be annotated, then report against
criteria. If step one works, a complete replication initiator is FUNCTIONAL and is simply
not dark - it needs no list to keep it out. If step one is broken, a hand-maintained list of
family names is the wrong repair, and a dangerous one: an exclusion by name removes the
protein from the record entirely, so its removal is invisible and unauditable. That is not
theoretical - nine of the seventy-three names on that list did not exist in Pfam-A at all,
so those rules had never fired, and nobody could tell from any output.

The curated family table is gone. Pfam-A 38.2 holds 30,134 families, of which 67 mention
replication in their description and 42 mention conjugation; the list named 16 and 15, and
nine of its 73 names did not exist in Pfam-A at all. Functional labels now come from the
tools themselves, in results/08_protein_labels/protein_labels.tsv - and they LABEL a
protein that stays in the table rather than excluding it. Nothing excludes anything by name.

Nothing is deleted here (P5). Proteins are flagged, and the flags are counted.
"""
import collections
import csv

import _ctx  # noqa: F401

from plasmidann.controls import control_recall

cfg = snakemake.params.gate

# ------------------------------------------------------------------------------------
# Read the cascade output.
# ------------------------------------------------------------------------------------
rows = []
with open(snakemake.input.prot, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        rows.append(r)

# ------------------------------------------------------------------------------------
# Positive control: reviewed Swiss-Prot proteins of known function, spiked into the query
# set at S2c and carried through every tier exactly as a real protein is.
#
# Independent by construction. An earlier version used proteins the cascade itself had
# labelled as backbone, which could not detect the failure that matters - a protein the
# cascade MISSED never enters a self-drawn control set.
# ------------------------------------------------------------------------------------
CONTROL_PREFIX = "CTRL_"
control_rows = [r for r in rows if r["seq_id"].startswith(CONTROL_PREFIX)]

if not control_rows and cfg.get("require_control_set", True):
    raise SystemExit(
        "S5 QUALITY GATE: no control proteins found in the annotation table.\n"
        "The controls are spiked in at S2c and must reach S3. A gate that silently skips "
        "itself is not a gate - check that results/03_dereplication/cascade_input.faa was used as the "
        "cascade query set.")

recall = control_recall(control_rows) if control_rows else 0.0
failed = [r["seq_id"] for r in control_rows if r["functional_class"] != "FUNCTIONAL"]

# WHICH TIER resolved each control, not just whether one did.
#
# The controls are Swiss-Prot proteins and one tier searches swissprot.dmnd, so they
# self-hit there at essentially perfect identity. That is a real annotation path and not a bug, but
# it makes overall recall a weak test: the gate would pass even if the curated Pfam tiers
# were completely broken, because the Swiss-Prot tier would rescue every control on its own.
#
# Reporting resolution per tier turns a nearly-trivial pass into a diagnostic. If controls
# only ever resolve at the Swiss-Prot tier, T1 and T2 have a recall problem that overall
# recall hides. Which tier that is comes from the cascade configuration (T3 in an earlier
# layout, T4 now), so the report names each tier's database rather than assuming one.
by_tier = collections.Counter(r.get("annot_tier") or "UNRESOLVED" for r in control_rows)
shallow = sum(n for t, n in by_tier.items() if t in ("T1", "T2"))
shallow_fraction = round(shallow / len(control_rows), 4) if control_rows else 0.0

# ------------------------------------------------------------------------------------
# Negative control: decoys built at S2d from real plasmid CDS, searched by every tier
# under the same thresholds as everything else (spec section 58.2).
#
# Expected behaviour is DARK. A decoy classed FUNCTIONAL is a false positive of the
# annotation cascade - the cascade named something that is not a protein - and since the
# deliverable is the complement of what the cascade could name, that rate is what tells a
# reader how clean the complement is.
#
# Reported, never halting. See the two-controls note at the top.
# ------------------------------------------------------------------------------------
DECOY_PREFIX = "DECOY_"
decoy_rows = [r for r in rows if r["seq_id"].startswith(DECOY_PREFIX)]
decoy_named = [r for r in decoy_rows if r["functional_class"] == "FUNCTIONAL"]
decoy_fpr = (round(len(decoy_named) / len(decoy_rows), 4) if decoy_rows else "")

# Broken down by construction, because they fail for different reasons. A shuffled decoy
# that gets named was named on COMPOSITION alone; a reverse-complement decoy that gets
# named is the shadow-ORF artefact the QC stage exists to flag. One rate hides which.
by_class = collections.Counter(
    "shuffled" if "_shuf_" in r["seq_id"] else "reverse_complement" for r in decoy_rows)
named_by_class = collections.Counter(
    "shuffled" if "_shuf_" in r["seq_id"] else "reverse_complement" for r in decoy_named)
# And by the tier that named each one: a decoy named at nr is a false positive of the free
# text tier specifically, and nothing else in the gate measures that tier.
named_by_tier = collections.Counter(r.get("annot_tier") or "" for r in decoy_named)

# ------------------------------------------------------------------------------------
# Artefact flags from S2b and edge-partial ORFs are the other two exclusions.
# ------------------------------------------------------------------------------------
artefact_ids = set()
with open(snakemake.input.artefact, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["artefact_flag"] == "1":
            artefact_ids.add(r["seq_id"])

cols = ["seq_id", "is_artefact", "target_eligible", "exclusion_reason"]
n_eligible = 0
# The dark set under two definitions (spec section 25). Ours counts a protein whose only
# homologues are themselves unnamed (UNCHARACTERIZED_HOMOLOG) as dark; FESNov
# (Rodriguez del Rio et al. 2024, Nature 626:377) calls a family unknown only when it has
# no homologue at all, which here is functional_class NONE. Both counts are reported.
n_dark_strict = 0
with open(snakemake.output.flags, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in rows:
        sid = r["seq_id"]
        if sid.startswith(CONTROL_PREFIX) or sid.startswith(DECOY_PREFIX):
            continue          # instrumentation, not a candidate
        reasons = []
        if sid in artefact_ids:
            reasons.append("artefact")
        # Only proteins nothing could name are screening candidates at all. A protein the
        # selection did not search (S2s) is neither named nor dark, and says so.
        if r["functional_class"] == "NOT_SEARCHED":
            reasons.append("not_searched")
        elif r["functional_class"] not in ("UNCHARACTERIZED_HOMOLOG", "NONE"):
            reasons.append("annotated")
        eligible = int(not reasons)
        n_eligible += eligible
        n_dark_strict += int(eligible and r["functional_class"] == "NONE")
        w.writerow({"seq_id": sid,
                    "is_artefact": int(sid in artefact_ids),
                    "target_eligible": eligible,
                    "exclusion_reason": ",".join(reasons)})

summary = (f"proteins={len(rows)} artefact={len(artefact_ids)} "
           f"target_eligible={n_eligible} dark_no_homologue={n_dark_strict} "
           f"control_n={len(control_rows)} "
           f"control_recall={recall} decoy_n={len(decoy_rows)} "
           f"decoy_false_positive_rate={decoy_fpr}")
with open(snakemake.output.report, "w") as out:
    out.write(summary + "\n")
    out.write(f"min_control_recall={cfg['min_control_recall']}\n")
    out.write("\ncontrols resolved per tier (controls self-hit at the swissprot tier):\n")
    sources = snakemake.params.tier_sources
    for tier, n in sorted(by_tier.items()):
        out.write(f"  {tier}\t{sources.get(tier, '')}\t{n}\n")
    out.write(f"resolved at T1 or T2 (curated Pfam path)\t{shallow}/{len(control_rows)}"
              f"\t{shallow_fraction}\n")
    out.write("NOTE: the controls are Swiss-Prot proteins, so they are named at or before the "
              "swissprot tier and never reach nr (skip_if_named_by). Label accuracy at the nr "
              "tier is therefore NOT validated by the positive controls; the decoys named "
              "per tier below are its only false-positive measurement.\n")
    out.write(f"\ndark set: {n_eligible} target-eligible (no informative name: "
              f"UNCHARACTERIZED_HOMOLOG or NONE), of which {n_dark_strict} have no homologue "
              "at all (NONE; the FESNov definition)\n")
    out.write("\nnegative controls (decoys): expected DARK, reported not halting\n")
    if not decoy_rows:
        out.write("  NONE REACHED THE GATE - the cascade has no false-positive "
                  "measurement for this run\n")
    for construction in sorted(by_class):
        n = by_class[construction]
        named = named_by_class.get(construction, 0)
        out.write(f"  {construction}\t{named}/{n}\t{round(named / n, 4)}\n")
    if decoy_rows:
        out.write("  named per tier (false-positive rate of each tier over all decoys):\n")
        for tier in snakemake.params.tier_sources:
            named = named_by_tier.get(tier, 0)
            out.write(f"    {tier}\t{snakemake.params.tier_sources[tier]}\t{named}/"
                      f"{len(decoy_rows)}\t{round(named / len(decoy_rows), 4)}\n")
    if decoy_named:
        out.write(f"\ndecoys named FUNCTIONAL ({len(decoy_named)}):\n")
        for r in decoy_named[:50]:
            out.write(f"  {r['seq_id']}\t{r.get('annot_tier', '')}\t"
                      f"{r.get('annot_label', '')}\n")
    if failed:
        out.write(f"\ncontrols NOT classed FUNCTIONAL ({len(failed)}):\n")
        for sid in failed[:50]:
            out.write(f"  {sid}\n")
print(summary)

# The gate. Failing here stops the run rather than producing a target list nobody should
# act on.
if recall < cfg["min_control_recall"]:
    raise SystemExit(
        f"S5 QUALITY GATE FAILED: positive control recall {recall} is below the required "
        f"{cfg['min_control_recall']}.\n\n"
        f"{len(failed)} of {len(control_rows)} reviewed Swiss-Prot proteins of KNOWN "
        "function were not classed FUNCTIONAL by this cascade. That means annotation is "
        "losing proteins it should catch, so 'unannotated' carries no information and any "
        "target list built on it would send unannotatable noise to the bench.\n\n"
        f"First failures: {', '.join(failed[:10])}\n\n"
        "NOTE: this gate is about proteins that SHOULD have been annotated. Genuinely "
        "novel proteins - integron cassettes, defence accessories, small ORFs - are "
        "supposed to come out unannotated and are never a reason to halt.\n"
        f"See {snakemake.output.report} for the full list.")
