"""S5: which proteins are screening candidates at all.

A protein is target-eligible when the cascade gave it no informative name
(UNCHARACTERIZED_HOMOLOG or NONE), it was searched, and it is not artefact-flagged. No
protein is excluded by family name, and none is deleted: every protein gets a row with its
exclusion reasons. The report counts the eligible proteins and, of those, the ones with no
homologue at all.
"""
import csv

import _ctx  # noqa: F401

with open(snakemake.input.prot, newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

artefact_ids = set()
with open(snakemake.input.artefact, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["artefact_flag"] == "1":
            artefact_ids.add(r["seq_id"])

cols = ["seq_id", "is_artefact", "target_eligible", "exclusion_reason"]
n_eligible = 0
# The dark set under two definitions. Ours counts a protein whose only
# homologues are themselves unnamed (UNCHARACTERIZED_HOMOLOG) as dark; FESNov
# (Rodriguez del Rio et al. 2024, Nature 626:377) calls a family unknown only when it has
# no homologue at all, which here is functional_class NONE. Both counts are reported.
n_dark_strict = 0
with open(snakemake.output.flags, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in rows:
        sid = r["seq_id"]
        reasons = []
        if sid in artefact_ids:
            reasons.append("artefact")
        # Only proteins nothing could name are screening candidates at all. A protein the
        # cascade selection did not search is neither named nor dark, and says so.
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
           f"target_eligible={n_eligible} dark_no_homologue={n_dark_strict}")
with open(snakemake.output.report, "w") as out:
    out.write(summary + "\n")
    out.write(f"\ndark set: {n_eligible} target-eligible (no informative name: "
              f"UNCHARACTERIZED_HOMOLOG or NONE), of which {n_dark_strict} have no homologue "
              "at all (NONE; the FESNov definition)\n")
print(summary)
