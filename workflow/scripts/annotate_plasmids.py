"""S4: the primary deliverable - every ORF on every plasmid, with its annotation.

This is a join, not a decision: nothing is classified or filtered here. It expands the
per-unique-protein annotation from S3 back out over every ORF that shares that sequence,
and attaches the S2b artefact flags.

Three things travel with each row that did not exist in v1:

  spans_origin          the ORF was reconstructed across the cut point of a circular
                        plasmid, so it runs start..plasmid_length then 1..end. A naive
                        end - start is negative for these and is a bug.
  artefact_flag         AntiFam or low-complexity evidence that this is not a protein.
                        Flagged, never removed (P5) - the exclusion happens at target
                        selection and stays countable.
  dark_covered_fraction how much of the protein uninformative hits cover, which is the
                        axis that discriminates within the dark set where
                        annot_completeness is constant NONE by construction.

Every row also carries the thresholds that produced its classification, so any downstream
table can be traced back to the numbers that made it (design principle P4).
"""
import _ctx  # noqa: F401
import csv

# Per-unique-protein annotation from the cascade.
with open(snakemake.input.prot, newline="") as fh:
    annot = {r["seq_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

# Artefact flags from S2b, keyed the same way.
with open(snakemake.input.artefact, newline="") as fh:
    artefact = {r["seq_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

# orf_id -> seq_id. One unique protein may correspond to many ORFs; dereplication is a
# compute optimisation and every ORF must reappear here.
orf_to_seq = {}
for line in open(snakemake.input.map):
    sid, ids = line.rstrip("\n").split("\t")
    for oid in ids.split(","):
        orf_to_seq[oid] = sid

cols = [
    # where the ORF is
    "orf_id", "plasmid_id", "start", "end", "strand", "partial", "spans_origin",
    # what it is
    "annot_tier", "annot_label", "functional_class", "homology_depth",
    "annot_qcov", "annot_tcov", "annot_evalue", "n_informative_hits",
    # how much of it is accounted for
    "explained_fraction", "annot_completeness", "meets_min_explained",
    # what the dark evidence says
    "dark_covered_fraction", "dark_completeness", "dark_evidence", "n_dark_databases",
    "uninformative_labels", "uninformative_tiers",
    # is it a protein at all
    "artefact_flag", "antifam_family", "artefact_reason",
    # provenance
    "thr_min_coverage", "thr_min_explained", "thr_narrow_at",
]

n = n_missing = 0
with open(snakemake.input.index, newline="") as fh, \
        open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in csv.DictReader(fh, delimiter="\t"):
        seq_id = orf_to_seq.get(r["orf_id"])
        if seq_id is None:
            # Dereplication is asserted lossless upstream, so this cannot happen; if it
            # ever does, it means the index and the map were built from different runs.
            n_missing += 1
            continue
        a = annot.get(seq_id, {})
        f = artefact.get(seq_id, {})
        # ORF-level columns win over protein-level ones where names collide (start, end).
        row = {k: r.get(k, a.get(k, f.get(k))) for k in cols}
        w.writerow(row)
        n += 1

assert n_missing == 0, f"{n_missing} ORFs had no protein mapping - index and map disagree"
assert n > 0, "no ORFs written - check that S1 and S2 produced output"

with open(snakemake.log[0], "w") as log:
    log.write(f"annotated {n} ORFs across {len(annot)} unique proteins\n")
print(f"annotated {n} ORFs across {len(annot)} unique proteins")
