"""Every ORF on every plasmid, with the annotation of its unique protein.

A join, not a decision: nothing is classified or filtered here. The per-protein cascade
annotation and the ORF QC artefact flags are expanded over every ORF that shares the
protein sequence; the ORF index supplies position, strand, partial and spans_origin (an ORF
reconstructed across the origin runs start..plasmid length, then 1..end). Artefact-flagged
proteins are kept and flagged. Each row carries the coverage, explained-fraction and
narrow-hit thresholds of its classification.
"""
import csv

import _ctx  # noqa: F401

# Per-unique-protein annotation from the cascade.
with open(snakemake.input.prot, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    annot = {r["seq_id"]: r for r in reader}
    known = set(reader.fieldnames)

# Artefact flags from ORF QC, keyed the same way.
with open(snakemake.input.artefact, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    artefact = {r["seq_id"]: r for r in reader}
    known |= set(reader.fieldnames)

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
    # 11, or 4 where pyrodigal's meta mode chose the Mycoplasma code
    "translation_table",
    # what it is, and where that came from (self, representative, plasmidscope,
    # not_searched, artefact_antifam; see cascade_resolve.py)
    "annot_source", "annot_representative",
    "annot_tier", "annot_label", "functional_class", "homology_depth",
    "annot_qcov", "annot_tcov", "annot_evalue", "n_informative_hits",
    # how much of it is accounted for
    "explained_fraction", "annot_completeness",
    # what the dark evidence says
    "dark_covered_fraction", "dark_completeness", "dark_evidence", "n_dark_databases",
    "uninformative_labels", "uninformative_tiers",
    # is it a protein at all
    "artefact_flag", "antifam_family", "artefact_reason",
    # provenance
    "thr_min_coverage", "thr_narrow_at",
]

n = n_missing = 0
with open(snakemake.input.index, newline="") as fh, \
        open(snakemake.output[0], "w", newline="") as out:
    reader = csv.DictReader(fh, delimiter="\t")
    # A column no input carries would otherwise be written empty on every row.
    missing = [c for c in cols if c not in known | set(reader.fieldnames)]
    if missing:
        raise SystemExit(f"annotate_plasmids: no input carries the columns {missing}")
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in reader:
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
assert n > 0, "no ORFs written - check that ORF calling and dereplication produced output"

with open(snakemake.log[0], "w") as log:
    log.write(f"annotated {n} ORFs across {len(annot)} unique proteins\n")
print(f"annotated {n} ORFs across {len(annot)} unique proteins")
