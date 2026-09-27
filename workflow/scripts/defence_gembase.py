"""Rule defence_gembase, between the two DefenseFinder phases: propagate component hits,
prune, and write genomic order.

Inputs: phase 1 component hits per unique protein, the protein map (seq_id -> orf_ids) and
the ORF index. A hit applies to every ORF sharing the sequence. Plasmids with no component
ORF are left out, since they cannot satisfy any model's quorum; on the test run this kept 41
of 100 plasmids holding 83% of the ORFs, so the saving is about 1.2x. Every ORF of a kept
plasmid is written, not only the component ORFs, because MacSyFinder counts the intervening
genes.

Outputs: the gembase FASTA for phase 2, and the map from gembase_id to orf_id and plasmid.
"""
import csv

import _ctx  # noqa: F401

from plasmidann.defence import candidate_plasmids, gembase_records, propagate_components

with open(snakemake.input.components, newline="") as fh:
    # A NOT_RUN table has one row without a seq_id; it gives no component.
    component_seqs = {r["seq_id"] for r in csv.DictReader(fh, delimiter="\t") if r["seq_id"]}

orf_to_seq = {}
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, ids = line.rstrip("\n").split("\t")
        for oid in ids.split(","):
            orf_to_seq[oid] = sid

component_orfs = propagate_components(component_seqs, orf_to_seq)
keep = candidate_plasmids(component_orfs)

n_written = 0
with open(snakemake.input.index, newline="") as fh, \
        open(snakemake.output.faa, "w") as out, \
        open(snakemake.output.map, "w", newline="") as mapout:
    w = csv.writer(mapout, delimiter="\t")
    w.writerow(["gembase_id", "orf_id", "plasmid_id"])
    for gid, plasmid, orf in gembase_records(csv.DictReader(fh, delimiter="\t"), keep):
        out.write(f">{gid}\n{orf['seq']}\n")
        w.writerow([gid, orf["orf_id"], plasmid])
        n_written += 1

print(f"component ORFs={len(component_orfs)} candidate plasmids={len(keep)} "
      f"proteins written in genomic order={n_written}")
assert n_written or not keep, "candidate plasmids selected but no ORFs written"
