"""S8a phase 1.5: propagate component labels, prune, and write genomic order.

Three things happen here, and the middle one is the reason the split works at all.

PROPAGATE   A component hit on a unique protein applies to every ORF that shares the
            sequence. This is what makes it legitimate to have searched only the
            dereplicated set.

PRUNE       A plasmid carrying no component cannot satisfy any model's quorum, so it never
            needs the expensive ordered representation. Measured on the test run it kept
            41 of 100 plasmids but 83% of their ORFs - the plasmids it removes are the
            small ones - so the saving in what reaches MacSyFinder is ~1.2x, not the 4-5x
            once estimated.

ORDER       ORFs are written in genomic order under MacSyFinder's gembase naming, so one
            file can hold many replicons and each is still treated separately. Origin-
            spanning genes sort to the front, because on a circular replicon the gene
            straddling the cut precedes position 1 - and 94% of these plasmids are circular.
"""
import collections
import csv

import _ctx  # noqa: F401

from plasmidann.defence import candidate_plasmids, gembase_id, order_orfs, propagate_components

# component hits on unique proteins, from phase 1
component_of_seq = {}
with open(snakemake.input.components, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        component_of_seq[r["seq_id"]] = r["component"]

# orf_id -> seq_id
orf_to_seq = {}
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, ids = line.rstrip("\n").split("\t")
        for oid in ids.split(","):
            orf_to_seq[oid] = sid

orf_labels = propagate_components(component_of_seq, orf_to_seq)
keep = candidate_plasmids(orf_labels)

# Collect every ORF of every candidate plasmid - not only the labelled ones. MacSyFinder
# counts INTERVENING genes, so the unlabelled neighbours are what inter_gene_max_space is
# measured across. Omitting them would make every system look artificially compact.
by_plasmid = collections.defaultdict(list)
with open(snakemake.input.index, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["plasmid_id"] in keep:
            by_plasmid[r["plasmid_id"]].append(r)

n_written = 0
with open(snakemake.output.faa, "w") as out, \
        open(snakemake.output.map, "w", newline="") as mapout:
    w = csv.writer(mapout, delimiter="\t")
    w.writerow(["gembase_id", "orf_id", "plasmid_id", "position", "component"])
    for plasmid in sorted(by_plasmid):
        for position, orf in enumerate(order_orfs(by_plasmid[plasmid]), start=1):
            gid = gembase_id(plasmid, position)
            out.write(f">{gid}\n{orf['seq']}\n")
            w.writerow([gid, orf["orf_id"], plasmid, position,
                        orf_labels.get(orf["orf_id"], "")])
            n_written += 1

print(f"labelled ORFs={len(orf_labels)} candidate plasmids={len(keep)} "
      f"proteins written in genomic order={n_written}")
assert n_written or not keep, "candidate plasmids selected but no ORFs written"
