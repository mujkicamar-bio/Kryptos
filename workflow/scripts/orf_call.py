"""S1: call every ORF on every plasmid, with circular-origin repair.

TWO DEFECTS FROM v1 ARE ADDRESSED HERE.

1. The 0-byte FASTA. v1 declared an output `faa`, opened it, and never wrote a single
   byte to it - all 600 files were empty, and because `dereplicate` declared them as
   inputs but only ever read `orf_index.tsv`, the DAG reported success. The fix is to
   REMOVE the output rather than populate it: the protein sequence already travels in the
   `seq` column of orf_index.tsv, and the FASTA that downstream stages actually consume is
   results/s2/unique_proteins.faa. Deleting the unused output removes the whole class of
   failure instead of the symptom.

2. Genes broken by linearising a circle. See src/darkorf/circular.py for the full
   explanation; in short, 93.9% of these plasmids are closed molecules written as lines,
   and 160,375 ORFs (1.12 per plasmid) are fragments created by the cut rather than by
   biology. Those fragments are indistinguishable from novel dark proteins downstream.
"""
import _ctx  # noqa: F401
import csv

import pyrodigal
from darkorf.circular import is_circular, overlap_for, resolve_origin_genes

min_aa = snakemake.params.min_orf_aa

# Pyrodigal's min_gene counts the STOP codon, so a protein of exactly `min_aa` residues
# needs (min_aa + 1) * 3 nucleotides. v1 passed min_aa * 3, which set the real floor one
# residue lower than declared - the reviewer's "min_orf_aa=20 is off by one, the floor is
# 19 aa". The correction changes the ORF count by ~0.14%; it is made because a declared
# threshold that does not mean what it says cannot be reasoned about, not because the
# ORFs matter numerically.
min_gene_nt = (min_aa + 1) * 3

gene_finder = pyrodigal.GeneFinder(
    meta=True,                              # no training set: every plasmid is called alone
    min_gene=min_gene_nt,
    min_edge_gene=min(min_gene_nt, 60),     # edge genes on genuinely linear molecules
)

# Topology decides whether a record gets origin repair. Unknown topology is treated as
# linear: extending a genuinely linear molecule would fabricate a junction that does not
# exist and could invent a chimeric gene across the two ends.
topology = {}
with open(snakemake.input.master, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        topology[row["plasmid_id"]] = row.get("topology", "")

n_records = n_genes = n_partial = n_origin = n_dropped_dup = 0

with open(snakemake.output.tsv, "w", newline="") as tsv:
    writer = csv.writer(tsv, delimiter="\t")
    writer.writerow(["plasmid_id", "start", "end", "strand", "partial",
                     "spans_origin", "seq"])

    def call(plasmid_id, chunks):
        """Call genes for one record and write them out."""
        global n_records, n_genes, n_partial, n_origin, n_dropped_dup
        if plasmid_id is None:
            return
        sequence = "".join(chunks)
        length = len(sequence)
        if not length:
            return
        n_records += 1

        circular = is_circular(topology.get(plasmid_id))
        if circular:
            # Append the head to the tail so a gene straddling the cut becomes contiguous.
            search_seq = sequence + sequence[:overlap_for(length)]
        else:
            search_seq = sequence

        raw = [{"start": g.begin, "end": g.end, "strand": g.strand,
                "partial": int(g.partial_begin or g.partial_end),
                "seq": g.translate().rstrip("*")}
               for g in gene_finder.find_genes(search_seq)]

        if circular:
            before = len(raw)
            genes = resolve_origin_genes(raw, original_length=length)
            n_dropped_dup += before - len(genes)
        else:
            genes = [dict(g, origin_spanning=False) for g in raw]

        for g in genes:
            n_genes += 1
            n_partial += g["partial"]
            n_origin += g["origin_spanning"]
            writer.writerow([plasmid_id, g["start"], g["end"], g["strand"],
                             g["partial"], int(g["origin_spanning"]), g["seq"]])

    current_id, buffer = None, []
    for line in open(snakemake.input.fasta):
        if line[0] == ">":
            call(current_id, buffer)
            current_id, buffer = line[1:].split()[0], []
        else:
            buffer.append(line.strip())
    call(current_id, buffer)

# An empty output here is always a bug, never a legitimate result: every shard contains
# plasmids and every plasmid contains genes. v1's silent 0-byte outputs are the reason
# this assertion exists.
assert n_genes > 0, (
    f"{snakemake.output.tsv}: no genes called from {n_records} records - "
    "check the shard is not empty and that pyrodigal is the expected version")

print(f"records={n_records} genes={n_genes} partial={n_partial} "
      f"spans_origin={n_origin} redundant_dropped={n_dropped_dup}")
