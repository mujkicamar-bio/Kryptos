"""S1: call every ORF on every plasmid, with circular-origin repair.

TWO DEFECTS FROM v1 ARE ADDRESSED HERE.

1. The 0-byte FASTA. v1 declared an output `faa`, opened it, and never wrote a single
   byte to it - all 600 files were empty, and because `dereplicate` declared them as
   inputs but only ever read `orf_index.tsv`, the DAG reported success. The fix is to
   REMOVE the output rather than populate it: the protein sequence already travels in the
   `seq` column of orf_index.tsv, and the FASTA that downstream stages actually consume is
   results/03_dereplication/unique_proteins.faa. Deleting the unused output removes the whole class of
   failure instead of the symptom.

2. Genes broken by linearising a circle. See src/darkorf/circular.py for the full
   explanation; in short, 93.9% of these plasmids are closed molecules written as lines,
   and 160,375 ORFs (1.12 per plasmid) are fragments created by the cut rather than by
   biology. Those fragments are indistinguishable from novel dark proteins downstream.
"""
import csv
import multiprocessing

import _ctx  # noqa: F401

from darkorf import genecall
from plasmidann.fasta import iter_fasta

min_aa = snakemake.params.min_orf_aa

# Pyrodigal's min_gene counts the STOP codon, so a protein of exactly `min_aa` residues
# needs (min_aa + 1) * 3 nucleotides. v1 passed min_aa * 3, which set the real floor one
# residue lower than declared - the reviewer's "min_orf_aa=20 is off by one, the floor is
# 19 aa". The correction changes the ORF count by ~0.14%; it is made because a declared
# threshold that does not mean what it says cannot be reasoned about, not because the
# ORFs matter numerically.
min_gene_nt = (min_aa + 1) * 3

topology = {}
with open(snakemake.input.master, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        topology[row["plasmid_id"]] = row.get("topology", "")

n_records = n_genes = n_partial = n_origin = n_dropped_dup = n_table4 = 0

records = ((rec, topology.get(rec[0], "")) for rec in iter_fasta([snakemake.input.fasta]))
with open(snakemake.output.tsv, "w", newline="") as tsv, \
        multiprocessing.Pool(snakemake.threads, initializer=genecall.configure,
                             initargs=(min_gene_nt,)) as pool:
    writer = csv.writer(tsv, delimiter="\t")
    writer.writerow(["plasmid_id", "start", "end", "strand", "partial",
                     "spans_origin", "translation_table", "seq"])
    # imap keeps input order, so the output is deterministic whatever the pool size.
    for plasmid_id, genes, dropped in pool.imap(genecall.call_record, records, chunksize=64):
        n_records += 1
        n_dropped_dup += dropped
        for g in genes:
            n_genes += 1
            n_partial += g["partial"]
            n_origin += g["origin_spanning"]
            n_table4 += g["translation_table"] == 4
            writer.writerow([plasmid_id, g["start"], g["end"], g["strand"],
                             g["partial"], int(g["origin_spanning"]),
                             g["translation_table"], g["seq"]])

# An empty output here is always a bug, never a legitimate result: every plasmid contains
# genes. v1's silent 0-byte outputs are the reason this assertion exists.
assert n_genes > 0, (
    f"{snakemake.output.tsv}: no genes called from {n_records} records - "
    "check the input is not empty and that pyrodigal is the expected version")

print(f"records={n_records} genes={n_genes} partial={n_partial} "
      f"spans_origin={n_origin} redundant_dropped={n_dropped_dup} "
      # Genes called with table 4 (TGA read as Trp): an open issue, counted so it is seen.
      f"translation_table_4={n_table4}")
