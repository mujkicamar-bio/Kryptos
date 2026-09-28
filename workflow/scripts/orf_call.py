"""Rule orf_call: call ORFs with pyrodigal (meta mode) on every analysis-set record, with
circular-origin repair (darkorf.circular) for closed molecules; writes orfs.tsv."""
import csv
import multiprocessing

import _ctx  # noqa: F401

from darkorf import genecall
from plasmidann.fasta import iter_fasta

min_aa = snakemake.params.min_orf_aa

# min_gene counts the stop codon, so a protein of min_aa residues needs (min_aa + 1) * 3 nt.
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
    writer.writerow(["plasmid_id", "start", "end", "strand", "partial", "partial_begin",
                     "partial_end", "spans_origin", "translation_table", "seq"])
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
                             g["partial"], g["partial_begin"], g["partial_end"],
                             int(g["origin_spanning"]),
                             g["translation_table"], g["seq"]])

# An empty output is never a legitimate result: every plasmid contains genes.
assert n_genes > 0, (
    f"{snakemake.output.tsv}: no genes called from {n_records} records - "
    "check the input is not empty and that pyrodigal is the expected version")

print(f"records={n_records} genes={n_genes} partial={n_partial} "
      f"spans_origin={n_origin} redundant_dropped={n_dropped_dup} "
      # Genes called with table 4 (TGA read as Trp), counted so they are seen.
      f"translation_table_4={n_table4}")
