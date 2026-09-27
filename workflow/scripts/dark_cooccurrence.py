"""Dark protein sequences that travel together, per Mash lineage (plasmidann.cooccurrence).

Inputs: dark_families.tsv (its `members` are the dark sequences of the dark families),
the protein map and plasmid_lineage.tsv. Output dark_cooccurrence.tsv, one row per pair of
sequences together in >= min_lineages_together lineages (seq_a < seq_b), sorted by p-value.
"""
import collections
import csv

import _ctx  # noqa: F401

from plasmidann.cooccurrence import cooccurrence

cfg = snakemake.params.cooccurrence

dark = set()
with open(snakemake.input.families, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        dark.update(filter(None, r["members"].split(",")))

plasmid_seqs = collections.defaultdict(set)
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, orf_ids = line.rstrip("\n").split("\t")
        if sid in dark:
            for orf in orf_ids.split(","):
                plasmid_seqs[orf.rsplit("|", 1)[0]].add(sid)

with open(snakemake.input.lineage, newline="") as fh:
    lineage_of = {r["plasmid_id"]: r["plasmid_lineage_cluster"]
                  for r in csv.DictReader(fh, delimiter="\t")}
n_lineages = len(set(lineage_of.values()))

rows, stats = cooccurrence(plasmid_seqs, lineage_of, n_lineages,
                           cfg["min_lineages_together"])

COLS = ["seq_a", "seq_b", "n_lineages_a", "n_lineages_b", "n_lineages_together",
        "n_lineages_total", "fraction_of_a", "fraction_of_b", "expected_together",
        "p_value", "q_value"]
with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=COLS, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"dark sequences {len(dark)}, {stats['n_sequences']} on {len(plasmid_seqs)} plasmids "
      f"in {n_lineages} lineages; enumerated (in >= {cfg['min_lineages_together']} "
      f"lineages): {stats['n_enumerated']}; pair occurrences on plasmids: "
      f"{stats['n_pair_occurrences']}; pairs tested: {stats['n_tested']}; reported: "
      f"{len(rows)}; q <= {cfg['fdr']}: {sum(r['q_value'] <= cfg['fdr'] for r in rows)}")
