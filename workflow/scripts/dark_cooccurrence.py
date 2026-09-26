"""S8g: dark families that travel together, per lineage.

Together = both families have a member ORF on the same plasmid; unit = the Stage 6 lineage;
family = the dark family at the primary resolution, through its dark members. The test,
its minimum and the multiple-testing correction are in plasmidann.cooccurrence.

Output dark_cooccurrence.tsv, one row per TESTED pair (family_a < family_b), sorted by
p-value. Pairs together in fewer than min_lineages_together lineages are not written: they
were not tested, and every other pair of families is in no lineage together at all.
"""
import collections
import csv

import _ctx  # noqa: F401

from plasmidann.cooccurrence import cooccurrence

cfg = snakemake.params.cooccurrence

family_of_seq = {}
with open(snakemake.input.families, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        for member in r["members"].split(","):
            family_of_seq[member] = r["family_id"]

plasmid_families = collections.defaultdict(set)
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, orf_ids = line.rstrip("\n").split("\t")
        fid = family_of_seq.get(sid)
        if fid is not None:
            for orf in orf_ids.split(","):
                plasmid_families[orf.rsplit("|", 1)[0]].add(fid)

with open(snakemake.input.lineage, newline="") as fh:
    lineage_of = {r["plasmid_id"]: r["plasmid_lineage_cluster"]
                  for r in csv.DictReader(fh, delimiter="\t")}
n_lineages = len(set(lineage_of.values()))

rows = cooccurrence(plasmid_families, lineage_of, n_lineages, cfg["min_lineages_together"])

COLS = ["family_a", "family_b", "n_lineages_a", "n_lineages_b", "n_lineages_together",
        "n_lineages_total", "fraction_of_a", "fraction_of_b", "expected_together",
        "p_value", "q_value", "status"]
with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=COLS, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

per_plasmid = collections.Counter(len(f) for f in plasmid_families.values())
n_pairs_on_plasmids = sum(n * m * (m - 1) // 2 for m, n in per_plasmid.items())
print(f"dark families {len(set(family_of_seq.values()))} on {len(plasmid_families)} "
      f"plasmids in {n_lineages} lineages; family pairs sharing a plasmid, counted per "
      f"plasmid: {n_pairs_on_plasmids}; tested (together in >= "
      f"{cfg['min_lineages_together']} lineages): {len(rows)}; "
      f"q <= {cfg['fdr']}: {sum(r['q_value'] <= cfg['fdr'] for r in rows)}")
