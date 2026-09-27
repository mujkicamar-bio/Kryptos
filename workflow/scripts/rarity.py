"""Stage 14: rarity labels per family, and the dark-family rarefaction curve.

Inputs: recurrence.tsv (the Stage 7 counts), dark_families.tsv, the protein map and
small_plasmids.txt. Outputs:

  family_rarity.tsv            the rarity labels of each family (plasmidann.rarity) with
                               the counts and the lineage thresholds behind them
  dark_family_rarefaction.tsv  dark families discovered against small-plasmid records
                               sampled; every small plasmid is on the axis, with a dark
                               family or without, and a family is discovered on a plasmid
                               through its dark small-plasmid members
"""
import collections
import csv

import _ctx  # noqa: F401

from plasmidann.rarity import rarefaction, rarity_labels, saturation

cfg = snakemake.params.rarity

COLS = ["family_id", "rarity_labels", "independent_plasmid_cluster_count",
        "unique_plasmid_count", "MOB_count", "host_count", "genus_count",
        "n_plasmids_with_species", "rare_max_lineages", "widespread_min_lineages"]

counts = collections.Counter()
n_families = 0
with open(snakemake.input.recurrence, newline="") as fh, \
        open(snakemake.output.rarity, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=COLS, delimiter="\t")
    writer.writeheader()
    for family in csv.DictReader(fh, delimiter="\t"):
        labels = rarity_labels(family, cfg)
        counts.update(labels or ["(none)"])
        n_families += 1
        writer.writerow({
            "family_id": family["family_id"],
            "rarity_labels": ",".join(labels),
            "independent_plasmid_cluster_count":
                family.get("independent_plasmid_cluster_count", ""),
            "unique_plasmid_count": family.get("unique_plasmid_count", ""),
            "MOB_count": family.get("MOB_count", ""),
            "host_count": family.get("host_count", ""),
            "genus_count": family.get("genus_count", ""),
            "n_plasmids_with_species": family.get("n_plasmids_with_species", ""),
            "rare_max_lineages": cfg["rare_max_lineages"],
            "widespread_min_lineages": cfg["widespread_min_lineages"],
        })

print(f"rarity: {n_families} families")
for label, n in counts.most_common():
    print(f"  {label:<26} {n}")

orfs_of_seq = collections.defaultdict(list)
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, orf_ids = line.rstrip("\n").split("\t")
        orfs_of_seq[sid] = orf_ids.split(",")

small_plasmids = {l.strip() for l in open(snakemake.input.small_ids) if l.strip()}
plasmid_families = {p: set() for p in small_plasmids}
with open(snakemake.input.dark_families, newline="") as fh:
    for family in csv.DictReader(fh, delimiter="\t"):
        for member in filter(None, family["small_members"].split(",")):
            for orf_id in orfs_of_seq.get(member, ()):
                plasmid = orf_id.rsplit("|", 1)[0]
                if plasmid in plasmid_families:
                    plasmid_families[plasmid].add(family["family_id"])

curve = rarefaction(plasmid_families, n_replicates=cfg["rarefaction_replicates"],
                    seed=snakemake.params.seed)

with open(snakemake.output.rarefaction, "w", newline="") as out:
    writer = csv.DictWriter(
        out, fieldnames=["n_plasmids", "mean_families", "min_families", "max_families",
                         "n_replicates"], delimiter="\t")
    writer.writeheader()
    writer.writerows(curve)

if curve:
    final = curve[-1]
    print(f"rarefaction: {final['n_plasmids']} small plasmids -> "
          f"{final['mean_families']} dark families "
          f"(replicate range {final['min_families']}-{final['max_families']})")
    gained = saturation(curve)
    if gained == "":
        print("             saturation undefined: fewer than three curve points, or no "
              "gain over the first step")
    else:
        print(f"             final slope is {gained} of the initial slope "
              "(families gained per plasmid added)")
        if gained > 0.2:
            print("             the curve is still climbing: the dark family count is a "
                  "LOWER BOUND, and more plasmids would reveal more families")
        else:
            print("             the curve has flattened: more plasmids of this kind would "
                  "add few new dark families")
else:
    print("rarefaction: no small plasmids - the curve is empty")
