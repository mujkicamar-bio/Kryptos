"""Stage 14: rarity labels per family, and the dark-family rarefaction curve (55).

Two outputs, because they answer two questions:

  family_rarity.tsv            what kind of distribution does each family have
  dark_family_rarefaction.tsv  has the collection saturated, or would more plasmids keep
                               revealing new dark families

See src/plasmidann/rarity.py for why breadth is counted in independent lineages rather than
plasmid records, and why the curve is averaged over replicate orderings.
"""
import collections
import csv

import _ctx  # noqa: F401

from plasmidann.rarity import RARITY_VERSION, rarefaction, rarity_labels, saturation

cfg = snakemake.params.rarity

# ------------------------------------------------------------------------------------
# Rarity labels, from the Stage 7 distribution counts.
# ------------------------------------------------------------------------------------
COLS = ["family_id", "rarity_labels", "independent_plasmid_cluster_count",
        "unique_plasmid_count", "MOB_count", "host_count", "genus_count",
        "n_plasmids_with_species", "rarity_version", "rare_max_lineages", "widely_conserved_min_lineages"]

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
            # The thresholds travel with every row: a label is meaningless without the
            # number that produced it, and section 54 requires them to be explicit.
            "rarity_version": RARITY_VERSION,
            "rare_max_lineages": cfg["rare_max_lineages"],
            "widely_conserved_min_lineages": cfg["widely_conserved_min_lineages"],
        })

print(f"rarity: {n_families} families")
for label, n in counts.most_common():
    print(f"  {label:<26} {n}")

# ------------------------------------------------------------------------------------
# Rarefaction: dark families discovered against SMALL plasmids sampled (section 55).
#
# The x-axis is every small plasmid, including the ones that carry no dark family: those
# are part of what was sampled, and leaving them out (58 of 100 plasmids were kept on the
# test set) made each step look richer than a real sample of small plasmids is. A family
# is discovered on a plasmid through its dark small-plasmid members (small_members), since
# the large plasmids are not on the axis.
# ------------------------------------------------------------------------------------
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
    print(f"             final slope is {gained} of the initial slope "
          "(families gained per plasmid added)")
    # Descriptive, never a target. Section 55: the ~1,000-candidate figure is an
    # experimental-budget objective, not a biological assumption.
    if gained != "" and float(gained) > 0.2:
        print("             the curve is still climbing: the dark family count is a "
              "LOWER BOUND, and more plasmids would reveal more families")
    else:
        print("             the curve has flattened: more plasmids of this kind would "
              "add few new dark families")
else:
    print("rarefaction: no dark families on any plasmid - the curve is empty")
