"""Stage 14: rarity labels per family, and the dark-family rarefaction curve.

Inputs: recurrence.tsv (the Stage 7 counts), dark_families.tsv, the protein map,
small_plasmids.txt and the Stage 6 lineages. Outputs:

  family_rarity.tsv            the rarity labels of each family (plasmidann.rarity) with
                               the counts and the lineage thresholds behind them; the
                               WIDESPREAD threshold is computed per resolution, from the
                               lineage counts of that resolution's families
  dark_family_rarefaction.tsv  dark families discovered against lineages sampled; every
                               lineage holding a small plasmid is on the axis, with a dark
                               family or without, and a family is discovered in a lineage
                               through its dark members on the lineage's small plasmids
"""
import collections
import csv

import _ctx  # noqa: F401

from plasmidann.rarity import rarefaction, rarity_labels, saturation, widespread_threshold

cfg = snakemake.params.rarity

COLS = ["family_id", "rarity_labels", "independent_plasmid_cluster_count",
        "unique_plasmid_count", "MOB_count", "host_count", "genus_count",
        "n_plasmids_with_species", "rare_max_lineages", "widespread_percentile",
        "widespread_min_lineages"]

# First pass: the lineage counts per resolution, from which the WIDESPREAD threshold is
# measured. Families of different resolutions are different units and are not pooled.
lineage_counts = collections.defaultdict(list)
with open(snakemake.input.recurrence, newline="") as fh:
    for family in csv.DictReader(fh, delimiter="\t"):
        lineage_counts[family["family_resolution"]].append(
            int(family["independent_plasmid_cluster_count"] or 0))
widespread = {res: widespread_threshold(c, cfg["widespread_percentile"])
              for res, c in lineage_counts.items()}
for res, t in sorted(widespread.items()):
    n = sum(1 for c in lineage_counts[res] if c > 0)
    print(f"WIDESPREAD at {res}: >= {t} lineages, the {cfg['widespread_percentile']}th "
          f"percentile of {n} families with a measured lineage count")
    if t is not None and t <= cfg["rare_max_lineages"]:
        print(f"  the threshold is at most rare_max_lineages ({cfg['rare_max_lineages']}), "
              "so a family at it is both RARE and WIDESPREAD")

counts = collections.Counter()
n_families = 0
with open(snakemake.input.recurrence, newline="") as fh, \
        open(snakemake.output.rarity, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=COLS, delimiter="\t")
    writer.writeheader()
    for family in csv.DictReader(fh, delimiter="\t"):
        threshold = widespread[family["family_resolution"]]
        labels = rarity_labels(family, {**cfg, "widespread_min_lineages": threshold})
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
            "widespread_percentile": cfg["widespread_percentile"],
            "widespread_min_lineages": "" if threshold is None else threshold,
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
with open(snakemake.input.lineage, newline="") as fh:
    lineage_of = {r["plasmid_id"]: r["plasmid_lineage_cluster"]
                  for r in csv.DictReader(fh, delimiter="\t")}
missing = sorted(small_plasmids - lineage_of.keys())
if missing:
    raise SystemExit(f"rarity: {len(missing)} small plasmid(s) have no lineage, "
                     f"e.g. {missing[:3]}")
lineage_families = {lineage_of[p]: set() for p in small_plasmids}
with open(snakemake.input.dark_families, newline="") as fh:
    for family in csv.DictReader(fh, delimiter="\t"):
        for member in filter(None, family["small_members"].split(",")):
            for orf_id in orfs_of_seq.get(member, ()):
                plasmid = orf_id.rsplit("|", 1)[0]
                if plasmid in small_plasmids:
                    lineage_families[lineage_of[plasmid]].add(family["family_id"])

curve = rarefaction(lineage_families, n_replicates=cfg["rarefaction_replicates"],
                    seed=snakemake.params.seed)

with open(snakemake.output.rarefaction, "w", newline="") as out:
    writer = csv.DictWriter(
        out, fieldnames=["n_lineages", "mean_families", "min_families", "max_families",
                         "n_replicates"], delimiter="\t")
    writer.writeheader()
    writer.writerows(curve)

if curve:
    final = curve[-1]
    print(f"rarefaction: {final['n_lineages']} lineages with a small plasmid -> "
          f"{final['mean_families']} dark families "
          f"(replicate range {final['min_families']}-{final['max_families']})")
    gained = saturation(curve)
    if gained == "":
        print("             saturation undefined: fewer than three curve points, or no "
              "gain over the first step")
    else:
        print(f"             final slope is {gained} of the initial slope "
              "(families gained per lineage added)")
        if gained > 0.2:
            print("             the curve is still climbing: the dark family count is a "
                  "LOWER BOUND, and more lineages would reveal more families")
        else:
            print("             the curve has flattened: more lineages of this kind would "
                  "add few new dark families")
else:
    print("rarefaction: no small plasmids - the curve is empty")
