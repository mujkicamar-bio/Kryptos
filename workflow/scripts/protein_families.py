"""The family table, with what the cascade found in each family.

Inputs: one cluster file per configured resolution (close, intermediate and broad in
config/targets.yaml), written by protein_clustering.py over every unique protein, annotated
or not; protein_annotation.tsv, the dark ids, small_plasmids.txt, the protein map and the
clonal registry (MOB clusters). A family is a sequence cluster and nothing else:
size and distribution are attributes, never filters, and a cluster of one is kept and
labelled ORPHAN.
family_id = <resolution>:<representative>, so an id does not depend on cluster order.

protein_families.tsv has one row per family holding at least one small-plasmid protein;
families of large-plasmid proteins alone are counted in the log. Counts are over all
members:

  family_size, n_orfs           unique proteins, and the gene copies behind them
  dark_member_count             members in the dark set
  annotated_member_count        members named by the cascade or PlasmidScope
                                (functional_class other than NONE, UNCHARACTERIZED_HOMOLOG
                                and NOT_SEARCHED); a family is as bright as its brightest
                                member (Durairaj et al., Nature 2023, 622:646)
  n_unnamed_excluded            searched members named by nothing that the dark set
                                excludes (low-complexity flag, partial-only)
  n_not_searched                members the cascade did not search (cascade_selection),
                                AntiFam artefacts among them
  percentage_dark_in_family, dark_only (a dark member and no named member), family_class
  n_small_members, n_large_members   members on a small / a large plasmid (a protein on
                                both counts in both)
  scope                         small_only or mixed (a member on a large plasmid), each
                                _known or _unknown (a named member or none)
  known_from                    small | large | both | '': where the named members sit

The distribution counts of a family (plasmid records, lineages, hosts, MOB clusters,
habitats) are computed in one place, recurrence.tsv (rule recurrence).

Within a family holding an unexplained small-plasmid protein, every member not annotated by
Tier 0 or flagged by AntiFam was searched by the cascade, so at the primary resolution known
and unknown mean the same on both sides.

dark_families.tsv is the subset at the primary resolution with a dark small-plasmid member,
listing only the dark members (and the dark small-plasmid members apart).
"""
import collections
import csv
import pathlib

import _ctx  # noqa: F401

from darkorf import ids

cfg = snakemake.params.clustering
# Classes that name nothing, or were never looked at. Everything else is a named member.
UNNAMED = {"NONE", "UNCHARACTERIZED_HOMOLOG", "NOT_SEARCHED"}

# The dark set is read, not recomputed, so the dark fraction uses the same definition.
dark = {line.strip() for line in open(snakemake.input.dark_ids) if line.strip()}

functional_class = {}
with open(snakemake.input.prot, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        functional_class[row["seq_id"]] = row["functional_class"]

small_plasmids = {line.strip() for line in open(snakemake.input.small_ids) if line.strip()}

seq_to_plasmids = collections.defaultdict(set)
seq_to_orf_count = collections.Counter()
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, orf_ids = line.rstrip("\n").split("\t")
        members = orf_ids.split(",")
        seq_to_orf_count[sid] = len(members)
        for orf_id in members:
            seq_to_plasmids[sid].add(orf_id.rsplit("|", 1)[0])
on_small = {sid for sid, ps in seq_to_plasmids.items() if ps & small_plasmids}
on_large = {sid for sid, ps in seq_to_plasmids.items() if ps - small_plasmids}

with open(snakemake.input.registry, newline="") as fh:
    mob_of = {r["plasmid_id"]: r["mob_cluster"] for r in csv.DictReader(fh, delimiter="\t")}


COLS = [
    "family_id", "family_resolution", "representative",
    "family_size", "n_orfs", "dark_member_count", "annotated_member_count",
    "n_unnamed_excluded",
    "percentage_dark_in_family", "family_class", "dark_only",
    "n_small_members", "n_large_members", "n_not_searched", "scope", "known_from",
    "members",
]

# dark_families.tsv holds families that contain a dark member, not only dark-only ones:
# a partly dark family still holds dark proteins to measure, and its named members say
# what the family does.
primary = cfg["primary"]
dark_rows = []

# protein_clustering writes one families_<resolution>_cluster.tsv per resolution.
cluster_file = {pathlib.Path(c).name: c for c in snakemake.input.clusters}
assert len(cluster_file) == len(cfg["resolutions"]), (cluster_file, cfg["resolutions"])

summary = []
with open(snakemake.output.families, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=COLS, delimiter="\t")
    writer.writeheader()

    for resolution in sorted(cfg["resolutions"]):
        clusters = cluster_file[f"families_{resolution}_cluster.tsv"]
        members = collections.defaultdict(list)
        with open(clusters) as fh:
            for line in fh:
                rep, mem = line.rstrip("\n").split("\t")
                members[rep].append(mem)

        n_dark_only = n_orphan = n_large_only = 0
        for rep, mem in sorted(members.items()):
            small = [m for m in mem if m in on_small]
            if not small:
                n_large_only += 1
                continue
            named = [m for m in mem if functional_class.get(m, "NOT_SEARCHED") not in UNNAMED]
            # A named protein found on small AND large plasmids is known from both sides.
            known_from = ({"small" for m in named if m in on_small}
                          | {"large" for m in named if m in on_large})
            n_not_searched = sum(
                1 for m in mem if functional_class.get(m, "NOT_SEARCHED") == "NOT_SEARCHED")
            plasmids = set()
            n_orfs = 0
            for m in mem:
                plasmids |= seq_to_plasmids.get(m, set())
                n_orfs += seq_to_orf_count.get(m, 0)

            mobs = {mob_of[p] for p in plasmids if mob_of.get(p)}

            n_dark = sum(1 for m in mem if m in dark)
            # A searched, unnamed protein the dark set excludes (low-complexity flag,
            # partial-only) is neither dark nor annotated, and is counted apart.
            n_annotated = len(named)
            n_unnamed_excluded = len(mem) - n_dark - n_annotated - n_not_searched
            pct_dark = round(100.0 * n_dark / len(mem), 2)

            # Compared on the counts, not the rounded percentage: 99.996% rounds to 100.0.
            # A family of excluded or unsearched members alone has no dark member and is
            # not dark-only.
            dark_only = int(n_annotated == 0 and n_dark > 0)
            n_dark_only += dark_only
            family_class = "ORPHAN" if len(mem) == 1 else "FAMILY"
            scope = (("mixed" if any(m in on_large for m in mem) else "small_only")
                     + ("_known" if named else "_unknown"))
            known_from = ("both" if len(known_from) == 2 else
                          known_from.pop() if known_from else "")
            n_orphan += family_class == "ORPHAN"
            n_large = sum(1 for m in mem if m in on_large)

            writer.writerow({
                "family_id": ids.family_id(resolution, rep),
                "family_resolution": resolution,
                "representative": rep,
                "family_size": len(mem),
                "n_orfs": n_orfs,
                "dark_member_count": n_dark,
                "annotated_member_count": n_annotated,
                "n_unnamed_excluded": n_unnamed_excluded,
                "percentage_dark_in_family": pct_dark,
                "family_class": family_class,
                "dark_only": dark_only,
                "n_small_members": len(small),
                "n_large_members": n_large,
                "n_not_searched": n_not_searched,
                "scope": scope,
                "known_from": known_from,
                "members": ",".join(mem),
            })

            small_dark = [m for m in small if m in dark]
            # A dark family is one with a dark SMALL-plasmid member: the dark proteins of
            # large plasmids matter here only as relatives of those.
            if resolution == primary and small_dark:
                dark_members = [m for m in mem if m in dark]
                dark_rows.append({
                    "family_id": ids.family_id(resolution, rep),
                    # The representative must be a dark small-plasmid member: rule
                    # structure_search searches this sequence structurally, and searching
                    # any other member would spend the ProstT5 budget on a protein outside
                    # the study.
                    "representative": rep if rep in small_dark else small_dark[0],
                    "n_members": len(mem),
                    "n_orfs": n_orfs,
                    "n_plasmids": len(plasmids),
                    "n_mob_clusters": len(mobs),
                    "family_class": family_class,
                    "dark_member_count": n_dark,
                    "annotated_member_count": n_annotated,
                    "percentage_dark_in_family": pct_dark,
                    "dark_only": dark_only,
                    "n_small_members": len(small),
                    "n_large_members": n_large,
                    "scope": scope,
                    "known_from": known_from,
                    # Only the dark members: the downstream stages measure the DARK
                    # proteins' evolution and context, and including annotated members
                    # here would silently widen every one of those measurements. The
                    # stages report small-plasmid members alone and all members side by
                    # side, so both lists are kept.
                    "members": ",".join(dark_members),
                    "small_members": ",".join(small_dark),
                })

        summary.append((resolution, len(members) - n_large_only, n_orphan,
                        n_dark_only, n_large_only))

DARK_COLS = ["family_id", "representative", "n_members", "n_orfs", "n_plasmids",
             "n_mob_clusters", "family_class", "dark_member_count",
             "annotated_member_count", "percentage_dark_in_family", "dark_only",
             "n_small_members", "n_large_members", "scope", "known_from", "members",
             "small_members"]
with open(snakemake.output.dark_families, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=DARK_COLS, delimiter="\t")
    writer.writeheader()
    writer.writerows(dark_rows)

print(f"{'resolution':<14}{'families':>10}{'orphans':>10}{'dark_only':>11}"
      f"{'large_only_not_written':>24}")
for resolution, n_fam, n_orphan, n_dark_only, n_large in summary:
    print(f"{resolution:<14}{n_fam:>10}{n_orphan:>10}{n_dark_only:>11}{n_large:>24}")

print(f"dark families at {primary}: {len(dark_rows)} with a dark small-plasmid member, "
      f"of which {sum(r['dark_only'] for r in dark_rows)} are 100% dark; by scope: "
      + " ".join(f"{k}={v}" for k, v in sorted(
          collections.Counter(r["scope"] for r in dark_rows).items())))
