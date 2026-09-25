"""Stage 5: the family table - EVERY unique protein's families, with what the cascade found.

The clustering itself runs before the cascade (S2f, protein_clustering.py), because the
selection of what the cascade searches is made on the families. This stage reads those
clusters and adds the annotation, the distribution and the small/large scope.

WHY EVERY PROTEIN AND NOT ONLY THE DARK ONES

An earlier version clustered the dark set alone. That cannot produce the table the
specification asks for. Section 31.2 requires

    family_size, dark_member_count, annotated_member_count, percentage_dark_in_family

and section 32 DERIVES a dark-only family as percentage_dark_in_family = 100%. All four
need the annotated members to be present in the clustering. With a dark-only input every
family is trivially 100% dark, the statistic carries no information, and a dark protein
sitting in a family of well-annotated homologs - which is a strong and interesting
observation - is indistinguishable from one that is genuinely alone.

Clustering everything also gives the project what it actually wants from a family table:
you can look up any protein, annotated or not, and see its relatives.

A FAMILY IS A SEQUENCE CLUSTER. NOTHING ELSE.

Section 31: "family membership is a descriptive sequence relationship. It does not
automatically establish function." An earlier design defined a family as ">= 3 members from
>= 2 distinct MOB clusters", which made the family COUNT a function of metadata
completeness rather than of sequence, and silently deleted lineage-restricted families when
a novel system confined to one plasmid lineage may be exactly what is worth finding. Size
and distribution are ATTRIBUTES computed here, never filters. Clusters of one are retained
and labelled ORPHAN.

THREE RESOLUTIONS (section 31.4)

close, intermediate and broad, with thresholds declared in config/targets.yaml rather than
implied by a default. One protein therefore belongs to three families, and the table says
which resolution each row is for. `primary` names the one downstream stages read.

FAMILY IDS ARE CONTENT-DERIVED (section 5.4)

    family_id = <resolution>:<representative_protein_id>

NOT an ordinal. The previous version numbered families F0000001, F0000002, ... in cluster
order, so inserting one protein anywhere in the collection renumbered every family after
it, and no family id could be compared between two runs. Section 31 states the requirement
directly: "family IDs must not depend on result ordering."

SMALL PLASMIDS AND LARGE ONES

The study is about small plasmids (01_analysis_set/small_plasmids.txt). Every plasmid is
clustered, and a family is written only when it holds at least one protein that occurs on
a small plasmid; the rest are counted in the log. All counts are over ALL members, so no
member is lost; the split is in its own columns:

    n_small_members   members that occur on at least one small plasmid
    n_large_members   members that occur on at least one large plasmid (a protein on
                      both counts in both)
    n_not_searched    members the cascade did not search (cascade_selection): they are
                      neither dark nor annotated
    scope             small_only_known | small_only_unknown | mixed_known | mixed_unknown
    known_from        small | large | both | '' - where the named members sit; a named
                      protein on small and large plasmids counts on both sides

mixed means a member on a large plasmid. known means a member the cascade or PlasmidScope
named (functional_class other than NONE, UNCHARACTERIZED_HOMOLOG and NOT_SEARCHED): a family
is as bright as its brightest member (Durairaj et al., Nature 2023, 622:646). Within every
family holding an unexplained small-plasmid protein, all members were annotated by
the same cascade (cascade_selection), so at the primary resolution known and unknown mean
the same on both sides.

n_members AND n_orfs ARE DIFFERENT NUMBERS

Clustering runs on the dereplicated set, so a protein whose sequence is identical on two
hundred plasmids is ONE member: it clusters alone and is labelled ORPHAN. That is a
defensible definition - it is not a family of divergent homologs - but `n_members: 1` alone
reads as "seen once", and a reader could not tell a genuine singleton from one of the most
widely carried proteins in the collection. n_orfs is the number of gene copies behind the
cluster.
"""
import collections
import csv

import _ctx  # noqa: F401

from darkorf import ids, status

cfg = snakemake.params.clustering
# Classes that name nothing, or were never looked at. Everything else is a named member.
UNNAMED = {"NONE", "UNCHARACTERIZED_HOMOLOG", "NOT_SEARCHED"}

# ------------------------------------------------------------------------------------
# Which proteins are dark, so the per-family dark fraction can be computed.
#
# Read from the cascade's own classification rather than recomputed: a family's dark
# fraction must mean the same thing as the dark set itself, and two definitions of dark
# would diverge the moment one threshold moved.
# ------------------------------------------------------------------------------------
dark = {line.strip() for line in open(snakemake.input.dark_ids) if line.strip()}

functional_class = {}
with open(snakemake.input.prot, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        functional_class[row["seq_id"]] = row["functional_class"]

small_plasmids = {line.strip() for line in open(snakemake.input.small_ids) if line.strip()}

# ------------------------------------------------------------------------------------
# Plasmid provenance per protein. Distribution is measured over INDEPENDENT units, not raw
# counts: spec section 2.6 - "a protein appearing many times in sequence databases does not
# mean that the observations are biologically independent".
# ------------------------------------------------------------------------------------
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

plasmid_meta = {}
with open(snakemake.input.registry, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        plasmid_meta[row["plasmid_id"]] = row

# Stage 6 lineage clusters: how many INDEPENDENT plasmid lineages a family occurs on.
# Distinct from family_plasmid_count, which counts records, and from family_MOB_count,
# which counts relaxase types - spec section 33.2 and design principle 2.6.
lineage_of = {}
with open(snakemake.input.lineage, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        lineage_of[row["plasmid_id"]] = row["plasmid_lineage_cluster"]


COLS = [
    "family_id", "family_resolution", "representative",
    # section 31.2
    "family_size", "n_orfs", "dark_member_count", "annotated_member_count",
    # unnamed members the dark set excludes (artefact, partial-only): neither dark nor named
    "n_unnamed_excluded",
    "percentage_dark_in_family", "family_class", "dark_only",
    # section 31.3
    "family_plasmid_count", "family_host_count",
    "family_genus_count", "family_MOB_count", "family_plasmid_lineage_count",
    "family_plasmid_lineage_status", "family_habitat_count",
    "family_small_plasmid_count",
    "n_small_members", "n_large_members", "n_not_searched", "scope", "known_from",
    "members",
]

# The DERIVED dark-family table (spec section 32). Downstream dark analysis - dN/dS,
# structure, context - operates on families that CONTAIN dark members, at the primary
# resolution only, and carries the dark_only flag with them.
#
# Containing-dark rather than dark-only, deliberately. Section 32 makes dark-only a
# DESCRIPTIVE label, not a filter, and a family that is 60% dark still holds dark proteins
# whose evolution and context are worth measuring - with the advantage that its annotated
# members say what the family does. Restricting to 100% dark would discard exactly the
# families where a dark protein is most interpretable.
primary = cfg["primary"]
dark_rows = []

summary = []
with open(snakemake.output.families, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=COLS, delimiter="\t")
    writer.writeheader()

    for resolution, clusters in zip(sorted(cfg["resolutions"]),
                                    sorted(snakemake.input.clusters)):
        # protein_clustering wrote families_<resolution>_cluster.tsv; sorted names pair
        # with sorted resolutions, and the assertion keeps it that way.
        assert f"families_{resolution}_cluster.tsv" in clusters, (resolution, clusters)
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

            meta = [plasmid_meta.get(p, {}) for p in plasmids]
            species = {m.get("species") for m in meta if m.get("species")}
            genera = {m.get("genus") for m in meta if m.get("genus")}
            mobs = {m.get("mob_cluster") for m in meta if m.get("mob_cluster")}
            habitats = {m.get("hab_top") for m in meta if m.get("hab_top")}
            lineages = {lineage_of[p] for p in plasmids if p in lineage_of}

            n_dark = sum(1 for m in mem if m in dark)
            # Members something NAMED. An unnamed protein the dark set excludes (an artefact,
            # a partial-only protein) is neither dark nor annotated; it is counted apart
            # rather than inflating the annotated count, which it used to.
            n_annotated = len(named)
            n_unnamed_excluded = len(mem) - n_dark - n_annotated - n_not_searched
            pct_dark = round(100.0 * n_dark / len(mem), 2)

            # Section 32: a dark-only family has no named member. Compared on the COUNT
            # rather than on the rounded percentage, because 99.996% rounds to 100.0 and a
            # family with one annotated member is not dark-only.
            dark_only = int(n_annotated == 0)
            n_dark_only += dark_only
            family_class = "ORPHAN" if len(mem) == 1 else "FAMILY"
            scope = (("mixed" if any(m in on_large for m in mem) else "small_only")
                     + ("_known" if named else "_unknown"))
            known_from = ("both" if len(known_from) == 2 else
                          known_from.pop() if known_from else "")
            n_orphan += family_class == "ORPHAN"

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
                "family_plasmid_count": len(plasmids),
                # The host is the species the plasmid was recovered from.
                "family_host_count": len(species),
                "family_genus_count": len(genera),
                "family_MOB_count": len(mobs),
                "family_plasmid_lineage_count": len(lineages),
                "family_plasmid_lineage_status": (
                    status.SUCCESS if lineages else status.NO_HIT),
                "family_habitat_count": len(habitats),
                "family_small_plasmid_count": len(plasmids & small_plasmids),
                "n_small_members": len(small),
                "n_large_members": sum(1 for m in mem if m in on_large),
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
                    # The representative must be a dark small-plasmid member: S8d searches
                    # this sequence structurally, and searching any other member would
                    # spend the ProstT5 budget on a protein outside the study.
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
                    "n_large_members": sum(1 for m in mem if m in on_large),
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

if not summary:
    raise SystemExit("protein_families: no resolutions configured - clustering.resolutions "
                     "is empty, so no family table was produced.")
