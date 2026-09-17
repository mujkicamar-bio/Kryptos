"""Stage 5: cluster EVERY unique protein into families with MMseqs2 (spec section 31).

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

n_members AND n_orfs ARE DIFFERENT NUMBERS

Clustering runs on the dereplicated set, so a protein whose sequence is identical on two
hundred plasmids is ONE member: it clusters alone and is labelled ORPHAN. That is a
defensible definition - it is not a family of divergent homologs - but `n_members: 1` alone
reads as "seen once", and a reader could not tell a genuine singleton from one of the most
widely carried proteins in the collection. n_orfs is the number of gene copies behind the
cluster.
"""
import _ctx  # noqa: F401
import collections
import csv
import pathlib
import subprocess

from darkorf import ids, status

cfg = snakemake.params.clustering
out_dir = pathlib.Path(snakemake.output.families).parent

# ------------------------------------------------------------------------------------
# Which proteins are dark, so the per-family dark fraction can be computed.
#
# Read from the cascade's own classification rather than recomputed: a family's dark
# fraction must mean the same thing as the dark set itself, and two definitions of dark
# would diverge the moment one threshold moved.
# ------------------------------------------------------------------------------------
dark = {line.strip() for line in open(snakemake.input.dark_ids) if line.strip()}

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

plasmid_meta = {}
with open(snakemake.input.registry, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        plasmid_meta[row["plasmid_id"]] = row


def genus_of(species):
    """The genus is the first token of a binomial.

    Crude and deliberately so: the alternative is a taxonomy lookup, and the count this
    feeds is a breadth DESCRIPTION rather than a classification. A species string that is
    not a binomial contributes itself, which over-counts genera rather than silently
    merging unrelated organisms.
    """
    return (species or "").split()[0] if species else ""


COLS = [
    "family_id", "family_resolution", "representative",
    # section 31.2
    "family_size", "n_orfs", "dark_member_count", "annotated_member_count",
    "percentage_dark_in_family", "family_class", "dark_only",
    # section 31.3
    "family_plasmid_count", "family_host_count", "family_species_count",
    "family_genus_count", "family_MOB_count", "family_plasmid_lineage_count",
    "family_plasmid_lineage_status", "family_habitat_count",
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

    for resolution, thresholds in sorted(cfg["resolutions"].items()):
        tmp = out_dir / f"mmseqs_tmp_{resolution}"
        prefix = str(out_dir / f"families_{resolution}")

        subprocess.run(
            f"mmseqs easy-cluster {snakemake.input.faa} {prefix} {tmp} "
            f"--min-seq-id {thresholds['min_seq_id']} -c {thresholds['coverage']} "
            f"--cov-mode {cfg['cov_mode']} --cluster-mode {cfg['cluster_mode']} "
            f"--cluster-reassign --threads {snakemake.threads} -v 1",
            shell=True, check=True)

        # easy-cluster writes <prefix>_cluster.tsv as representative<TAB>member.
        members = collections.defaultdict(list)
        with open(f"{prefix}_cluster.tsv") as fh:
            for line in fh:
                rep, mem = line.rstrip("\n").split("\t")
                members[rep].append(mem)

        n_dark_only = n_orphan = 0
        for rep, mem in sorted(members.items()):
            plasmids = set()
            n_orfs = 0
            for m in mem:
                plasmids |= seq_to_plasmids.get(m, set())
                n_orfs += seq_to_orf_count.get(m, 0)

            meta = [plasmid_meta.get(p, {}) for p in plasmids]
            species = {m.get("species") for m in meta if m.get("species")}
            mobs = {m.get("mob_cluster") for m in meta if m.get("mob_cluster")}
            habitats = {m.get("hab_top") for m in meta if m.get("hab_top")}

            n_dark = sum(1 for m in mem if m in dark)
            n_annotated = len(mem) - n_dark
            pct_dark = round(100.0 * n_dark / len(mem), 2)

            # Section 32: a dark-only family is exactly 100% dark. Compared on the COUNT
            # rather than on the rounded percentage, because 99.996% rounds to 100.0 and a
            # family with one annotated member is not dark-only.
            dark_only = int(n_annotated == 0)
            n_dark_only += dark_only
            family_class = "ORPHAN" if len(mem) == 1 else "FAMILY"
            n_orphan += family_class == "ORPHAN"

            writer.writerow({
                "family_id": ids.family_id(resolution, rep),
                "family_resolution": resolution,
                "representative": rep,
                "family_size": len(mem),
                "n_orfs": n_orfs,
                "dark_member_count": n_dark,
                "annotated_member_count": n_annotated,
                "percentage_dark_in_family": pct_dark,
                "family_class": family_class,
                "dark_only": dark_only,
                "family_plasmid_count": len(plasmids),
                # The host IS the organism the plasmid was recovered from, which this
                # collection records as `species`. Reported under both names because
                # section 31.3 asks for both and they are the same measurement here;
                # collapsing them would silently drop a field the specification names.
                "family_host_count": len(species),
                "family_species_count": len(species),
                "family_genus_count": len({genus_of(s) for s in species if genus_of(s)}),
                "family_MOB_count": len(mobs),
                # Plasmid lineage clusters are Stage 6 (spec section 33) and do not exist
                # yet. NOT_RUN rather than 0: "no lineages counted" and "lineage clustering
                # has not been performed" are different statements (section 2.9).
                "family_plasmid_lineage_count": "",
                "family_plasmid_lineage_status": status.NOT_RUN,
                "family_habitat_count": len(habitats),
                "members": ",".join(mem),
            })

            if resolution == primary and n_dark:
                dark_members = [m for m in mem if m in dark]
                dark_rows.append({
                    "family_id": ids.family_id(resolution, rep),
                    # The representative must be a DARK member: S8d searches this sequence
                    # structurally, and searching an annotated representative would spend
                    # the ProstT5 budget on a protein that is not part of the dark set.
                    "representative": rep if rep in dark else dark_members[0],
                    "n_members": len(mem),
                    "n_orfs": n_orfs,
                    "n_plasmids": len(plasmids),
                    "n_mob_clusters": len(mobs),
                    "family_class": family_class,
                    "dark_member_count": n_dark,
                    "annotated_member_count": n_annotated,
                    "percentage_dark_in_family": pct_dark,
                    "dark_only": dark_only,
                    # Only the dark members: the downstream stages measure the DARK
                    # proteins' evolution and context, and including annotated members
                    # here would silently widen every one of those measurements.
                    "members": ",".join(dark_members),
                })

        summary.append((resolution, len(members), n_orphan, n_dark_only))

DARK_COLS = ["family_id", "representative", "n_members", "n_orfs", "n_plasmids",
             "n_mob_clusters", "family_class", "dark_member_count",
             "annotated_member_count", "percentage_dark_in_family", "dark_only", "members"]
with open(snakemake.output.dark_families, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=DARK_COLS, delimiter="\t")
    writer.writeheader()
    writer.writerows(dark_rows)

print(f"{'resolution':<14}{'families':>10}{'orphans':>10}{'dark_only':>11}")
for resolution, n_fam, n_orphan, n_dark_only in summary:
    print(f"{resolution:<14}{n_fam:>10}{n_orphan:>10}{n_dark_only:>11}")

print(f"dark families at {primary}: {len(dark_rows)} containing at least one dark member, "
      f"of which {sum(r['dark_only'] for r in dark_rows)} are 100% dark")

if not summary:
    raise SystemExit("protein_families: no resolutions configured - clustering.resolutions "
                     "is empty, so no family table was produced.")
