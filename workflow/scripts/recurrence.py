"""Distribution and recurrence of each family, counted over independent units.

Inputs: protein_families.tsv, the protein map, the clonal registry, plasmid_lineage.tsv
and the master table. Output recurrence.tsv, one row per family, with counts kept apart
because a protein on two thousand plasmid records may be one plasmid deposited two
thousand times:

  plasmid_occurrence_count           gene copies
  unique_plasmid_count               distinct plasmid records
  independent_plasmid_cluster_count  distinct lineages, the denominator of a
                                     recurrence claim
  host_count, genus_count            distinct observed host species and genera
                                     (plasmidann.hosts)
  predicted_host_range_*             MOB-suite's predicted host ranges, never a host
  MOB_count                          distinct MOB-suite clusters
  habitat_count                      distinct habitats (hab_top)
  database_source_count              source databases of the records; provenance, never a
                                     denominator
"""
import collections
import csv

import _ctx  # noqa: F401

from darkorf import status

seq_to_plasmids = collections.defaultdict(set)
seq_to_orf_count = collections.Counter()
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, orf_ids = line.rstrip("\n").split("\t")
        members = orf_ids.split(",")
        seq_to_orf_count[sid] = len(members)
        for orf_id in members:
            seq_to_plasmids[sid].add(orf_id.rsplit("|", 1)[0])

meta_of = {}
with open(snakemake.input.registry, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        meta_of[row["plasmid_id"]] = row

lineage_of = {}
with open(snakemake.input.lineage, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        lineage_of[row["plasmid_id"]] = row["plasmid_lineage_cluster"]

# Source databases per plasmid: the master table's `sources`, delimited by ',' or ';'.
sources_of = {}
with open(snakemake.input.master, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        raw = (row.get("sources") or "").strip()
        if raw:
            sources_of[row["plasmid_id"]] = {
                s.strip() for s in raw.replace(";", ",").split(",") if s.strip()}


COLS = [
    "family_id", "family_resolution", "representative",
    "plasmid_occurrence_count", "unique_plasmid_count",
    "independent_plasmid_cluster_count", "independent_cluster_status",
    "host_count", "genus_count",
    # n_plasmids_with_host counts plasmids with a host named to the genus or finer, and
    # host_count_status is NOT_MEASURED when there is none. With genus-level hosts only,
    # the status is SUCCESS and host_count is 0 while genus_count is not.
    "n_plasmids_with_host", "n_plasmids_with_species", "host_count_status",
    # MOB-suite's predicted host range over EVERY plasmid, hosted or not: a separate
    # measurement at any rank, never counted as a host (clonal_registry)
    "n_plasmids_with_predicted_range", "predicted_host_range_count",
    "predicted_host_ranges",
    "MOB_count", "habitat_count",
    "database_source_count",
]

n_rows = 0
max_inflation = (0.0, "")
with open(snakemake.input.families, newline="") as fh, \
        open(snakemake.output.tsv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=COLS, delimiter="\t")
    writer.writeheader()

    for family in csv.DictReader(fh, delimiter="\t"):
        members = [m for m in family["members"].split(",") if m]
        plasmids = set()
        occurrences = 0
        for member in members:
            plasmids |= seq_to_plasmids.get(member, set())
            occurrences += seq_to_orf_count.get(member, 0)

        meta = [meta_of.get(p, {}) for p in plasmids]
        species = {m.get("species") for m in meta if m.get("species")}
        genera = {m.get("genus") for m in meta if m.get("genus")}
        n_with_host = sum(1 for m in meta if m.get("genus"))
        predicted = [m.get("predicted_host_range") for m in meta
                     if m.get("predicted_host_range")]
        mobs = {m.get("mob_cluster") for m in meta if m.get("mob_cluster")}
        habitats = {m.get("hab_top") for m in meta if m.get("hab_top")}
        lineages = {lineage_of[p] for p in plasmids if p in lineage_of}
        sources = set().union(*(sources_of.get(p, set()) for p in plasmids)) \
            if plasmids else set()

        writer.writerow({
            "family_id": family["family_id"],
            "family_resolution": family.get("family_resolution", ""),
            "representative": family.get("representative", ""),
            "plasmid_occurrence_count": occurrences,
            "unique_plasmid_count": len(plasmids),
            "independent_plasmid_cluster_count": len(lineages),
            # No member plasmid in the lineage table: independence was not measured.
            "independent_cluster_status": status.SUCCESS if lineages else status.NOT_RUN,
            "host_count": len(species),
            "genus_count": len(genera),
            "n_plasmids_with_host": n_with_host,
            # SINGLE_HOST (rule rarity) needs every plasmid named to the species.
            "n_plasmids_with_species": sum(1 for m in meta if m.get("species")),
            "host_count_status": status.SUCCESS if n_with_host else status.NOT_MEASURED,
            "n_plasmids_with_predicted_range": len(predicted),
            "predicted_host_range_count": len(set(predicted)),
            # ';' between ranges: a single range may itself list several phyla with ','.
            "predicted_host_ranges": ";".join(sorted(set(predicted))),
            "MOB_count": len(mobs),
            "habitat_count": len(habitats),
            "database_source_count": len(sources),
        })
        n_rows += 1

        if lineages:
            inflation = len(plasmids) / len(lineages)
            if inflation > max_inflation[0]:
                max_inflation = (inflation, family["family_id"])

print(f"recurrence: {n_rows} families")
# The largest record-to-lineage ratio of any family: how far record counts overstate
# independence in this run.
if max_inflation[1]:
    print(f"            largest record-to-lineage ratio {max_inflation[0]:.1f}x "
          f"({max_inflation[1]}) - plasmid records overstate independence by this much")

if not n_rows:
    raise SystemExit("recurrence: no families read - every distribution count would be "
                     "absent, and the family table is upstream of this stage.")
