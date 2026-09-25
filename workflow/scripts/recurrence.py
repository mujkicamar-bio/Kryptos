"""Stage 7: distribution and recurrence, counted over independent units (spec section 34).

THE ONE RULE THIS STAGE EXISTS TO ENFORCE

Section 34.2, last line: "Database record counts must never be treated as independent
biological observations." Section 2.6 says the same as a design principle. A protein on two
thousand plasmid records may be one clinical plasmid deposited two thousand times, and every
claim about how widespread it is inherits that error unless the counts are kept apart.

So this stage reports SEVEN different counts for each family and never collapses them:

    plasmid_occurrence_count          gene copies. The largest number, and the least
                                      informative on its own.
    unique_plasmid_count              distinct plasmid records.
    independent_plasmid_cluster_count distinct Stage 6 lineages. THIS is the denominator a
                                      recurrence claim needs.
    host_count / species_count        distinct host species (binomials).
    genus_count                       distinct host genera, including hosts named only to
                                      the genus ("Acidovorax sp.").
    MOB_count                         distinct relaxase types. Breadth of mobility, not of
                                      evolution - section 33 keeps these separate.
    habitat_count                     distinct environments.

Section 33.2's worked example is the shape to expect: 2,143 occurrences, 8 MOB clusters, 47
independent clusters. The three numbers describe different biological properties and a
reader who is shown only the first will draw the wrong conclusion.

DATABASE RECURRENCE (section 34.1)

    database_source_count   how many source databases contributed the records
    database_record_count   the raw record count

These are PROVENANCE, not biology. They are reported because section 34.1 asks for them and
because a reader should be able to see that a number is large for a database reason; they
are never the denominator of anything.

WHY THIS IS A SEPARATE STAGE AND NOT A COLUMN ON THE FAMILY TABLE

The family table describes composition - what is in this family. This describes
distribution - where it has been seen, and how independently. Both are per family, and it
would be easy to merge them; keeping them apart is what stops the independent count being
quietly replaced by the occurrence count when someone needs "a number for how common this
is".
"""
import collections
import csv

import _ctx  # noqa: F401

from darkorf import status

# ------------------------------------------------------------------------------------
# Per-protein provenance: which plasmids, and how many gene copies.
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

# ------------------------------------------------------------------------------------
# Per-plasmid metadata and lineage.
# ------------------------------------------------------------------------------------
meta_of = {}
with open(snakemake.input.registry, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        meta_of[row["plasmid_id"]] = row

lineage_of = {}
with open(snakemake.input.lineage, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        lineage_of[row["plasmid_id"]] = row["plasmid_lineage_cluster"]

# Source databases per plasmid, for the provenance counts. The master table records these
# as a delimited `sources` string; absence is absence and contributes nothing.
sources_of = {}
with open(snakemake.input.master, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        raw = (row.get("sources") or "").strip()
        if raw:
            sources_of[row["plasmid_id"]] = {
                s.strip() for s in raw.replace(";", ",").split(",") if s.strip()}


COLS = [
    "family_id", "family_resolution", "representative",
    # section 34, the seven biological counts
    "plasmid_occurrence_count", "unique_plasmid_count",
    "independent_plasmid_cluster_count", "independent_cluster_status",
    "host_count", "species_count", "genus_count",
    # how many of the family's plasmids have a recorded host; host_count is NOT_MEASURED,
    # not 0, when none has (clonal_registry, plasmidann.hosts)
    "n_plasmids_with_host", "n_plasmids_with_species", "host_count_status",
    # MOB-suite's predicted host range over EVERY plasmid, hosted or not: a separate
    # measurement at any rank, never counted as a host (clonal_registry)
    "n_plasmids_with_predicted_range", "predicted_host_range_count",
    "predicted_host_ranges",
    "MOB_count", "habitat_count",
    # section 34.1, provenance - never a denominator
    "database_record_count", "database_source_count",
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
            # A family whose plasmids are all missing from the lineage table has not been
            # measured for independence; reporting 0 would read as "no independent
            # lineages", which is a much stronger and quite different claim (section 2.9).
            "independent_cluster_status": status.SUCCESS if lineages else status.NOT_RUN,
            # host and species are the same measurement in this collection; section 34
            # names both, and collapsing them would drop a field the spec asks for.
            "host_count": len(species),
            "species_count": len(species),
            "genus_count": len(genera),
            "n_plasmids_with_host": n_with_host,
            # SINGLE_HOST (Stage 14) needs every plasmid named to the species.
            "n_plasmids_with_species": sum(1 for m in meta if m.get("species")),
            "host_count_status": status.SUCCESS if n_with_host else status.NOT_MEASURED,
            "n_plasmids_with_predicted_range": len(predicted),
            "predicted_host_range_count": len(set(predicted)),
            # ';' between ranges: a single range may itself list several phyla with ','.
            "predicted_host_ranges": ";".join(sorted(set(predicted))),
            "MOB_count": len(mobs),
            "habitat_count": len(habitats),
            "database_record_count": len(plasmids),
            "database_source_count": len(sources),
        })
        n_rows += 1

        if lineages:
            inflation = len(plasmids) / len(lineages)
            if inflation > max_inflation[0]:
                max_inflation = (inflation, family["family_id"])

print(f"recurrence: {n_rows} families")
# The headline diagnostic: how far raw plasmid counts overstate independence for the worst
# family in the run. A ratio near 1 means records and lineages agree; a large one means
# some family's apparent prevalence is mostly redeposition.
if max_inflation[1]:
    print(f"            largest record-to-lineage ratio {max_inflation[0]:.1f}x "
          f"({max_inflation[1]}) - plasmid records overstate independence by this much")

if not n_rows:
    raise SystemExit("recurrence: no families read - every distribution count would be "
                     "absent, and the family table is upstream of this stage.")
