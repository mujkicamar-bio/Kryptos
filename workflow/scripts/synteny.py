"""Stage 9: synteny and context conservation across a cluster's lineages.

See src/plasmidann/synteny.py for what the measurements are, why they are kept apart, and
why they are counted over plasmid lineages (Stage 6) rather than occurrences. This script
builds each cluster's occurrence list from the annotation table and the protein map, tags
every occurrence with its plasmid's lineage, then measures.

It reads the SAME neighbourhood definition as Stage 8 - plasmidann.context.directons and
flanks - rather than restating it. Two definitions of "neighbour" in one pipeline would
diverge the moment one window changed, and the two stages' numbers would stop being
comparable without anything saying so.

TWO LEVELS (synteny.levels)

    close         a gene is its close cluster (90% identity); neighbours are named by
                  their close cluster: "the same gene"
    intermediate  a gene is its protein family (clustering.primary); neighbours are named
                  by their family: "the same family"

Each level gives one row per cluster holding a dark small-plasmid member - the rule
protein_families.py uses for dark_families.tsv, applied at that level's resolution - and
the occurrences measured are those of the cluster's DARK members. At the primary level the
derived set must equal dark_families.tsv, so the two definitions cannot drift silently.
Neighbour identity at the close level is strict, so close-level conservation is expected to
be at most the family-level value for the same arrangement. MMseqs2 clusterings at
different thresholds are not nested, so a close row lists the primary families of its dark
members in intermediate_family_ids.

A NEIGHBOUR IS IDENTIFIED BY ITS CLUSTER, NOT ITS LABEL

Every gene has one, dark genes included, and it means the same on small and large
plasmids, whose genes were not all annotated alike.
modal_left, modal_right and modal_synteny therefore hold cluster ids, and left and right
are upstream and downstream on the gene's own strand. operon_like still asks for a
FUNCTIONAL partner, as Stage 8 does: an uncharacterised partner says nothing about what the
arrangement is for.

TWO SETS OF OCCURRENCES

Every measurement is reported over all occurrences of the cluster's dark members, and again
with the prefix small_ over the occurrences on small plasmids alone. On a small plasmid the
+-3 window often reaches every other gene of the molecule (a quarter of small-plasmid ORFs
on the test set), so small_lineage_neighborhood_conservation there measures conservation of
gene content rather than of local order; the immediate-neighbour measurements do not depend
on the window.
"""
import collections
import csv
import pathlib

import _ctx  # noqa: F401

from darkorf import ids
from darkorf.circular import is_circular
from plasmidann.context import directons, flanks
from plasmidann.synteny import conservation

cfg = snakemake.params.context
levels = snakemake.params.synteny["levels"]
min_lineages = snakemake.params.synteny["min_lineages"]
primary = snakemake.params.primary
if primary not in levels:
    raise SystemExit(f"synteny.levels {levels} must include clustering.primary ({primary}): "
                     "the family table reads the primary-level rows")

# One cluster file per level, in the order of synteny.levels.
cluster_of = {}
for level, path in zip(levels, snakemake.input.clusters, strict=True):
    if pathlib.Path(path).name != f"families_{level}_cluster.tsv":
        raise SystemExit(f"synteny: cluster file {path} does not match level {level}")
    cluster_of[level] = {}
    last_rep = None
    with open(path) as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")
            if rep != last_rep:
                # The file is grouped by representative, so one id string serves the group.
                last_rep, cluster_id = rep, ids.family_id(level, rep)
            cluster_of[level][mem] = cluster_id

dark = {l.strip() for l in open(snakemake.input.dark_ids) if l.strip()}
small_plasmids = {l.strip() for l in open(snakemake.input.small_ids) if l.strip()}
with open(snakemake.input.lineage, newline="") as fh:
    lineage_of = {r["plasmid_id"]: r["plasmid_lineage_cluster"]
                  for r in csv.DictReader(fh, delimiter="\t")}
circular = set()
with open(snakemake.input.registry, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if is_circular(r.get("topology")):
            circular.add(r["plasmid_id"])
with open(snakemake.input.lengths, newline="") as fh:
    length_of = {r["plasmid_id"]: int(r["length_bp"])
                 for r in csv.DictReader(fh, delimiter="\t")}

# ------------------------------------------------------------------------------------
# Every ORF per plasmid, with the label the cascade gave it.
# ------------------------------------------------------------------------------------
by_plasmid = collections.defaultdict(list)
class_of = {}
with open(snakemake.input.annotation, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        by_plasmid[r["plasmid_id"]].append({
            "orf_id": r["orf_id"], "start": int(r["start"]), "end": int(r["end"]),
            "strand": 1 if r["strand"] in ("1", "+") else -1})
        class_of[r["orf_id"]] = r.get("functional_class") or ""

# Stage 6 assigns every plasmid of the analysis set a lineage. A plasmid without one means
# the inputs come from different runs; guessing its own lineage would inflate independence.
missing = sorted(set(by_plasmid) - set(lineage_of))
if missing:
    raise SystemExit(f"synteny: {len(missing)} plasmids have no lineage in "
                     f"{snakemake.input.lineage}, e.g. {missing[:3]}")

orfs_of_seq = collections.defaultdict(list)
seq_of_orf = {}
on_small = set()
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, orf_ids = line.rstrip("\n").split("\t")
        for orf_id in orf_ids.split(","):
            orfs_of_seq[sid].append(orf_id)
            seq_of_orf[orf_id] = sid
            if orf_id.rsplit("|", 1)[0] in small_plasmids:
                on_small.add(sid)

# ------------------------------------------------------------------------------------
# Ordered neighbour ORFs per ORF, computed once per plasmid and named per level only for
# the occurrences measured: two translated copies for every ORF of the full run would
# double the dominant memory cost for rows that are never read.
#
# Left and right are kept apart, since their conservation is measured separately, and are
# relative to the gene: upstream and downstream on its own strand, so that the same
# arrangement written in the opposite orientation is the same arrangement. On a circular
# plasmid the window and the directons wrap across the origin
# (plasmidann.context.flanks, plasmidann.context.directons).
# ------------------------------------------------------------------------------------
context_of = {}
window = cfg["neighbourhood_window"]
for plasmid_id, genes in by_plasmid.items():
    units = directons(genes, max_gap=cfg["max_operon_gap"],
                      circular=plasmid_id in circular, length=length_of.get(plasmid_id))
    unit_of = {oid: i for i, unit in enumerate(units) for oid in unit}
    flanks_of = flanks(genes, window, circular=plasmid_id in circular)

    for gene in genes:
        orf_id = gene["orf_id"]
        # Nearest-first on both sides, so [0] is the immediate neighbour.
        left, right = flanks_of[orf_id]
        if gene["strand"] == -1:
            left, right = right, left

        unit = units[unit_of[orf_id]] if orf_id in unit_of else [orf_id]
        # Operon-like means sharing a transcriptional unit with a FUNCTIONAL partner,
        # the same definition context_features uses (operon_with_annotated). An
        # uncharacterised partner says nothing about what the arrangement is for.
        operon = any(class_of.get(p) == "FUNCTIONAL" for p in unit if p != orf_id)
        context_of[orf_id] = (tuple(left), tuple(right), operon, lineage_of[plasmid_id])


def occurrence(orf_id, level):
    """One occurrence for plasmidann.synteny.conservation, neighbours named at `level`."""
    left, right, operon, lineage = context_of[orf_id]
    name = cluster_of[level]
    return {"left": [name.get(seq_of_orf.get(o, ""), "") for o in left],
            "right": [name.get(seq_of_orf.get(o, ""), "") for o in right],
            "operon": operon, "lineage": lineage}


MEASURES = ["n_occurrences", "context_recurrence", "n_lineages", "n_lineages_discordant",
            "lineage_left_conservation", "lineage_right_conservation",
            "lineage_neighborhood_conservation", "lineage_operon_like_conservation",
            "lineage_synteny_conservation", "modal_left", "modal_right", "modal_synteny",
            "status"]
COLS = ["family_id", "level", "intermediate_family_ids", "synteny_min_lineages",
        *MEASURES, *(f"small_{m}" for m in MEASURES)]

with open(snakemake.input.families, newline="") as fh:
    dark_families = {r["family_id"] for r in csv.DictReader(fh, delimiter="\t")}

with open(snakemake.output.tsv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=COLS, delimiter="\t")
    writer.writeheader()

    for level in levels:
        # The dark members of every cluster that holds a dark small-plasmid member.
        dark_members = collections.defaultdict(list)
        for sid in sorted(dark):
            if sid in cluster_of[level]:
                dark_members[cluster_of[level][sid]].append(sid)
        measured = sorted(c for c, mem in dark_members.items()
                          if any(m in on_small for m in mem))
        if level == primary and set(measured) != dark_families:
            raise SystemExit(
                f"synteny: the {level} clusters holding a dark small-plasmid member differ "
                f"from dark_families ({len(set(measured) - dark_families)} only here, "
                f"{len(dark_families - set(measured))} only in {snakemake.input.families})")

        n_success = n_conserved = 0
        for cluster_id in measured:
            occurrences, small = [], []
            for member in dark_members[cluster_id]:
                for orf_id in orfs_of_seq.get(member, ()):
                    if orf_id in context_of:
                        occ = occurrence(orf_id, level)
                        occurrences.append(occ)
                        if orf_id.rsplit("|", 1)[0] in small_plasmids:
                            small.append(occ)

            result = conservation(occurrences, min_lineages)
            result.update({f"small_{k}": v
                           for k, v in conservation(small, min_lineages).items()})
            result.update({
                "family_id": cluster_id, "level": level, "synteny_min_lineages": min_lineages,
                "intermediate_family_ids": ",".join(sorted(
                    {cluster_of[primary].get(m, "") for m in dark_members[cluster_id]}))})
            writer.writerow({k: result.get(k, "") for k in COLS})
            if result["status"] == "SUCCESS":
                n_success += 1
                n_conserved += result["lineage_synteny_conservation"] >= 0.9

        print(f"synteny: {level}: {len(measured)} clusters, {n_success} measured over "
              f">= {min_lineages} lineages, {n_conserved} with lineage synteny "
              "conservation >= 0.9")
