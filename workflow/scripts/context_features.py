"""S8c: genomic context per ORF, aggregated to one row of rates per family.

WHAT IT MEASURES

For every dark family, the fraction of the plasmids carrying it on which a member ORF

    cons_defence                 sits inside a defence system (DefenseFinder)
    cons_integron                sits inside an integron cassette array (IntegronFinder)
    cons_is_element              sits inside an IS element (ISEScan)
    cons_annotated_neighbour     has a FUNCTIONAL gene within the +-3 neighbourhood
    cons_operon_with_annotated   shares a directon with a FUNCTIONAL gene
    cons_two_gene_operon         is one of exactly two genes in its directon

These are descriptive rates, not tests. There is no background and no enrichment: the
context enrichment test, its stratified background and the label grouping it ran on were
removed on 2026-09-25. What each dark ORF's neighbours are is still described per family by
synteny (S9).

A rate of 0 is a measurement - the family was examined and the feature was absent - not a
missing value. The rates are not corrected for plasmid size: on a six-gene plasmid a +-3
neighbourhood is the whole molecule, so cons_annotated_neighbour is high there by
construction.

THE UNIT IS THE PLASMID

Not the family member. Members of a family are homologs on plasmids that are frequently
near-identical, so counting members makes sequencing effort look like evidence. Clonal
redundancy between distinct plasmids remains uncorrected; see docs/PARAMETER_PROVENANCE.md.
"""
import collections
import csv

import _ctx  # noqa: F401

from darkorf.circular import is_circular
from plasmidann.context import directons, neighbourhood, overlapping_islands

cfg = snakemake.params.context
FEATURES = ("defence", "integron", "is_element", "annotated_neighbour",
            "operon_with_annotated", "two_gene_operon")

# ------------------------------------------------------------------------------------
# Every ORF on every plasmid, with whatever the cascade named it.
# ------------------------------------------------------------------------------------
by_plasmid = collections.defaultdict(list)
class_of = {}
with open(snakemake.input.annotation, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        by_plasmid[r["plasmid_id"]].append({
            "orf_id": r["orf_id"], "start": int(r["start"]), "end": int(r["end"]),
            "strand": 1 if r["strand"] in ("1", "+") else -1})
        class_of[r["orf_id"]] = r.get("functional_class") or "NONE"

# ------------------------------------------------------------------------------------
# Which ORFs correspond to which unique protein.
#
# The reverse direction matters: the per-family loop below previously rescanned the whole
# 9.3M-entry forward map per family, which is O(families x ORFs) - measured at 0.37-0.39 s
# per family, projecting to 10 hours to 3.6 days single-core for ~50k families.
# ------------------------------------------------------------------------------------
orfs_of_seq = collections.defaultdict(list)
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, ids = line.rstrip("\n").split("\t")
        for oid in ids.split(","):
            orfs_of_seq[sid].append(oid)

# ------------------------------------------------------------------------------------
# Islands: defence systems, integron cassette arrays and IS elements, as intervals per
# plasmid.
#
# These are system-level calls rather than per-protein labels, so they stay intervals: a
# dark ORF INSIDE a defence system is a different statement from one merely beside a
# defence component.
# ------------------------------------------------------------------------------------
islands = collections.defaultdict(list)
defence_orfs = set()
with open(snakemake.input.defence, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get("orf_id"):
            defence_orfs.add(r["orf_id"])

with open(snakemake.input.integrons, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        try:
            islands[r["plasmid_id"]].append(
                {"name": "integron", "start": int(r["start"]), "end": int(r["end"])})
        except (ValueError, KeyError):
            continue

with open(snakemake.input.is_elements, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        islands[r["plasmid_id"]].append(
            {"name": "is_element", "start": int(r["start"]), "end": int(r["end"])})

for pid, genes in by_plasmid.items():
    for g in genes:
        if g["orf_id"] in defence_orfs:
            islands[pid].append({"name": "defence", "start": g["start"], "end": g["end"]})

# ------------------------------------------------------------------------------------
# Per-ORF context: the islands it sits inside, and its annotated neighbours and directon
# partners.
# ------------------------------------------------------------------------------------
# The neighbour window and the directons wrap across the origin of a circular plasmid
# (context.flanks, context.directons).
circular = set()
with open(snakemake.input.master, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if is_circular(r.get("topology")):
            circular.add(r["plasmid_id"])
with open(snakemake.input.lengths, newline="") as fh:
    length_of = {r["plasmid_id"]: int(r["length_bp"])
                 for r in csv.DictReader(fh, delimiter="\t")}

context_of = {}
for pid, genes in by_plasmid.items():
    units = directons(genes, max_gap=cfg["max_operon_gap"], circular=pid in circular,
                      length=length_of.get(pid))
    unit_of = {oid: i for i, unit in enumerate(units) for oid in unit}
    plasmid_islands = islands.get(pid, [])

    for g in genes:
        oid = g["orf_id"]
        # EVERY island, not the first one found: an ORF inside a defence system that also
        # sits in a cassette array is both.
        ctx = {island["name"] for island in overlapping_islands(g, plasmid_islands)}

        if any(class_of.get(n) == "FUNCTIONAL"
               for n in neighbourhood(genes, oid, window=cfg["neighbourhood_window"],
                                      circular=pid in circular)):
            ctx.add("annotated_neighbour")

        # Directon membership with at least one annotated partner: the dark ORF is
        # predicted to be co-transcribed with something we understand, which is a stronger
        # claim than adjacency.
        unit = units[unit_of[oid]] if oid in unit_of else [oid]
        partners = [p for p in unit if p != oid]
        if partners and any(class_of.get(p) == "FUNCTIONAL" for p in partners):
            ctx.add("operon_with_annotated")
        if len(unit) == 2:
            ctx.add("two_gene_operon")

        context_of[oid] = ctx

# ------------------------------------------------------------------------------------
# Per family: the fraction of its plasmids on which each feature is present.
# ------------------------------------------------------------------------------------
n_families = 0
with open(snakemake.output.families, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["family_id", "n_units"] + [f"cons_{f}" for f in FEATURES])
    with open(snakemake.input.families, newline="") as fh:
        for fam in csv.DictReader(fh, delimiter="\t"):
            # The unit is the plasmid. orf_id is '<plasmid_id>|<ordinal>'.
            per_plasmid = collections.defaultdict(set)
            for member in fam["members"].split(","):
                for oid in orfs_of_seq.get(member, ()):
                    per_plasmid[oid.rsplit("|", 1)[0]] |= context_of.get(oid, set())
            n_units = len(per_plasmid)
            if not n_units:
                continue
            n_families += 1
            w.writerow([fam["family_id"], n_units] + [
                round(sum(f in present for present in per_plasmid.values()) / n_units, 6)
                for f in FEATURES])

print(f"context: {n_families} families over {len(by_plasmid)} plasmids")
