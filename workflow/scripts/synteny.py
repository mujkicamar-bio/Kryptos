"""Stage 9: synteny and context conservation across a family's occurrences (section 42).

See src/plasmidann/synteny.py for what the six measurements are and why they are kept
apart. This script builds each family's occurrence list from the annotation table and the
protein map, then measures.

It reads the SAME neighbourhood definition as Stage 8 - plasmidann.context.directons and
neighbourhood - rather than restating it. Two definitions of "neighbour" in one pipeline
would diverge the moment one window changed, and the two stages' numbers would stop being
comparable without anything saying so.

A NEIGHBOUR IS IDENTIFIED BY ITS FAMILY, NOT ITS LABEL

Two neighbours are "the same" when their proteins share a family (S2f, the primary
clustering). Every gene
has one - dark genes included, which a label comparison dropped as empty - and it means the
same on small and large plasmids, whose genes were not all annotated alike. modal_left,
modal_right and modal_synteny therefore hold family ids, and left and right are upstream and
downstream on the gene's own strand. operon_like still asks for a FUNCTIONAL partner, as
Stage 8 does: an uncharacterised partner says nothing about what the arrangement is for.

TWO SETS OF OCCURRENCES

Every measurement is reported over all occurrences of the family's dark members, and again
with the prefix small_ over the occurrences on small plasmids alone.
"""
import collections
import csv

import _ctx  # noqa: F401

from darkorf import ids
from darkorf.circular import is_circular
from plasmidann.context import directons, flanks
from plasmidann.synteny import conservation

cfg = snakemake.params.context

family_of_seq = {}
with open(snakemake.input.clusters) as fh:
    for line in fh:
        rep, mem = line.rstrip("\n").split("\t")
        family_of_seq[mem] = ids.family_id(snakemake.params.primary, rep)
small_plasmids = {l.strip() for l in open(snakemake.input.small_ids) if l.strip()}
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

orfs_of_seq = collections.defaultdict(list)
family_of_orf = {}
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, orf_ids = line.rstrip("\n").split("\t")
        for orf_id in orf_ids.split(","):
            orfs_of_seq[sid].append(orf_id)
            family_of_orf[orf_id] = family_of_seq.get(sid, "")

# ------------------------------------------------------------------------------------
# Ordered neighbours per ORF, computed once per plasmid.
#
# Split into left and right around the ORF rather than returned as one list: section 42's
# left_neighbor_conservation and right_neighbor_conservation are separate measurements and
# cannot be recovered from a combined neighbourhood.
#
# LEFT AND RIGHT ARE RELATIVE TO THE GENE, upstream and downstream of its own strand, not
# the record's coordinates. The same arrangement written in the opposite orientation was
# otherwise counted as a different one (3 of 17 measured test families changed). On a
# circular plasmid the window and the directons wrap across the origin
# (plasmidann.context.flanks, plasmidann.context.directons).
# ------------------------------------------------------------------------------------
context_of = {}
window = cfg["neighbourhood_window"]
for plasmid_id, genes in by_plasmid.items():
    units = directons(genes, max_gap=cfg["max_operon_gap"],
                      circular=plasmid_id in circular, length=length_of.get(plasmid_id))
    unit_of = {oid: i for i, unit in enumerate(units) for oid in unit}

    for gene in genes:
        orf_id = gene["orf_id"]
        # Nearest-first on both sides, so [0] is the immediate neighbour.
        left, right = flanks(genes, orf_id, window, circular=plasmid_id in circular)
        if gene["strand"] == -1:
            left, right = right, left

        unit = units[unit_of[orf_id]] if orf_id in unit_of else [orf_id]
        partners = [p for p in unit if p != orf_id]
        context_of[orf_id] = {
            "left": [family_of_orf.get(o, "") for o in left],
            "right": [family_of_orf.get(o, "") for o in right],
            # Operon-like means sharing a transcriptional unit with a FUNCTIONAL partner,
            # the same definition context_features uses (operon_with_annotated). An
            # uncharacterised partner says nothing about what the arrangement is for.
            "operon": any(class_of.get(p) == "FUNCTIONAL" for p in partners),
        }

MEASURES = ["n_occurrences", "context_recurrence",
            "left_neighbor_conservation", "right_neighbor_conservation",
            "neighborhood_conservation", "operon_like_conservation",
            "synteny_conservation", "modal_left", "modal_right", "modal_synteny", "status"]
COLS = ["family_id", *MEASURES, *(f"small_{m}" for m in MEASURES)]

n_rows = n_conserved = 0
with open(snakemake.input.families, newline="") as fh, \
        open(snakemake.output.tsv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=COLS, delimiter="\t")
    writer.writeheader()

    for family in csv.DictReader(fh, delimiter="\t"):
        occurrences, small = [], []
        for member in family["members"].split(","):
            for orf_id in orfs_of_seq.get(member, ()):
                if orf_id in context_of:
                    occurrences.append(context_of[orf_id])
                    if orf_id.rsplit("|", 1)[0] in small_plasmids:
                        small.append(context_of[orf_id])

        result = conservation(occurrences)
        result.update({f"small_{k}": v for k, v in conservation(small).items()})
        result["family_id"] = family["family_id"]
        writer.writerow({k: result.get(k, "") for k in COLS})
        n_rows += 1
        if result["status"] == "SUCCESS" and result["synteny_conservation"] != "" \
                and float(result["synteny_conservation"]) >= 0.9:
            n_conserved += 1

print(f"synteny: {n_rows} families, {n_conserved} with synteny conservation >= 0.9")
# Conservation is computed over OCCURRENCES. A family on forty redepositions of one plasmid
# shows perfect synteny from one biological event, so this number must be read against the
# independent-lineage count from Stage 7 - which is why both tables carry family_id.
print("         read these against independent_plasmid_cluster_count from Stage 7: "
      "conservation over occurrences is not conservation over lineages")
