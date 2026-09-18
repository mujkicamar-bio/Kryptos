"""Stage 9: synteny and context conservation across a family's occurrences (section 42).

See src/plasmidann/synteny.py for what the six measurements are and why they are kept
apart. This script builds each family's occurrence list from the annotation table and the
protein map, then measures.

It reads the SAME neighbourhood definition as Stage 8 - plasmidann.context.directons and
neighbourhood - rather than restating it. Two definitions of "neighbour" in one pipeline
would diverge the moment one window changed, and the two stages' numbers would stop being
comparable without anything saying so.
"""
import _ctx  # noqa: F401
import collections
import csv

from plasmidann.context import directons, neighbourhood
from plasmidann.synteny import conservation

cfg = snakemake.params.context

# ------------------------------------------------------------------------------------
# Every ORF per plasmid, with the label the cascade gave it.
# ------------------------------------------------------------------------------------
by_plasmid = collections.defaultdict(list)
label_of = {}
with open(snakemake.input.annotation, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        by_plasmid[r["plasmid_id"]].append({
            "orf_id": r["orf_id"], "start": int(r["start"]), "end": int(r["end"]),
            "strand": 1 if r["strand"] in ("1", "+") else -1})
        label_of[r["orf_id"]] = r.get("annot_label") or ""

orfs_of_seq = collections.defaultdict(list)
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, orf_ids = line.rstrip("\n").split("\t")
        for orf_id in orf_ids.split(","):
            orfs_of_seq[sid].append(orf_id)

# ------------------------------------------------------------------------------------
# Ordered neighbours per ORF, computed once per plasmid.
#
# Split into left and right around the ORF rather than returned as one list: section 42's
# left_neighbor_conservation and right_neighbor_conservation are separate measurements and
# cannot be recovered from a combined neighbourhood.
# ------------------------------------------------------------------------------------
context_of = {}
window = cfg["neighbourhood_window"]
for plasmid_id, genes in by_plasmid.items():
    ordered = sorted(genes, key=lambda g: (g["start"], g["end"]))
    positions = {g["orf_id"]: i for i, g in enumerate(ordered)}
    units = directons(genes, max_gap=cfg["max_operon_gap"])
    unit_of = {oid: i for i, unit in enumerate(units) for oid in unit}

    for gene in ordered:
        orf_id = gene["orf_id"]
        i = positions[orf_id]
        # Nearest-first on both sides, so [0] is the immediate neighbour.
        left = [label_of.get(g["orf_id"], "") for g in ordered[max(0, i - window):i]][::-1]
        right = [label_of.get(g["orf_id"], "") for g in ordered[i + 1:i + 1 + window]]

        unit = units[unit_of[orf_id]] if orf_id in unit_of else [orf_id]
        partners = [p for p in unit if p != orf_id]
        context_of[orf_id] = {
            "left": left,
            "right": right,
            # Operon-like means sharing a transcriptional unit with a NAMED partner. An
            # unnamed partner says nothing about what the arrangement is for.
            "operon": any(label_of.get(p) for p in partners),
        }

COLS = ["family_id", "n_occurrences", "context_recurrence",
        "left_neighbor_conservation", "right_neighbor_conservation",
        "neighborhood_conservation", "operon_like_conservation",
        "synteny_conservation", "modal_left", "modal_right", "modal_synteny", "status"]

n_rows = n_conserved = 0
with open(snakemake.input.families, newline="") as fh, \
        open(snakemake.output.tsv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=COLS, delimiter="\t")
    writer.writeheader()

    for family in csv.DictReader(fh, delimiter="\t"):
        occurrences = []
        for member in family["members"].split(","):
            for orf_id in orfs_of_seq.get(member, ()):
                if orf_id in context_of:
                    occurrences.append(context_of[orf_id])

        result = conservation(occurrences)
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
