"""The deliverable: one complete annotation table, and one complete family table.

WHAT THIS IS FOR

The pipeline's output is not a shortlist. It is every annotation the run could produce, in a
form you can open, sort and filter yourself, so that selecting candidates is YOUR step and
not a decision baked into a Snakemake rule.

Two files, because there are two natural units:

  annotation_complete.csv     one row per ORF - 9.3M of them. What the cascade called it,
                              how much of it that explained, what the dark evidence says,
                              its orthology terms, and - where the ORF is dark - the
                              evidence assembled for its family.
  dark_families_complete.csv  one row per dark family. The selection surface: every piece
                              of evidence S6-S8 produced, side by side.

CSV, not TSV, because these are the files that get opened in a spreadsheet. Every field is
quoted by csv.writer where it needs to be, which matters: a DIAMOND stitle is free text and
routinely contains commas.

NOTHING IS FILTERED AND NOTHING IS RANKED HERE

Every ORF appears, including artefact-flagged ones. Every dark family appears, including
ORPHANs and families that failed every test. Absence of evidence is written as an explicit
status - TOO_FEW_MEMBERS, NO_DIVERGENCE, NO_SIGNAL - never as a blank that reads as a failed
test. Choosing what to do with all that is the report's job, and the report is you.
"""
import _ctx  # noqa: F401
import collections
import csv

from plasmidann.targets import (reality_lines, darkness_state, reality_thresholds,
                                check_reality_config)

cfg = snakemake.params.prioritisation
evo_cfg = snakemake.params.evolution
check_reality_config(cfg, evo_cfg)
THRESHOLDS = reality_thresholds(cfg, evo_cfg)


def index(path, key):
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r[key]] = r
    return out


families = index(snakemake.input.families, "family_id")
evolution = index(snakemake.input.evolution, "family_id")
recheck = index(snakemake.input.recheck, "family_id")
def index_context(path):
    """Reduce the LONG family_context table to one row per family.

    S8c writes one row per (family, category) because the label vocabulary is open. The
    report needs a per-family view, and it takes the strongest association by q-value -
    the only one of the three numbers that distinguishes a category conserved across three
    plasmids from one conserved across three hundred.

    Every category the family carries is also kept as cons_<category>, so a reader can
    still ask for one by name.
    """
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            row = out.setdefault(r["family_id"], {"family_id": r["family_id"]})
            row[f"cons_{r['category']}"] = r["observed_rate"]
            if r["status"] != "SUCCESS":
                continue
            try:
                q = float(r["q_value"])
            except (TypeError, ValueError):
                continue
            if "_best_q" not in row or q < row["_best_q"]:
                row.update({"_best_q": q, "top_hypothesis": r["category"],
                            "top_subcategories": r["subcategories"],
                            "top_conservation": r["observed_rate"],
                            "top_enrichment": r["enrichment"],
                            "top_q_value": r["q_value"]})
    return out


context = index_context(snakemake.input.context)
structure = index(snakemake.input.structure, "seq_id")
orthology = index(snakemake.input.orthology, "seq_id")

# Which unique protein each ORF is, and which family each unique protein belongs to.
seq_of_orf = {}
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, ids = line.rstrip("\n").split("\t")
        for oid in ids.split(","):
            seq_of_orf[oid] = sid

family_of_seq = {}
for fid, fam in families.items():
    for member in fam["members"].split(","):
        family_of_seq[member] = fid

# ------------------------------------------------------------------------------------
# The family table: the selection surface.
# ------------------------------------------------------------------------------------
FAMILY_COLS = [
    "family_id", "representative", "family_class",
    "n_members", "n_orfs", "n_plasmids", "n_mob_clusters",
    # evidence that it is a real protein
    "reality_n", "reality_lines", "reality_lines_implied",
    "dnds_median", "dnds_min", "dnds_status", "under_purifying_selection", "n_pairs",
    "rnacode_p", "rnacode_p_antisense", "rnacode_status", "coding_signal",
    "consensus_hit", "consensus_label", "collectively_novel",
    # what it might do
    "darkness_state", "structural_match", "structural_description", "structure_evalue",
    "top_hypothesis", "top_subcategories", "top_conservation", "top_enrichment",
    "top_q_value", "high_confidence",
    # The named context columns that survive the curated list's deletion. defence and
    # integron are ISLAND membership, called by DefenseFinder and IntegronFinder, and are
    # categories in their own right whatever the label grouping does. backbone_adjacent
    # and ta_candidate are gone with the 73-name list that defined them; the equivalent
    # claim now lives in top_hypothesis, which names whatever the tools actually called
    # the neighbour.
    "cons_defence", "cons_integron",
    "cons_annotated_neighbour", "cons_operon_with_annotated", "cons_two_gene_operon",
]

family_rows = {}
with open(snakemake.output.families, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=FAMILY_COLS, extrasaction="ignore")
    w.writeheader()
    for fid, fam in families.items():
        evo = evolution.get(fid, {})
        ctx = context.get(fid, {})
        struct = structure.get(fam["representative"], {})
        record = {"dnds_status": evo.get("dnds_status", ""),
                  "dnds_median": evo.get("dnds_median", ""),
                  "n_mob_clusters": fam["n_mob_clusters"],
                  "n_members": fam["n_members"],
                  "structural_match": struct.get("target", "")}
        n, fired, implied = reality_lines(record, THRESHOLDS)
        row = {**fam, **evo, **recheck.get(fid, {}), **ctx,
               "reality_n": n,
               "reality_lines": "+".join(fired) or "none",
               "reality_lines_implied": "+".join(implied) or "none",
               "darkness_state": darkness_state(record),
               "structural_match": struct.get("target", ""),
               "structural_description": struct.get("target_description", ""),
               "structure_evalue": struct.get("evalue", "")}
        family_rows[fid] = row
        w.writerow(row)

# ------------------------------------------------------------------------------------
# The ORF table: everything, with the family evidence joined on where it exists.
# ------------------------------------------------------------------------------------
CARRIED = ["family_id", "reality_n", "reality_lines", "darkness_state", "dnds_median",
           "dnds_status", "coding_signal", "collectively_novel",
           "top_hypothesis", "top_conservation", "structural_match",
           "structural_description"]

n_orfs = n_dark = 0
with open(snakemake.input.annotation, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    # The ORF table has no family columns of its own, so CARRIED is joined on by name.
    # Any collision would silently overwrite an annotation column, so it is refused.
    clash = set(CARRIED) & set(reader.fieldnames)
    if clash:
        raise SystemExit(
            f"the family columns {sorted(clash)} already exist in the annotation table; "
            "joining them would overwrite annotation with family evidence")
    cols = (list(reader.fieldnames)
            + ["seq_id", "cog_category", "kegg_pathways", "preferred_name",
               "eggnog_description"] + CARRIED)
    with open(snakemake.output.annotation, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in reader:
            n_orfs += 1
            sid = seq_of_orf.get(r["orf_id"], "")
            fid = family_of_seq.get(sid)
            orth = orthology.get(sid, {})
            row = {**r, "seq_id": sid,
                   "cog_category": orth.get("cog_category", ""),
                   "kegg_pathways": orth.get("kegg_pathways", ""),
                   "preferred_name": orth.get("preferred_name", ""),
                   "eggnog_description": orth.get("eggnog_description", "")}
            if fid:
                n_dark += 1
                fam = family_rows[fid]
                row.update({c: fam.get(c, "") for c in CARRIED})
                row["family_id"] = fid
            w.writerow(row)

print(f"annotation_complete.csv: {n_orfs} ORFs ({n_dark} in a dark family)")
print(f"dark_families_complete.csv: {len(families)} families, nothing filtered, "
      "nothing ranked")
