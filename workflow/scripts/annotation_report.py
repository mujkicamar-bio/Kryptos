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
                              of evidence the run produced, side by side - annotation,
                              evolution, context, structure, distribution (Stage 7),
                              synteny (Stage 9), rarity (Stage 14) and the Stage 15
                              evidence dimensions.

CSV, not TSV, because these are the files that get opened in a spreadsheet. Every field is
quoted by csv.writer where it needs to be, which matters: a DIAMOND stitle is free text and
routinely contains commas.

NOTHING IS FILTERED AND NOTHING IS RANKED HERE

Every ORF appears, including artefact-flagged ones. Every dark family appears, including
ORPHANs and families that failed every test. Absence of evidence is written as an explicit
status - TOO_FEW_MEMBERS, NO_DIVERGENCE, NO_SIGNAL - never as a blank that reads as a failed
test. Choosing what to do with all that is the report's job, and the report is you.

Spec section 76 draws the boundary as a list of columns that must NOT be here:
candidate_score, novelty_score, experimental_rank, top_1000. Section 2.3 gives the reason,
and it is not tidiness - two proteins with the same composite can be entirely different
bets, one with overwhelming evidence that it is a real protein and no idea what it does,
the other with a sharp hypothesis resting on almost nothing. Those demand different
experiments, and a single number destroys the distinction.

evidence_dimension_count counts DISTINCT MEASUREMENTS present, which is why it belongs
here where a score does not. supporting_observations_count is reported beside it and
labelled, in the column name itself, as not independent.
"""
import _ctx  # noqa: F401
import collections
import csv

from plasmidann import integration
from plasmidann.evidence import reality_lines, darkness_state, reality_thresholds

# check_reality_config went with Layer C: it validated that min_reality_lines - a SELECTION
# parameter - was reachable. Nothing selects here, so there is no such parameter to check.
evo_cfg = snakemake.params.evolution
THRESHOLDS = reality_thresholds(evo_cfg)


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
recurrence = index(snakemake.input.recurrence, "family_id")
synteny = index(snakemake.input.synteny, "family_id")
rarity = index(snakemake.input.rarity, "family_id")

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
    "top_q_value",
    # The named context columns that survive the curated list's deletion. defence and
    # integron are ISLAND membership, called by DefenseFinder and IntegronFinder, and are
    # categories in their own right whatever the label grouping does. backbone_adjacent
    # and ta_candidate are gone with the 73-name list that defined them; the equivalent
    # claim now lives in top_hypothesis, which names whatever the tools actually called
    # the neighbour.
    "cons_defence", "cons_integron",
    "cons_annotated_neighbour", "cons_operon_with_annotated", "cons_two_gene_operon",
    # Stage 7 (section 34): SEVEN counts, never collapsed. A family on forty copies of one
    # redeposited plasmid is one observation, and reading only the first of these numbers
    # is how a reader concludes otherwise.
    "plasmid_occurrence_count", "unique_plasmid_count",
    "independent_plasmid_cluster_count", "independent_cluster_status",
    "host_count", "species_count", "genus_count", "MOB_count", "habitat_count",
    "database_record_count", "database_source_count",
    # Stage 9 (section 42): six conservation measurements, kept apart because they fail
    # apart - a conserved left neighbour with a variable right one is a real arrangement
    # that a single averaged context score would hide.
    "context_recurrence", "n_occurrences", "left_neighbor_conservation",
    "right_neighbor_conservation", "neighborhood_conservation",
    "operon_like_conservation", "synteny_conservation", "modal_left", "modal_right",
    "modal_synteny", "synteny_status",
    # Stage 14 (section 54): descriptors, not a ranking. RARE is not better than
    # WIDELY_CONSERVED. The version travels because a label's definition can change.
    "rarity_labels", "rarity_version",
    # Stage 15 (sections 56-57): dimensions counted, never scored.
    "evidence_dimensions_present", "evidence_dimension_count",
    "supporting_observations_count", "supporting_observations_are_not_independent",
    "functional_hypothesis", "functional_hypothesis_support",
]

# Stage 7 and Stage 14 both report the distribution counts, because each stage needs them.
# The report takes them from Stage 7, which is where they are computed; taking rarity's
# copies as well would put the same number in the row twice under one name, and whichever
# was merged last would win silently if the two ever disagreed.
RARITY_COLS = ("rarity_labels", "rarity_version")

# S9 writes a bare `status`. Every stage does, which is exactly why it cannot be merged
# under that name: the family row already carries dnds_status and independent_cluster_status
# and a third would overwrite by accident rather than by decision.
SYNTENY_COLS = ("n_occurrences", "context_recurrence", "left_neighbor_conservation",
                "right_neighbor_conservation", "neighborhood_conservation",
                "operon_like_conservation", "synteny_conservation", "modal_left",
                "modal_right", "modal_synteny")

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
        rec = recurrence.get(fid, {})
        syn = synteny.get(fid, {})
        rar = rarity.get(fid, {})
        row = {**fam, **evo, **recheck.get(fid, {}), **ctx,
               # Explicit column lists rather than a dict merge: recurrence carries its own
               # `representative` and synteny its own `status`, and merging wholesale would
               # overwrite the family's representative and one of the row's other statuses
               # without anything failing.
               **{c: rec.get(c, "") for c in FAMILY_COLS
                  if c in ("plasmid_occurrence_count", "unique_plasmid_count",
                           "independent_plasmid_cluster_count",
                           "independent_cluster_status", "host_count", "species_count",
                           "genus_count", "MOB_count", "habitat_count",
                           "database_record_count", "database_source_count")},
               **{c: syn.get(c, "") for c in SYNTENY_COLS},
               "synteny_status": syn.get("status", ""),
               **{c: rar.get(c, "") for c in RARITY_COLS},
               "reality_n": n,
               "reality_lines": "+".join(fired) or "none",
               "reality_lines_implied": "+".join(implied) or "none",
               "darkness_state": darkness_state(record),
               "structural_match": struct.get("target", ""),
               "structural_description": struct.get("target_description", ""),
               "structure_evalue": struct.get("evalue", "")}

        # Stage 15 reads the assembled row, so it sees exactly the evidence a reader sees.
        # Computing it from the source tables instead would let the two drift, and the
        # dimension count is a claim ABOUT this row.
        row.update(integration.evidence_summary(row))
        hypothesis, support = integration.functional_hypothesis(row)
        row["functional_hypothesis"] = hypothesis
        row["functional_hypothesis_support"] = support

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
