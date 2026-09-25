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
import collections
import csv

import _ctx  # noqa: F401

from plasmidann import integration
from plasmidann.context import overlapping_islands
from plasmidann.evidence import darkness_state, reality_lines, reality_thresholds

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
context = index(snakemake.input.context, "family_id")

# IS elements (S8e) as intervals per plasmid, for the per-ORF is_element column. The
# family view already carries them through the context table as cons_is_element.
is_elements = collections.defaultdict(list)
with open(snakemake.input.is_elements, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        is_elements[r["plasmid_id"]].append(
            {"start": int(r["start"]), "end": int(r["end"]), "family": r["family"]})
structure = index(snakemake.input.structure, "seq_id")
# Per ORF: the observed host of the ORF's own plasmid, and MOB-suite's predicted host range
# for it as a separate measurement (clonal_registry).
registry = index(snakemake.input.registry, "plasmid_id")
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
    # section 31.2
    "dark_member_count", "annotated_member_count", "percentage_dark_in_family", "dark_only",
    # small plasmids and large ones (spec section 13.3); the small_ columns below repeat a
    # measurement over the small-plasmid members or occurrences alone
    "n_small_members", "n_large_members", "scope", "known_from",
    # evidence that it is a real protein
    "reality_n", "reality_lines", "reality_lines_implied",
    "dnds_median", "dnds_min", "dnds_status", "under_purifying_selection", "n_pairs",
    "rnacode_p", "rnacode_p_antisense", "rnacode_status", "coding_signal",
    "small_n_aligned", "small_dnds_median", "small_dnds_min", "small_dnds_status",
    "small_under_purifying_selection", "small_n_pairs", "small_rnacode_p",
    "small_rnacode_p_antisense", "small_rnacode_status", "small_coding_signal",
    "consensus_status", "consensus_hit", "consensus_label", "collectively_novel",
    # what it might do
    "darkness_state", "structural_match", "structural_description", "structure_evalue",
    # S8c: the fraction of the family's plasmids on which a member has each context
    # feature. Descriptive rates; there is no enrichment test.
    "cons_defence", "cons_integron", "cons_is_element",
    "cons_annotated_neighbour", "cons_operon_with_annotated", "cons_two_gene_operon",
    # Stage 7 (section 34): SEVEN counts, never collapsed. A family on forty copies of one
    # redeposited plasmid is one observation, and reading only the first of these numbers
    # is how a reader concludes otherwise.
    "plasmid_occurrence_count", "unique_plasmid_count",
    "independent_plasmid_cluster_count", "independent_cluster_status",
    "host_count", "species_count", "genus_count", "n_plasmids_with_host",
    "n_plasmids_with_species", "host_count_status", "n_plasmids_with_predicted_range",
    "predicted_host_range_count",
    "predicted_host_ranges", "MOB_count", "habitat_count",
    "database_record_count", "database_source_count",
    # Stage 9 (section 42): six conservation measurements, kept apart because they fail
    # apart - a conserved left neighbour with a variable right one is a real arrangement
    # that a single averaged context score would hide.
    "context_recurrence", "n_occurrences", "left_neighbor_conservation",
    "right_neighbor_conservation", "neighborhood_conservation",
    "operon_like_conservation", "synteny_conservation", "modal_left", "modal_right",
    "modal_synteny", "synteny_status",
    "small_n_occurrences", "small_context_recurrence", "small_left_neighbor_conservation",
    "small_right_neighbor_conservation", "small_neighborhood_conservation",
    "small_operon_like_conservation", "small_synteny_conservation", "small_modal_left",
    "small_modal_right", "small_modal_synteny", "small_synteny_status",
    # Stage 14 (section 54): descriptors, not a ranking. RARE is not better than
    # WIDELY_CONSERVED. The version travels because a label's definition can change.
    "rarity_labels", "rarity_version",
    # Stage 15 (sections 56-57): dimensions counted, never scored.
    "evidence_dimensions_present", "evidence_dimension_count",
    "supporting_observations_count", "supporting_observations_are_not_independent",
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
                "modal_right", "modal_synteny", "small_n_occurrences",
                "small_context_recurrence", "small_left_neighbor_conservation",
                "small_right_neighbor_conservation", "small_neighborhood_conservation",
                "small_operon_like_conservation", "small_synteny_conservation",
                "small_modal_left", "small_modal_right", "small_modal_synteny")

# Every small_ measurement the evolution and synteny tables write must reach the report
# (spec section 13.3: each is reported twice). extrasaction="ignore" below drops anything
# not listed without a word - which is how six small_ evolution and five small_ synteny
# columns went missing - so the lists are checked against the tables' own headers here.
for path in (snakemake.input.evolution, snakemake.input.synteny):
    with open(path, newline="") as fh:
        header = next(csv.reader(fh, delimiter="\t"), [])
    missing = [c for c in header if c.startswith("small_") and c != "small_status"
               and c != "small_evidence_note" and c not in FAMILY_COLS]
    if missing:
        raise SystemExit(f"annotation_report: {path} writes {missing}, which the family "
                         "report does not carry. Add them to FAMILY_COLS.")

family_rows = {}
with open(snakemake.output.families, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=FAMILY_COLS, extrasaction="ignore")
    w.writeheader()
    for fid, fam in families.items():
        evo = evolution.get(fid, {})
        ctx = context.get(fid, {})
        struct = structure.get(fam["representative"], {})
        rec = recurrence.get(fid, {})
        record = {"dnds_status": evo.get("dnds_status", ""),
                  "dnds_median": evo.get("dnds_median", ""),
                  "independent_plasmid_cluster_count":
                      rec.get("independent_plasmid_cluster_count", ""),
                  "n_members": fam["n_members"],
                  "structural_match": struct.get("target", "")}
        n, fired, implied = reality_lines(record, THRESHOLDS)
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
                           "genus_count", "n_plasmids_with_host",
                           "n_plasmids_with_species", "host_count_status",
                           "n_plasmids_with_predicted_range",
                           "predicted_host_range_count", "predicted_host_ranges",
                           "MOB_count", "habitat_count",
                           "database_record_count", "database_source_count")},
               **{c: syn.get(c, "") for c in SYNTENY_COLS},
               "synteny_status": syn.get("status", ""),
               "small_synteny_status": syn.get("small_status", ""),
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

        family_rows[fid] = row
        w.writerow(row)

# ------------------------------------------------------------------------------------
# The ORF table: everything, with the family evidence joined on where it exists.
# ------------------------------------------------------------------------------------
CARRIED = ["family_id", "scope", "reality_n", "reality_lines", "darkness_state",
           "dnds_median",
           "dnds_status", "coding_signal", "collectively_novel",
           "structural_match", "structural_description"]

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
               "eggnog_description", "is_element", "host_species", "host_genus",
               "predicted_host_range"] + CARRIED)
    with open(snakemake.output.annotation, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in reader:
            n_orfs += 1
            sid = seq_of_orf.get(r["orf_id"], "")
            fid = family_of_seq.get(sid)
            orth = orthology.get(sid, {})
            host = registry.get(r["plasmid_id"], {})
            row = {**r, "seq_id": sid,
                   "host_species": host.get("species", ""),
                   "host_genus": host.get("genus", ""),
                   "predicted_host_range": host.get("predicted_host_range", ""),
                   "cog_category": orth.get("cog_category", ""),
                   "kegg_pathways": orth.get("kegg_pathways", ""),
                   "preferred_name": orth.get("preferred_name", ""),
                   "eggnog_description": orth.get("eggnog_description", ""),
                   # The IS families of every element the ORF overlaps; empty when none.
                   "is_element": ",".join(dict.fromkeys(
                       e["family"] for e in overlapping_islands(
                           {"start": int(r["start"]), "end": int(r["end"])},
                           is_elements.get(r["plasmid_id"], []))))}
            if fid:
                n_dark += 1
                fam = family_rows[fid]
                row.update({c: fam.get(c, "") for c in CARRIED})
                row["family_id"] = fid
            w.writerow(row)

print(f"annotation_complete.csv: {n_orfs} ORFs ({n_dark} in a dark family)")
print(f"dark_families_complete.csv: {len(families)} families, nothing filtered, "
      "nothing ranked")
