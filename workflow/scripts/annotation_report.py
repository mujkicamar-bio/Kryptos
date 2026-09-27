"""The deliverable: one complete annotation table, and one complete family table.

WHAT THIS IS FOR

The pipeline's output is not a shortlist. It is every annotation the run could produce, in a
form you can open, sort and filter yourself, so that selecting candidates is YOUR step and
not a decision baked into a Snakemake rule.

Two files, because there are two natural units:

  annotation_complete.csv     one row per ORF - 9.3M of them. What the cascade called it,
                              how much of it that explained, what the dark evidence says,
                              its orthology terms, the plasmid label databases' labels of
                              its protein (amr, metal, ta, conj_role, mge, antidefence), its
                              CONJScan system and its plasmid's mobility class, its
                              plasmid's geNomad phage-plasmid label, the close-level
                              synteny of its close cluster, and - where the ORF is dark -
                              the evidence assembled for its family.
  dark_families_complete.csv  one row per dark family. The selection surface: every piece
                              of evidence the run produced, side by side - annotation,
                              evolution, context, structure, distribution (Stage 7),
                              synteny (Stage 9), the dark sequences its members travel
                              with (S8g), rarity (Stage 14) and the Stage 15 evidence
                              dimensions.

WHAT CONTEXT HOLDS is not a column here: 12_context_and_structure/family_context_terms.tsv
is the long table of context terms per family (amr:, metal:, defence:, conj: ...), counted
over lineages. It carries no top term and no confidence; the family table carries the
descriptive rates (cons_*) only, beside the co-occurrence partners from
dark_cooccurrence.tsv, which are pairs of dark protein sequences rather than context terms.

CSV, not TSV, because these are the files that get opened in a spreadsheet. Every field is
quoted by csv.writer where it needs to be, which matters: a DIAMOND stitle is free text and
routinely contains commas.

NOTHING IS FILTERED AND NOTHING IS RANKED HERE

Every ORF appears, including artefact-flagged ones. Every dark family appears, including
ORPHANs and families that failed every test. Absence of evidence is written as an explicit
status - TOO_FEW_MEMBERS, TOO_SHORT, NO_SIGNAL - never as a blank that reads as a failed
test. Choosing what to do with all that is the report's job, and the report is you.

There is no score, rank or candidate column. evidence_dimension_count counts the distinct
measurements present (plasmidann.integration); supporting_observations_count is reported
beside it and labelled, in the column name itself, as not independent.
"""
import collections
import csv

import _ctx  # noqa: F401

from darkorf import status
from darkorf.ids import family_id as cluster_family_id
from plasmidann import integration, labeldb
from plasmidann.context import overlapping_islands
from plasmidann.cooccurrence import family_partners
from plasmidann.evidence import darkness_state, reality_lines, reality_thresholds

evo_cfg = snakemake.params.evolution
THRESHOLDS = reality_thresholds(evo_cfg)
cooc_cfg = snakemake.params.cooccurrence


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

# S8g: pairs of dark sequences that share a plasmid in more lineages than chance predicts.
# The pair table holds the reported pairs only; a family none of whose members is in one
# had no member together with another sequence in cooccurrence.min_lineages_together
# lineages.
with open(snakemake.input.cooccurrence, newline="") as fh:
    pairs = [{**r, "p_value": float(r["p_value"]), "q_value": float(r["q_value"]),
              "n_lineages_together": int(r["n_lineages_together"]),
              "fraction_of_a": float(r["fraction_of_a"]),
              "fraction_of_b": float(r["fraction_of_b"])}
             for r in csv.DictReader(fh, delimiter="\t")]

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
partners = family_partners(pairs, cooc_cfg["fdr"], family_of_seq)

# ------------------------------------------------------------------------------------
# The family table: the selection surface.
# ------------------------------------------------------------------------------------
FAMILY_COLS = [
    "family_id", "representative", "family_class",
    "n_members", "n_orfs", "n_plasmids", "n_mob_clusters",
    "dark_member_count", "annotated_member_count", "percentage_dark_in_family", "dark_only",
    # small plasmids and large ones; the small_ columns below repeat a
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
    "cons_conj",
    # Stage 7: SEVEN counts, never collapsed. A family on forty copies of one
    # redeposited plasmid is one observation, and reading only the first of these numbers
    # is how a reader concludes otherwise.
    "plasmid_occurrence_count", "unique_plasmid_count",
    "independent_plasmid_cluster_count", "independent_cluster_status",
    "host_count", "genus_count", "n_plasmids_with_host",
    "n_plasmids_with_species", "host_count_status", "n_plasmids_with_predicted_range",
    "predicted_host_range_count",
    "predicted_host_ranges", "MOB_count", "habitat_count",
    "database_source_count",
    # Stage 9: six conservation measurements, kept apart because they fail
    # apart - a conserved left neighbour with a variable right one is a real arrangement
    # that a single averaged context score would hide. Counted over Stage 6 lineages, one
    # vote per lineage; the family table carries the primary-level (family) rows.
    "context_recurrence", "n_occurrences", "n_lineages", "n_lineages_discordant",
    "lineage_left_conservation", "lineage_right_conservation",
    "lineage_neighborhood_conservation", "lineage_operon_like_conservation",
    "lineage_synteny_conservation", "modal_left", "modal_right", "modal_synteny",
    "synteny_status", "synteny_min_lineages",
    "small_n_occurrences", "small_context_recurrence", "small_n_lineages",
    "small_n_lineages_discordant", "small_lineage_left_conservation",
    "small_lineage_right_conservation", "small_lineage_neighborhood_conservation",
    "small_lineage_operon_like_conservation", "small_lineage_synteny_conservation",
    "small_modal_left", "small_modal_right", "small_modal_synteny", "small_synteny_status",
    # S8g: the dark sequences this family's dark members travel with. A partner sequence
    # counts when it shares a plasmid with a member in more lineages than chance predicts
    # (q <= fdr); the best partner (lowest q) is named whether significant or not, with its
    # q and the fraction of the member's lineages in which the two share a plasmid.
    # TOO_FEW_LINEAGES: no member shared a plasmid with another sequence in min_lineages
    # lineages, so no pair of it was reported.
    "cooccurrence_status", "n_cooccurring_partners", "top_cooccurring_partner",
    "top_cooccurring_partner_q", "top_cooccurring_partner_fraction",
    "cooccurrence_fdr", "cooccurrence_min_lineages",
    # Stage 14: descriptors, not a ranking. RARE is not better than
    # WIDESPREAD.
    "rarity_labels",
    # Stage 15: dimensions counted, never scored.
    "evidence_dimensions_present", "evidence_dimension_count",
    "supporting_observations_count", "supporting_observations_are_not_independent",
]

# Stage 7 and Stage 14 both report the distribution counts, because each stage needs them.
# The report takes them from Stage 7, which is where they are computed; taking rarity's
# copies as well would put the same number in the row twice under one name, and whichever
# was merged last would win silently if the two ever disagreed.
RARITY_COLS = ("rarity_labels",)

# S9 writes a bare `status`. Every stage does, which is exactly why it cannot be merged
# under that name: the family row already carries dnds_status and independent_cluster_status
# and a third would overwrite by accident rather than by decision.
SYNTENY_COLS = ("n_occurrences", "context_recurrence", "n_lineages", "n_lineages_discordant",
                "lineage_left_conservation", "lineage_right_conservation",
                "lineage_neighborhood_conservation", "lineage_operon_like_conservation",
                "lineage_synteny_conservation", "modal_left", "modal_right",
                "modal_synteny", "synteny_min_lineages", "small_n_occurrences",
                "small_context_recurrence", "small_n_lineages",
                "small_n_lineages_discordant", "small_lineage_left_conservation",
                "small_lineage_right_conservation",
                "small_lineage_neighborhood_conservation",
                "small_lineage_operon_like_conservation",
                "small_lineage_synteny_conservation", "small_modal_left",
                "small_modal_right", "small_modal_synteny")

# Every small_ measurement the evolution and synteny tables write must reach the report.
# extrasaction="ignore" below drops any column not listed without a word, so the lists are
# checked against the tables' own headers here.
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
                           "independent_cluster_status", "host_count",
                           "genus_count", "n_plasmids_with_host",
                           "n_plasmids_with_species", "host_count_status",
                           "n_plasmids_with_predicted_range",
                           "predicted_host_range_count", "predicted_host_ranges",
                           "MOB_count", "habitat_count",
                           "database_source_count")},
               **{c: syn.get(c, "") for c in SYNTENY_COLS},
               "synteny_status": syn.get("status", ""),
               "small_synteny_status": syn.get("small_status", ""),
               **partners.get(fid, {"n_cooccurring_partners": 0}),
               "cooccurrence_status": (status.SUCCESS if fid in partners
                                       else status.TOO_FEW_LINEAGES),
               "cooccurrence_fdr": cooc_cfg["fdr"],
               "cooccurrence_min_lineages": cooc_cfg["min_lineages_together"],
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
# Per ORF, the close-level (90% identity) synteny row of the ORF's close cluster, where
# Stage 9 measured one - that is, where the cluster holds a dark small-plasmid member.
CLOSE_COLS = {"close_family_id": "family_id", "close_n_lineages": "n_lineages",
              "close_lineage_synteny_conservation": "lineage_synteny_conservation",
              "close_modal_synteny": "modal_synteny", "close_synteny_status": "status"}
close_rows = {fid: r for fid, r in synteny.items() if r.get("level") == "close"}
close_of_seq = {}
with open(snakemake.input.clusters_close) as fh:
    for line in fh:
        rep, member = line.rstrip("\n").split("\t")
        cid = cluster_family_id("close", rep)
        if cid in close_rows:
            close_of_seq[member] = cid

# Per protein, the plasmid label databases' labels (S4d) by term type, as 'source:label':
# the ORF's own labels, in the vocabulary the context terms use (plasmidann.labeldb).
# AMRFinderPlus VIRULENCE and STRESS acid/heat elements have no term type and no column.
LABEL_COLS = {"amr": "amr_labels", "metal": "metal_labels", "ta": "ta_labels",
              "conj_role": "conj_role_labels", "mge": "mge_labels",
              "antidefence": "antidefence_labels"}
labels_of_seq = collections.defaultdict(lambda: collections.defaultdict(set))
with open(snakemake.input.labels_plasmid, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        prefix = labeldb.term_prefix(r)
        if prefix:
            labels_of_seq[r["seq_id"]][LABEL_COLS[prefix]].add(f"{r['source']}:{r['label']}")

# Per ORF, its CONJScan system and component (S8f); per plasmid, its mobility class.
conj_of_orf = collections.defaultdict(lambda: {"conj_system": set(), "conj_component": set()})
with open(snakemake.input.conjugation, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        conj_of_orf[r["orf_id"]]["conj_system"].add(r["system"])
        conj_of_orf[r["orf_id"]]["conj_component"].add(r["component"])
conj_class = {}
with open(snakemake.input.conjugation_class, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        conj_class[r["plasmid_id"]] = r["class"]
# Per plasmid, its phage-plasmid label from geNomad (S8h).
with open(snakemake.input.phage_plasmids, newline="") as fh:
    phage_plasmid = {r["plasmid_id"]: r["phage_plasmid"]
                     for r in csv.DictReader(fh, delimiter="\t")}

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
               "predicted_host_range"]
            + list(LABEL_COLS.values())
            + ["conj_system", "conj_component", "plasmid_conjscan_class",
               "plasmid_phage_plasmid"]
            + list(CLOSE_COLS) + CARRIED)
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
            # Several labels of one type, or several systems, are joined with '; '. Label
            # names contain commas (CARD families), and csv quoting keeps them intact.
            for col, values in labels_of_seq.get(sid, {}).items():
                row[col] = "; ".join(sorted(values))
            for col, values in conj_of_orf.get(r["orf_id"], {}).items():
                row[col] = "; ".join(sorted(values))
            row["plasmid_conjscan_class"] = conj_class.get(r["plasmid_id"], "")
            row["plasmid_phage_plasmid"] = phage_plasmid.get(r["plasmid_id"], "")
            close = close_rows.get(close_of_seq.get(sid, ""))
            if close:
                row.update({c: close[k] for c, k in CLOSE_COLS.items()})
            if fid:
                n_dark += 1
                fam = family_rows[fid]
                row.update({c: fam.get(c, "") for c in CARRIED})
                row["family_id"] = fid
            w.writerow(row)

print(f"annotation_complete.csv: {n_orfs} ORFs ({n_dark} in a dark family)")
print(f"dark_families_complete.csv: {len(families)} families, nothing filtered, "
      "nothing ranked")
