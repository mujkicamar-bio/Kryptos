"""S4c: one long table of every functional label every tool produced.

WHY THIS TABLE EXISTS

The question it serves is "is this protein replication, mobilisation, or conjugation", and
the answer cannot come from a hand-written list. The list this replaces held 73 Pfam family
names; Pfam-A 38.2 holds 30,134 families, of which 67 mention replication in their
description and 42 mention conjugation. It named 16 and 15, named no MobB and no MobD,
included nine names that do not exist in Pfam-A at all, and assigned every role with no
source.

So the labels are collected verbatim, with their kind and their provenance, and the
grouping into biological categories is derived afterwards from the vocabulary observed
here. That order matters: a category built from the labels the data actually contains can be
described in a methods section, and a category built from recollection cannot.

WHY IT IS LONG AND NOT WIDE

The vocabulary is open. Pfam-A 38.2 has 30,134 families, nr product names are unbounded,
and a protein carries a different number of labels from every source. One row per
(protein, source, kind, label) is the only shape that holds that without a column per
family, and it is the shape a grouping step reads naturally: select the distinct labels of
one kind, decide their categories, join back.

SEARCH-CLUSTER MEMBERS, CONTROLS AND DECOYS

The cascade searched only the representative of each 90% search cluster (S2s), so hits.tsv
holds representatives alone. Every label a representative's hits produced is written for
each of its members too, with the representative in `via_representative`: without that,
about 140,000 members of the selected families carried no Pfam, pharokka, Swiss-Prot or nr
label at all. The spiked controls (CTRL_) and decoys (DECOY_) are instrumentation, not
plasmid proteins, and contribute no labels.

THE PLASMID LABEL DATABASES (S4d)

label_databases.py searches TADB, BacMet, oriTDB, CARD, mobileOG-db, dbAPIS, Anti-CRISPRdb
and AMRFinderPlus and writes 08_protein_labels/protein_labels_plasmid.tsv; its rows are
merged here as their own kinds (plasmidann.labeldb.KIND), with seq_id as protein_id,
label_kind as kind, subject as accession and the database's tier or call (1, 2, Perfect,
Strict, the AMRFinderPlus method) as tier. sub_label is kept - AMRFinderPlus's element
type/subtype decides whether a gene is an amr or a metal context term - and is empty for
every other source. Identity and coverage stay in protein_labels_plasmid.tsv, whose
columns differ from the cascade's e-value and query coverage. The source 'card' names
both pharokka's CARD search (kinds card_gene_family, card_mechanism) and the direct search
of the CARD protein homolog models (kind card_amr_family); the kind tells them apart.

This step also writes 08_protein_labels/label_disagreements.tsv (labeldb.disagreements):
every protein on which two sources make incompatible statements - the Tier 0 gene symbol
against a database naming a gene, CARD against AMRFinderPlus, BacMet against AMRFinderPlus
and CARD, TADB against a DefenseFinder component, oriTDB against a CONJScan component. The
Tier 0 symbols are the KEGG gene symbols of PlasmidScope's KOs; the KEGG KO list is used for
that mapping only. The table changes no label.

DEDUPLICATION

A label seen several times for one protein - the same Pfam family hit by two tiers, the
same product name from ten nr subjects - is one statement, not ten. Rows are keyed on
(protein_id, source, kind, label, sub_label) and the best supporting statistics are kept, because an
unmerged table would let a widespread label outvote a rare one purely by copy number when
the categories are counted.
"""
import csv
import sys

import _ctx  # noqa: F401

from plasmidann import labeldb, labels, pfam_meta
from plasmidann.cascade import as_float
from plasmidann.decoys import DECOY_PREFIX

CONTROL_PREFIX = "CTRL_"

# Which database and version produced each source, for the provenance columns. A label
# without its database release cannot be reproduced, and a category built on it cannot be
# described in a paper.
_DATABASE = {
    "pfam": ("Pfam-A", "pfam_version"),
    "swissprot": ("NCBI swissprot", "swissprot_version"),
    "nr": ("NCBI ClusteredNR", "nr_version"),
    "eggnog": ("eggNOG", "eggnog_version"),
    # pharokka ships the phage families, CARD and VFDB as ONE versioned bundle
    # (data/refs/pharokka/VERSION_x_y_z) and does not expose the CARD or VFDB snapshot
    # dates separately, so all three cite the bundle version. The family table is PHROG v4.
    "pharokka": ("pharokka databases (PHROG v4)", "pharokka_db_version"),
    "card": ("pharokka databases (CARD)", "pharokka_db_version"),
    "vfdb": ("pharokka databases (VFDB)", "pharokka_db_version"),
}

# The database each plasmid label source searched; its release travels on every row.
_LABEL_DATABASE = {
    "tadb": "TADB", "bacmet": "BacMet", "oritdb": "oriTDB",
    "card": "CARD protein homolog models", "mobileog": "mobileOG-db", "dbapis": "dbAPIS",
    "acrdb": "Anti-CRISPRdb", "amrfinder": "AMRFinderPlus database",
}

params = snakemake.params
pfam = pfam_meta.load(snakemake.input.pfam_dat)

# (protein_id, source, kind, label) -> row. Identity is fixed on first sight and later
# sightings only improve the statistics, so a label's evidence is the strongest seen rather
# than the last row read.
rows = {}

# Representative -> the members that take its result (selection.tsv, role 'member').
members_of = {}
with open(snakemake.input.selection, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["role"] == "member":
            members_of.setdefault(r["search_representative"], []).append(r["seq_id"])


def add(protein_id, source, tier, entry, evalue="", coverage="", representative=""):
    """Admit one label, merging it with any previous sighting of the same statement."""
    if entry["kind"] not in labels.KINDS:
        # An undeclared kind is a statement nothing downstream can group. Failing here
        # names the kind; discovering it later looks like a missing category.
        raise SystemExit(
            f"protein_labels: undeclared label kind {entry['kind']!r} from source "
            f"{source!r}. Add it to plasmidann.labels.KINDS or stop emitting it.")
    database, version_key = _DATABASE[source]
    key = (protein_id, source, entry["kind"], entry["label"], "")
    existing = rows.get(key)
    if existing is None:
        rows[key] = {
            "protein_id": protein_id, "source": source, "tier": tier,
            "kind": entry["kind"], "label": entry["label"], "sub_label": "",
            "accession": entry.get("accession", ""),
            "evidence_evalue": evalue, "evidence_coverage": coverage,
            "database": database, "database_version": str(params.get(version_key, "")),
            "via_representative": representative,
        }
        return
    if as_float(evalue) < as_float(existing["evidence_evalue"]):
        existing["evidence_evalue"] = evalue
        existing["evidence_coverage"] = coverage
        existing["accession"] = entry.get("accession", "") or existing["accession"]


# ------------------------------------------------------------------------------------
# Cascade hits: Pfam families from the hmmer tiers, product names from the DIAMOND tiers.
# ------------------------------------------------------------------------------------
n_hits = 0
for path in snakemake.input.hits:
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            q = row["query"]
            if q.startswith(CONTROL_PREFIX) or q.startswith(DECOY_PREFIX):
                continue
            n_hits += 1
            # The source travels with the row, declared per tier in config/cascade.yaml.
            # It used to be looked up from the tier id here, which broke the moment a tier
            # was inserted. labels_from_hit refuses a row without one.
            for entry in labels.labels_from_hit(row, pfam=pfam):
                for pid, via in [(q, "")] + [(m, q) for m in members_of.get(q, ())]:
                    add(pid, row["source"], row.get("tier", ""), entry,
                        evalue=row.get("evalue", ""), coverage=row.get("coverage", ""),
                        representative=via)

# ------------------------------------------------------------------------------------
# Orthology: gene symbols and every controlled identifier eggNOG assigns.
# ------------------------------------------------------------------------------------
n_orth = 0
with open(snakemake.input.orthology, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        n_orth += 1
        for entry in labels.labels_from_orthology(row):
            add(row["seq_id"], "eggnog", "S4b", entry)

# ------------------------------------------------------------------------------------
# The plasmid label databases (S4d), already one row per statement.
# ------------------------------------------------------------------------------------
with open(snakemake.input.labels_plasmid, newline="") as fh:
    plasmid_labels = list(csv.DictReader(fh, delimiter="\t"))
for r in plasmid_labels:
    if r["label_kind"] not in labels.KINDS:
        raise SystemExit(
            f"protein_labels: undeclared label kind {r['label_kind']!r} from source "
            f"{r['source']!r} in {snakemake.input.labels_plasmid}.")
    rows[(r["seq_id"], r["source"], r["label_kind"], r["label"], r["sub_label"])] = {
        "protein_id": r["seq_id"], "source": r["source"], "tier": r["tier"],
        "kind": r["label_kind"], "label": r["label"], "sub_label": r["sub_label"],
        "accession": r["subject"], "evidence_evalue": "", "evidence_coverage": "",
        "database": _LABEL_DATABASE[r["source"]],
        "database_version": r["database_version"], "via_representative": ""}

cols = ["protein_id", "source", "tier", "kind", "label", "sub_label", "accession",
        "evidence_evalue", "evidence_coverage", "database", "database_version",
        "via_representative"]
with open(snakemake.output.tsv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    writer.writeheader()
    for key in sorted(rows):
        writer.writerow(rows[key])

# ------------------------------------------------------------------------------------
# Cross-source disagreements. Written beside the labels; no label changes because of them.
# ------------------------------------------------------------------------------------
def read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


with open(snakemake.input.orthology, newline="") as fh:
    tier0 = labeldb.tier0_symbols(csv.DictReader(fh, delimiter="\t"),
                                  labeldb.read_ko_symbols(snakemake.input.ko_list))
orf_to_protein = labeldb.read_protein_map(snakemake.input.map)
conflicts = labeldb.disagreements(
    plasmid_labels, tier0=tier0,
    defence=labeldb.by_protein(read_rows(snakemake.input.defence), orf_to_protein),
    conj=labeldb.by_protein(read_rows(snakemake.input.conjugation), orf_to_protein))
with open(snakemake.output.disagreements, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=labeldb.DISAGREEMENT_COLUMNS, delimiter="\t")
    writer.writeheader()
    writer.writerows(conflicts)

by_kind = {}
for (_, _, kind, _, _) in rows:
    by_kind[kind] = by_kind.get(kind, 0) + 1
distinct_labels = len({(k, l) for (_, _, k, l, _) in rows})

by_conflict = {}
for c in conflicts:
    by_conflict[c["conflict_type"]] = by_conflict.get(c["conflict_type"], 0) + 1
print(f"protein_labels: {len(conflicts)} cross-source disagreements on "
      f"{len({c['seq_id'] for c in conflicts})} proteins {by_conflict}")
print(f"protein_labels: read {n_hits} hits, {n_orth} orthology rows and "
      f"{len(plasmid_labels)} plasmid label database rows -> "
      f"{len(rows)} label rows on {len({k[0] for k in rows})} proteins, "
      f"{distinct_labels} distinct (kind, label) pairs")
for kind in sorted(by_kind):
    print(f"  {kind:<22} {by_kind[kind]}")

# An empty table here is always a bug: the cascade named a large fraction of the
# collection, and every named protein carries at least one label by construction.
if not rows:
    sys.exit("protein_labels: no labels extracted from any source - the functional "
             "grouping has no substrate. Check that hits.tsv carries a tier column and "
             "that orthology.tsv is not empty.")
