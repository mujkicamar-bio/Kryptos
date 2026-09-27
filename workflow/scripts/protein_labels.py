"""One long table of every functional label every tool produced, and their disagreements.

Inputs: the cascade hits (hits.tsv), the orthology table (eggNOG-mapper and PlasmidScope),
the plasmid label databases' table (protein_labels_plasmid.tsv) and the Pfam metadata.
Labels are copied verbatim with their kind and database release (plasmidann.labels); no
biological role is assigned here. One row per (protein_id, source, kind, label, sub_label):
a label seen several times for one protein keeps the row of its best e-value, so a label
is not counted once per supporting hit. The cascade searched one representative per 90%
search cluster; each member takes its representative's labels, named in
via_representative. Spiked controls and decoys contribute no labels. The plasmid label
database rows keep their sub_label and their tier or call (1, 2, Perfect, Strict, the
AMRFinderPlus method); for every other source sub_label is empty.

label_disagreements.tsv lists the cross-source conflicts found by labeldb.disagreements;
it changes no label.
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

# (protein_id, source, kind, label, sub_label) -> row.
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
        existing.update(tier=tier, evidence_evalue=evalue, evidence_coverage=coverage,
                        accession=entry.get("accession", "") or existing["accession"])


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
            add(row["seq_id"], "eggnog", row["orthology_source"], entry)

# ------------------------------------------------------------------------------------
# The plasmid label databases, already one row per statement.
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

# Every named protein carries at least one label, and the cascade names a large fraction
# of any real collection, so an empty table means an input is wrong.
if not rows:
    sys.exit("protein_labels: no labels extracted from any source - the functional "
             "grouping has no substrate. Check that hits.tsv carries a tier column and "
             "that orthology.tsv is not empty.")
