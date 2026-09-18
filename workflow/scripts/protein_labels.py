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

DEDUPLICATION

A label seen several times for one protein - the same Pfam family hit by two tiers, the
same product name from ten nr subjects - is one statement, not ten. Rows are keyed on
(protein_id, source, kind, label) and the best supporting statistics are kept, because an
unmerged table would let a widespread label outvote a rare one purely by copy number when
the categories are counted.
"""
import _ctx  # noqa: F401
import csv
import sys

from plasmidann import labels, pfam_meta

# Which database and version produced each source, for the provenance columns. A label
# without its database release cannot be reproduced, and a category built on it cannot be
# described in a paper.
_DATABASE = {
    "pfam": ("Pfam-A", "pfam_version"),
    "swissprot": ("NCBI swissprot", "swissprot_version"),
    "nr": ("NCBI nr", "nr_version"),
    "eggnog": ("eggNOG", "eggnog_version"),
}

params = snakemake.params
pfam = pfam_meta.load(snakemake.input.pfam_dat)

# (protein_id, source, kind, label) -> row. Identity is fixed on first sight and later
# sightings only improve the statistics, so a label's evidence is the strongest seen rather
# than the last row read.
rows = {}


def _as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("inf")


def add(protein_id, source, tier, entry, evalue="", coverage=""):
    """Admit one label, merging it with any previous sighting of the same statement."""
    if entry["kind"] not in labels.KINDS:
        # An undeclared kind is a statement nothing downstream can group. Failing here
        # names the kind; discovering it later looks like a missing category.
        raise SystemExit(
            f"protein_labels: undeclared label kind {entry['kind']!r} from source "
            f"{source!r}. Add it to plasmidann.labels.KINDS or stop emitting it.")
    database, version_key = _DATABASE[source]
    key = (protein_id, source, entry["kind"], entry["label"])
    existing = rows.get(key)
    if existing is None:
        rows[key] = {
            "protein_id": protein_id, "source": source, "tier": tier,
            "kind": entry["kind"], "label": entry["label"],
            "accession": entry.get("accession", ""),
            "evidence_evalue": evalue, "evidence_coverage": coverage,
            "database": database, "database_version": str(params.get(version_key, "")),
        }
        return
    if _as_float(evalue) < _as_float(existing["evidence_evalue"]):
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
            n_hits += 1
            # The source travels with the row, declared per tier in config/cascade.yaml.
            # It used to be looked up from the tier id here, which broke the moment a tier
            # was inserted. labels_from_hit refuses a row without one.
            for entry in labels.labels_from_hit(row, pfam=pfam):
                add(row["query"], row["source"], row.get("tier", ""), entry,
                    evalue=row.get("evalue", ""), coverage=row.get("coverage", ""))

# ------------------------------------------------------------------------------------
# Orthology: gene symbols and every controlled identifier eggNOG assigns.
# ------------------------------------------------------------------------------------
n_orth = 0
with open(snakemake.input.orthology, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        n_orth += 1
        for entry in labels.labels_from_orthology(row):
            add(row["seq_id"], "eggnog", "S4b", entry)

cols = ["protein_id", "source", "tier", "kind", "label", "accession",
        "evidence_evalue", "evidence_coverage", "database", "database_version"]
with open(snakemake.output.tsv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    writer.writeheader()
    for key in sorted(rows):
        writer.writerow(rows[key])

by_kind = {}
for (_, _, kind, _) in rows:
    by_kind[kind] = by_kind.get(kind, 0) + 1
distinct_labels = len({(k, l) for (_, _, k, l) in rows})

print(f"protein_labels: read {n_hits} hits and {n_orth} orthology rows -> "
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
