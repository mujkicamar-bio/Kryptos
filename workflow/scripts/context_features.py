"""S8c: genomic context per ORF, aggregated to one row of rates per family.

WHAT IT MEASURES

For every dark family, the fraction of the plasmids carrying it on which a member ORF

    cons_defence                 overlaps a component gene of a DefenseFinder system
    cons_conj                    overlaps a component gene of a CONJScan system
    cons_integron                overlaps an integron cassette array (IntegronFinder)
    cons_is_element              overlaps an IS element (ISEScan)
    cons_annotated_neighbour     has a FUNCTIONAL gene within the +-window neighbourhood
    cons_operon_with_annotated   shares a directon with a FUNCTIONAL gene
    cons_two_gene_operon         is one of exactly two genes in its directon

These are descriptive rates, with no background and no enrichment test. A rate of 0 means
the family was examined and the feature was absent. cons_defence and cons_conj are empty
when that stage recorded NOT_RUN, because nothing was searched. The rates are not corrected
for plasmid size: on a six-gene plasmid a +-3 neighbourhood is the whole molecule, so
cons_annotated_neighbour is high there by construction.

THE UNIT IS THE PLASMID

Not the family member. Members of a family are homologs on plasmids that are frequently
near-identical, so counting members makes sequencing effort look like evidence. Clonal
redundancy between distinct plasmids remains uncorrected; see docs/PARAMETER_PROVENANCE.md.

CONTEXT TERMS (family_context_terms.tsv)

A second, long table names WHAT the context holds: one row per family and term, such as
amr:<CARD family> or defence:<system>, counted per Stage 6 lineage rather than per plasmid,
so clonal redundancy is corrected there. The terms, the two neighbour rules and the
paralogue exclusion are in plasmidann.context_terms. Rows are written for the dark
families above, over their dark members as the rates are, AND for every other family at
the primary resolution over all its members (family_set 'known'): the known families are
the benchmark tools/calibrate_context.py measures precision on.
"""
import collections
import csv
import sys

import _ctx  # noqa: F401

from darkorf import status
from darkorf.circular import is_circular
from plasmidann.context import directons, flanks, overlapping_islands
from plasmidann.context_terms import (
    COLUMNS,
    family_term_rows,
    label_term,
    orf_term_sources,
    system_term,
    window_covers_plasmid,
)

cfg = snakemake.params.context
FEATURES = ("defence", "integron", "is_element", "annotated_neighbour",
            "operon_with_annotated", "two_gene_operon", "conj")
# A family's members column lists every member protein; at the intermediate resolution a
# large family exceeds csv's default field limit of 131,072 characters (~4,000 ids).
csv.field_size_limit(sys.maxsize)

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
# What the context terms are made of: each protein's labels, and each ORF's systems.
# ------------------------------------------------------------------------------------
terms_of_seq = collections.defaultdict(set)
with open(snakemake.input.labels, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        term = label_term(r["kind"], r["label"], r.get("sub_label", ""))
        if term:
            terms_of_seq[r["protein_id"]].add(term)

system_terms = collections.defaultdict(set)
defence_orfs, conj_orfs = set(), set()
not_run = set()
for path, prefix, orfs in ((snakemake.input.defence, "defence", defence_orfs),
                           (snakemake.input.conjugation, "conj", conj_orfs)):
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            # A NOT_RUN stage writes one row with that status and no ORF.
            if r["status"] == status.NOT_RUN:
                not_run.add(prefix)
            else:
                orfs.add(r["orf_id"])
                system_terms[r["orf_id"]].add(system_term(prefix, r["system"]))

# ------------------------------------------------------------------------------------
# Which ORFs correspond to which unique protein, indexed by protein so that each family
# looks up only its own members.
# ------------------------------------------------------------------------------------
orfs_of_seq = collections.defaultdict(list)
# Only ORFs that can be a term's source need their protein, to recognise a tandem
# paralogue: a labelled ORF or a system component. A full reverse map would be 9.3M entries.
label_terms, seq_of_source = {}, {}
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, ids = line.rstrip("\n").split("\t")
        for oid in ids.split(","):
            orfs_of_seq[sid].append(oid)
            if sid in terms_of_seq:
                label_terms[oid] = terms_of_seq[sid]
            if sid in terms_of_seq or oid in system_terms:
                seq_of_source[oid] = sid

# ------------------------------------------------------------------------------------
# Islands, as intervals per plasmid: integron elements, IS elements, and the component
# genes of defence and conjugation systems.
# ------------------------------------------------------------------------------------
islands = collections.defaultdict(list)
with open(snakemake.input.integrons, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        islands[r["plasmid_id"]].append(
            {"name": "integron", "start": int(r["start"]), "end": int(r["end"])})

with open(snakemake.input.is_elements, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        islands[r["plasmid_id"]].append(
            {"name": "is_element", "start": int(r["start"]), "end": int(r["end"])})

for pid, genes in by_plasmid.items():
    for g in genes:
        if g["orf_id"] in defence_orfs:
            islands[pid].append({"name": "defence", "start": g["start"], "end": g["end"]})
        if g["orf_id"] in conj_orfs:
            islands[pid].append({"name": "conj", "start": g["start"], "end": g["end"]})

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
# Per ORF, only where non-empty: its term sources, and whether its window is the plasmid.
term_sources, covers_plasmid = {}, set()
for pid, genes in by_plasmid.items():
    units = directons(genes, max_gap=cfg["max_operon_gap"], circular=pid in circular,
                      length=length_of.get(pid))
    unit_of = {oid: i for i, unit in enumerate(units) for oid in unit}
    flanks_of = flanks(genes, window=cfg["neighbourhood_window"], circular=pid in circular)
    plasmid_islands = islands.get(pid, [])

    for g in genes:
        oid = g["orf_id"]
        left, right = flanks_of[oid]
        near = left + right
        ctx = {island["name"] for island in overlapping_islands(g, plasmid_islands)}

        if any(class_of.get(n) == "FUNCTIONAL" for n in near):
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

        sources = orf_term_sources(oid, near, directon=unit, label_terms=label_terms,
                                   system_terms=system_terms)
        if sources:
            term_sources[oid] = sources
        if window_covers_plasmid(near, len(genes)):
            covers_plasmid.add(oid)

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
                "" if f in not_run else
                round(sum(f in present for present in per_plasmid.values()) / n_units, 6)
                for f in FEATURES])

print(f"context: {n_families} families over {len(by_plasmid)} plasmids")

# ------------------------------------------------------------------------------------
# Context terms per family and lineage: the dark families above, then every other family
# at the primary resolution as the calibration benchmark.
# ------------------------------------------------------------------------------------
with open(snakemake.input.lineage, newline="") as fh:
    lineage_of = {r["plasmid_id"]: r["plasmid_lineage_cluster"]
                  for r in csv.DictReader(fh, delimiter="\t")}


def primary_families():
    """(family_id, members) at the primary resolution, read afresh on each call."""
    with open(snakemake.input.all_families, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["family_resolution"] == snakemake.params.primary:
                yield r["family_id"], r["members"].split(",")


# The focal family of a source ORF, for the tandem-paralogue exclusion. Membership is
# the full primary family, annotated members included, not the dark-only member list.
source_seqs = set(seq_of_source.values())
family_of_seq = {}
for family_id, members in primary_families():
    for m in members:
        if m in source_seqs:
            family_of_seq[m] = family_id
family_of_orf = {oid: family_of_seq[sid] for oid, sid in seq_of_source.items()
                 if sid in family_of_seq}


def occurrences(members):
    return [(oid, term_sources.get(oid, ()), oid in covers_plasmid)
            for m in members for oid in orfs_of_seq.get(m, ())]


n_rows = collections.Counter()
with open(snakemake.output.terms, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t")
    w.writeheader()
    dark_ids = set()
    with open(snakemake.input.families, newline="") as fh:
        for fam in csv.DictReader(fh, delimiter="\t"):
            dark_ids.add(fam["family_id"])
            rows = family_term_rows(fam["family_id"], "dark",
                                    occurrences(fam["members"].split(",")),
                                    family_of_orf, lineage_of)
            w.writerows(rows)
            n_rows["dark"] += len(rows)
    for family_id, members in primary_families():
        if family_id in dark_ids:
            continue
        rows = family_term_rows(family_id, "known", occurrences(members),
                                family_of_orf, lineage_of)
        w.writerows(rows)
        n_rows["known"] += len(rows)

print(f"context terms: {n_rows['dark']} rows for dark families, {n_rows['known']} for "
      f"known families at {snakemake.params.primary}")
