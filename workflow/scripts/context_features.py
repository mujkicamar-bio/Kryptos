"""S8c: genomic context per ORF, aggregated to families and tested against a background.

WHAT CHANGED AND WHY

This stage used to describe a dark ORF's neighbourhood with six hand-picked features, two
of which - backbone_adjacent and ta_candidate - were decided by a hand-written list of 73
Pfam family names. Pfam-A 38.2 holds 30,134 families, of which 67 mention replication in
their description and 42 mention conjugation; the list named 16 and 15, named no MobB and
no MobD, and nine of its names did not exist in Pfam-A at all, so those entries had never
once matched anything and no output could have revealed it.

A neighbourhood is now described by the labels the tools actually produced for the
neighbours, read from results/08_protein_labels/protein_labels.tsv. Those labels are grouped into
categories through plasmidann.categories, which with no rules configured makes each
(kind, label) its own category. The biological grouping - replication, mobilisation,
conjugation - is derived later from the observed vocabulary and applied here as a
configuration change.

THE TABLE IS LONG, NOT WIDE

Six features fitted in six column pairs. An open vocabulary does not: there are 30,134 Pfam
families alone, and nr product names are unbounded. One row per (family, category) is the
only shape that holds it, and it is the shape the later grouping reads naturally.

CATEGORY AND SUBCATEGORY

The test runs on the flat CATEGORY. The SUBCATEGORY is carried beside it and never tested,
because the two want opposite granularity: every category tested is another hypothesis and
the false-discovery-rate correction weakens as their number grows, while a reader of the
table wants the finer distinction. Carrying both means a finer grouping can be re-cut from
the output without re-running anything.

THE STATISTICAL TRAP THIS STAGE EXISTS TO AVOID

On a 5 kb cryptic plasmid carrying six genes, a plus or minus three neighbourhood IS the
entire plasmid. Everything co-occurs with everything, and a raw co-occurrence frequency
would rank the smallest plasmids as the most informative when they are the least - which
would be catastrophic here, because small cryptic plasmids are a stratum of interest. So
every association is reported as enrichment over a background, with a Fisher exact
p-value and a Benjamini-Hochberg q-value.

THE BACKGROUND IS STRATIFIED (spec sections 52-53)

It used to be the whole corpus: one rate for the entire collection. That does not avoid
the trap above, it only moves it. A flat background under-corrects for small plasmids,
where the window is the molecule, and over-corrects for large ones - so the correction is
weakest exactly where the artefact is strongest.

A family is now compared against the plasmids in ITS OWN STRATA, banded on the covariates
configured in targets.yaml. The question becomes "is this association unusual for a
plasmid like this", which is the question being asked.

A stratum below background.min_stratum_size is too thin to estimate a rate from, and the
comparison falls back to the pooled collection. Section 52.2 makes that fallback something
that has to be visible: "No normalized prevalence field is valid without its normalization
definition." So every row carries normalization_method, normalization_version and a
background_definition naming the stratum and its size, with the raw counts beside them.

THE UNIT IS THE PLASMID

Not the family member. Members of a family are homologs on plasmids that are frequently
near-identical, so counting members makes sequencing effort look like evidence. Clonal
redundancy between distinct plasmids remains uncorrected; see docs/PARAMETER_PROVENANCE.md.
"""
import _ctx  # noqa: F401
import collections
import csv

from plasmidann import background, categories, enrich
from plasmidann.context import directons, neighbourhood, overlapping_islands

cfg = snakemake.params.context
bg_cfg = snakemake.params.background
COVARIATES = list(bg_cfg["covariates"])
MIN_STRATUM = int(bg_cfg["min_stratum_size"])
cmap = categories.CategoryMap(categories.load_rules(snakemake.params.get("categories")))

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
# Which ORFs correspond to which unique protein. BOTH directions, once.
#
# The reverse direction matters: the per-family loop below previously rescanned the whole
# 9.3M-entry forward map per family, which is O(families x ORFs) - measured at 0.37-0.39 s
# per family, projecting to 10 hours to 3.6 days single-core for ~50k families.
# ------------------------------------------------------------------------------------
seq_of_orf = {}
orfs_of_seq = collections.defaultdict(list)
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, ids = line.rstrip("\n").split("\t")
        for oid in ids.split(","):
            seq_of_orf[oid] = sid
            orfs_of_seq[sid].append(oid)

# ------------------------------------------------------------------------------------
# The categories of every unique protein, from the labels the tools produced.
#
# Read per protein, not per ORF: a protein identical on forty plasmids was labelled once.
# Each entry is a (category, subcategory) pair.
# ------------------------------------------------------------------------------------
categories_of_seq = collections.defaultdict(set)
n_label_rows = 0
with open(snakemake.input.labels, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        n_label_rows += 1
        categories_of_seq[r["protein_id"]] |= set(
            cmap.categories_for(r["kind"], r["label"]))

# ------------------------------------------------------------------------------------
# Islands: defence systems and integron cassette arrays, as intervals per plasmid.
#
# These are system-level calls rather than per-protein labels, so they stay intervals: a
# dark ORF INSIDE a defence system is a different statement from one merely beside a
# defence component.
# ------------------------------------------------------------------------------------
islands = collections.defaultdict(list)
defence_orfs = {}
with open(snakemake.input.defence, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get("orf_id"):
            # The system's own name from MacSyFinder's model, not a class of our invention.
            defence_orfs[r["orf_id"]] = (r.get("system") or "").rsplit("/", 1)[-1]

for f in snakemake.input.integrons:
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                islands[r["plasmid_id"]].append(
                    {"name": "integron", "start": int(r["start"]), "end": int(r["end"]),
                     "detail": r.get("integron_type", "")})
            except (ValueError, KeyError):
                continue

for pid, genes in by_plasmid.items():
    for g in genes:
        system = defence_orfs.get(g["orf_id"])
        if system:
            # The category stays the generic 'defence'; the SYSTEM's own name from
            # MacSyFinder's model travels as the subcategory. Making the system part of
            # the category would split one hypothesis into hundreds, each too rare to
            # test, and would break the defence_island stratum that reads 'defence'.
            islands[pid].append({"name": "defence", "start": g["start"],
                                 "end": g["end"], "detail": system})

# ------------------------------------------------------------------------------------
# Per-ORF context: the categories of its neighbours, its directon partners, and the
# islands it sits inside. Every entry is a (category, subcategory) pair.
# ------------------------------------------------------------------------------------
context_of = {}
for pid, genes in by_plasmid.items():
    units = directons(genes, max_gap=cfg["max_operon_gap"])
    unit_of = {oid: i for i, unit in enumerate(units) for oid in unit}
    plasmid_islands = islands.get(pid, [])

    for g in genes:
        oid = g["orf_id"]
        ctx = set()

        # EVERY island, not the first one found. Defence intervals are appended after
        # integron intervals, so returning the first made a dark ORF inside a defence
        # system that also sat in a cassette array label as 'integron' only.
        # The island's own detail - the defence system's name, the integron's class -
        # becomes the subcategory, so it is in the table without multiplying the number of
        # categories the test has to correct across.
        for island in overlapping_islands(g, plasmid_islands):
            ctx.add((island["name"], island.get("detail", "")))

        # What the neighbours are, in the tools' own words. No hand-assigned class: a
        # neighbour labelled MobA_MobL contributes the category of MobA_MobL, which with no
        # rules configured is 'pfam_family:MobA_MobL'.
        for n in neighbourhood(genes, oid, window=cfg["neighbourhood_window"]):
            ctx |= categories_of_seq.get(seq_of_orf.get(n, ""), set())
            if class_of.get(n) == "FUNCTIONAL":
                ctx.add(("annotated_neighbour", ""))

        # Directon membership with at least one annotated partner: the dark ORF is
        # predicted to be co-transcribed with something we understand, which is a stronger
        # claim than adjacency.
        unit = units[unit_of[oid]] if oid in unit_of else [oid]
        partners = [p for p in unit if p != oid]
        if partners and any(class_of.get(p) == "FUNCTIONAL" for p in partners):
            ctx.add(("operon_with_annotated", ""))
        # The geometry of a tight two-gene unit, recorded as geometry. It used to be read
        # as a toxin-antitoxin call whenever the partner appeared on the 73-name list;
        # whether the partner's category makes it one is a question for the grouping, not
        # for this stage, and the partner's categories are already in ctx for that to be
        # answered without re-running anything.
        if len(unit) == 2:
            ctx.add(("two_gene_operon", ""))

        context_of[oid] = ctx

# ------------------------------------------------------------------------------------
# Which stratum each plasmid belongs to (Stage 13, spec section 52).
#
# size_bp comes from the master table because it is the record's own length; gene_count is
# counted from THIS run's annotation, because that is the number the neighbourhood window
# was applied to. Taking gene count from metadata instead would band a plasmid by a figure
# the run never used.
# ------------------------------------------------------------------------------------
size_of = {}
with open(snakemake.input.master, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        size_of[r["plasmid_id"]] = r.get("size_bp")

stratum_of = {}
for pid, genes in by_plasmid.items():
    stratum_of[pid] = background.stratum_key(
        {"size_bp": size_of.get(pid), "n_genes": len(genes)}, COVARIATES)

# ------------------------------------------------------------------------------------
# Background counts, per PLASMID, both pooled and per stratum. A plasmid counts for a
# category if any of its ORFs has that category in context. Counting ORFs instead would
# let one gene-dense plasmid dominate.
# ------------------------------------------------------------------------------------
plasmids_with = collections.Counter()
stratum_with = collections.defaultdict(collections.Counter)
stratum_size = collections.Counter()

for pid, genes in by_plasmid.items():
    present = set()
    for g in genes:
        present |= context_of.get(g["orf_id"], set())
    stratum = stratum_of[pid]
    stratum_size[stratum] += 1
    # Counted on the CATEGORY alone. Counting (category, subcategory) pairs would make the
    # background depend on how finely the categories happen to be subdivided, so adding a
    # subcategory would silently change every enrichment value in the table.
    for category in {c for c, _ in present}:
        plasmids_with[category] += 1
        stratum_with[stratum][category] += 1

n_plasmids = len(by_plasmid)


def background_for(strata, category):
    """(count, total, method, definition) for one category against one family's strata.

    A family may span several strata - it is a family precisely because it recurs - so the
    reference population is the union of the strata its plasmids occupy. That keeps the
    family's own plasmids inside its background, which is what makes K >= k hold and what
    enrich.fisher_enrichment requires.

    Below MIN_STRATUM plasmids the stratum cannot support a rate, and the comparison falls
    back to the pooled collection. The fallback is returned, never applied quietly.
    """
    total = sum(stratum_size[st] for st in strata)
    if total < MIN_STRATUM:
        return (plasmids_with.get(category, 0), n_plasmids, background.METHOD_POOLED,
                background.describe_background(
                    COVARIATES, sorted(strata)[0] if strata else (), n_plasmids, True))
    count = sum(stratum_with[st].get(category, 0) for st in strata)
    # One stratum is named in full; several are summarised, because a family spanning
    # forty strata would otherwise write forty band labels into every row.
    if len(strata) == 1:
        described = background.describe_background(
            COVARIATES, next(iter(strata)), total, False)
    else:
        described = (f"stratified on {'+'.join(COVARIATES)}; union of {len(strata)} "
                     f"strata the family occupies; {total} plasmids")
    return count, total, background.METHOD_STRATIFIED, described


with open(snakemake.output.background, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["stratum", "covariates", "category", "n_units_with_category", "n_units",
                "background_rate", "usable"])
    # The pooled row stays, labelled as pooled. It is the fallback population, so a reader
    # checking a row that fell back needs to find it in this file.
    for category, count in sorted(plasmids_with.items()):
        w.writerow(["POOLED", "+".join(COVARIATES), category, count, n_plasmids,
                    round(count / n_plasmids, 6) if n_plasmids else 0.0, 1])
    for stratum in sorted(stratum_size):
        size = stratum_size[stratum]
        for category, count in sorted(stratum_with[stratum].items()):
            w.writerow(["/".join(stratum), "+".join(COVARIATES), category, count, size,
                        round(count / size, 6), int(size >= MIN_STRATUM)])

# ------------------------------------------------------------------------------------
# Per family: one row per category present in any member's context, tested against the
# background, with q-values corrected across the categories tested for that family.
# ------------------------------------------------------------------------------------
cols = ["family_id", "category", "subcategories", "n_units", "n_units_with_category",
        "observed_rate", "background_rate", "enrichment", "odds_ratio", "p_value",
        "q_value", "status",
        # Stage 13 provenance. Section 52.1 keeps the raw counts; section 52.2 makes the
        # three definition fields mandatory beside any normalised value.
        "raw_count", "raw_prevalence", "background_count", "background_total",
        "background_prevalence", "normalized_prevalence", "normalization_method",
        "normalization_version", "background_definition"]

n_rows = n_families = 0
with open(snakemake.output.families, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    with open(snakemake.input.families, newline="") as fh:
        for fam in csv.DictReader(fh, delimiter="\t"):
            n_families += 1
            seq_members = fam["members"].split(",")
            orfs = [o for m in seq_members for o in orfs_of_seq.get(m, ())]

            # The unit is the plasmid. orf_id is '<plasmid_id>|<ordinal>'.
            per_plasmid = collections.defaultdict(set)
            for oid in orfs:
                per_plasmid[oid.rsplit("|", 1)[0]] |= context_of.get(oid, set())
            n_units = len(per_plasmid)
            if not n_units:
                continue

            # Counts on the flat category; the subcategories seen for it are collected
            # beside the count and reported, never tested.
            counts = collections.Counter()
            subcategories = collections.defaultdict(set)
            for present in per_plasmid.values():
                for category, subcategory in present:
                    if subcategory:
                        subcategories[category].add(subcategory)
                for category in {c for c, _ in present}:
                    counts[category] += 1

            # The strata this family occupies. Its background is the union of them; see
            # background_for.
            strata = {stratum_of[pid] for pid in per_plasmid if pid in stratum_of}

            tested = []
            for category, k in sorted(counts.items()):
                K, N, method, definition = background_for(strata, category)
                tested.append((category, k, enrich.fisher_enrichment(
                    k=k, n=n_units, K=K, N=N),
                    background.normalise(k, n_units, K, N, method, definition)))

            # Corrected across the categories tested for THIS family. A family tested
            # against four hundred categories and one tested against three are not
            # comparable without it, and the open vocabulary makes that the normal case.
            q_values = enrich.benjamini_hochberg([r["p_value"] for _, _, r, _ in tested])

            for (category, k, result, norm), q in zip(tested, q_values):
                n_rows += 1
                row = {
                    "family_id": fam["family_id"], "category": category,
                    "subcategories": ",".join(sorted(subcategories.get(category, ()))),
                    "n_units": n_units, "n_units_with_category": k,
                    "observed_rate": result["observed_rate"],
                    "background_rate": result["background_rate"],
                    "enrichment": result["enrichment"],
                    "odds_ratio": result["odds_ratio"],
                    "p_value": result["p_value"], "q_value": q,
                    "status": result["status"]}
                row.update(norm)
                w.writerow(row)

grouping = "per label (no category rules configured)" if cmap.is_identity else "per category"
usable = sum(1 for st in stratum_size if stratum_size[st] >= MIN_STRATUM)
print(f"context: {n_label_rows} label rows, {len(plasmids_with)} categories over "
      f"{n_plasmids} plasmids, {n_rows} association rows for {n_families} families, "
      f"grouping is {grouping}")
# How much of the collection the stratified background can actually serve. If few strata
# clear the floor, most rows fell back to pooled and the stratification is nominal - which
# is a configuration problem (too many covariates) rather than a result.
print(f"         background stratified on {'+'.join(COVARIATES)}: {len(stratum_size)} "
      f"strata, {usable} at or above min_stratum_size={MIN_STRATUM}")
