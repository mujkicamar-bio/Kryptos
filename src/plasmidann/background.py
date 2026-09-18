"""Stage 13: explicit backgrounds for every prevalence and enrichment claim (section 52).

THE RULE THIS MODULE ENFORCES

Section 52.2, stated as an absolute: "No normalized prevalence field is valid without its
normalization definition." So every normalised number this module produces travels with
three fields - normalization_method, normalization_version, background_definition - and the
raw count and raw prevalence are retained beside it (section 52.1). A normalised value on
its own is not a result; it is a number whose meaning depends entirely on a choice the
reader cannot see.

WHY A FLAT CORPUS BACKGROUND IS NOT ENOUGH

Section 53 requires context enrichment to be compared against "an appropriate background of
non-dark ORFs or another explicitly justified reference population", and lists size-matched,
stratified and permutation approaches.

The reason is the trap this pipeline already documents at S8. On a 5 kb cryptic plasmid
carrying six genes, a plus or minus three neighbourhood IS the entire molecule, so
everything co-occurs with everything. A flat corpus background divides by one number for
the whole collection, which under-corrects for small plasmids and over-corrects for large
ones. Small cryptic plasmids are a stratum of interest here, so the error lands exactly
where it does most damage.

A stratified background compares a dark ORF only against ORFs on plasmids of similar size,
gene count and habitat - so the comparison answers "is this association unusual FOR A
PLASMID LIKE THIS", which is the question.

THE STRATA

Section 52 lists the covariates that matter: plasmid_length, gene_count, host_taxonomy,
MOB_group, plasmid_lineage_cluster, habitat, source_database, clonal redundancy. Which of
those a given run stratifies on is configuration, not a constant here, because adding
strata makes each one smaller and at some point the background is too thin to estimate.
`min_stratum_size` is the floor, and a stratum below it falls back to the pooled background
with the fallback RECORDED rather than applied silently.

LENGTH AND GENE COUNT ARE BANDED, NOT USED RAW

A continuous covariate cannot define a stratum. Bands are logarithmic because plasmid size
spans four orders of magnitude and a linear band would put a 2 kb and a 9 kb plasmid in
different strata while pooling everything above 100 kb.
"""
import collections
import math

# Bumped whenever the stratification or the arithmetic changes in a way that makes two runs
# incomparable. Recorded on every normalised row (section 52.2), so a table can always be
# matched to the code that produced it.
NORMALIZATION_VERSION = "1"

METHOD_STRATIFIED = "stratified_prevalence_ratio"
METHOD_POOLED = "pooled_prevalence_ratio"


def length_band(length, base=2.0):
    """A logarithmic size band label for a plasmid length.

    Logarithmic because plasmid size spans four orders of magnitude: a linear band wide
    enough to be useful at 200 kb would put every cryptic plasmid in one bucket, and one
    narrow enough for cryptic plasmids would give hundreds of near-empty strata at the top.

    The label is the band's lower bound in kb, so it reads as a size rather than an index.
    """
    try:
        length = float(length)
    except (TypeError, ValueError):
        return "unknown"
    if length <= 0:
        return "unknown"
    exponent = math.floor(math.log(length, base))
    return f"{base ** exponent / 1000:.3g}kb"


def gene_count_band(n_genes):
    """Gene count bands. Same reasoning as length, on the same scale."""
    try:
        n_genes = int(n_genes)
    except (TypeError, ValueError):
        return "unknown"
    if n_genes <= 0:
        return "unknown"
    return str(2 ** math.floor(math.log2(n_genes)))


def stratum_key(plasmid, covariates):
    """The stratum a plasmid belongs to, as a tuple over the configured covariates.

    A covariate the metadata does not carry contributes 'unknown' rather than being
    dropped. Dropping it would silently merge plasmids with missing metadata into the
    strata of plasmids that have it, which makes metadata completeness look like biology -
    the same error the family definition once made.
    """
    values = []
    for name in covariates:
        if name == "plasmid_length":
            values.append(length_band(plasmid.get("size_bp")))
        elif name == "gene_count":
            values.append(gene_count_band(plasmid.get("n_genes")))
        else:
            values.append(str(plasmid.get(name) or "unknown"))
    return tuple(values)


def build_background(plasmids, covariates):
    """Group plasmids into strata. Returns {stratum_key: [plasmid_id, ...]}.

    `plasmids` maps plasmid_id to its metadata dict.
    """
    strata = collections.defaultdict(list)
    for plasmid_id, meta in plasmids.items():
        strata[stratum_key(meta, covariates)].append(plasmid_id)
    return strata


def normalise(observed_count, observed_total, background_count, background_total,
              method, background_definition):
    """A prevalence ratio with the provenance section 52.2 requires.

    Returns raw_count and raw_prevalence (section 52.1, always retained) beside the
    normalised value and the three fields that define it. The normalised value is never
    returned alone.

    A zero background prevalence yields an empty normalised value rather than infinity: an
    infinite ratio is not a measurement, and a reader ranking on it would put every
    never-before-seen association above every well-supported one.
    """
    raw_prevalence = observed_count / observed_total if observed_total else 0.0
    background_prevalence = (background_count / background_total
                             if background_total else 0.0)
    normalised = (round(raw_prevalence / background_prevalence, 6)
                  if background_prevalence else "")
    return {
        # section 52.1 - always retained, whatever normalisation was applied
        "raw_count": observed_count,
        "raw_prevalence": round(raw_prevalence, 6),
        "background_count": background_count,
        "background_total": background_total,
        "background_prevalence": round(background_prevalence, 6),
        # section 52.2 - a normalised value is invalid without these three
        "normalized_prevalence": normalised,
        "normalization_method": method,
        "normalization_version": NORMALIZATION_VERSION,
        "background_definition": background_definition,
    }


def describe_background(covariates, stratum, n_plasmids, fell_back):
    """A human-readable background_definition string for one comparison.

    It names the covariates, the stratum's values and its size, so a reader can tell
    whether the comparison rests on forty plasmids or on four - and whether the stratified
    background was used at all.
    """
    if fell_back:
        return (f"pooled over all {n_plasmids} plasmids (stratum "
                f"{'/'.join(stratum)} below min_stratum_size)")
    return (f"stratified on {'+'.join(covariates)} = {'/'.join(stratum)}; "
            f"{n_plasmids} plasmids in stratum")
