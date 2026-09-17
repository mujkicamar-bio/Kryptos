"""S4b: orthology terms for the proteins the cascade DID name.

WHY A DARK-PROTEIN PIPELINE ANNOTATES THE KNOWN FRACTION HARDER

eggNOG-mapper is not a dark-hunting tier and must not be used as one: over Swiss-Prot and
nr it adds very little for finding homologues. It is here for the opposite reason.

S8 asks what a dark ORF's neighbours DO. The cascade answers in free text - "Aspartokinase",
"MULTISPECIES: hypothetical protein [Enterobacteriaceae]" - and free text cannot be
aggregated. Two neighbours in the same pathway, described differently by two databases, look
like two unrelated observations. FESNov's neighbourhood metric is defined over KEGG pathway
membership precisely because a pathway is a term you can count.

So annotating the known fraction properly is what makes the unknown fraction interpretable.
It is the annotated genes around a dark ORF that turn it into a hypothesis.

THE PLACEHOLDER TRAP

eggNOG-mapper writes a bare "-" in every field it has nothing for. Read literally that
becomes a COG category named "-" and a KEGG pathway named "-", and both then aggregate into
the neighbourhood composition as though they were real terms - producing a "-" pathway
shared by most of the corpus, which would be the single most enriched context feature in
the run.
"""

# eggNOG-mapper's placeholder for "no value". Anything equal to this is absence, not a term.
NO_VALUE = "-"


def _clean(value):
    value = (value or "").strip()
    return "" if value == NO_VALUE else value


def _clean_list(value):
    """Parse a comma-separated emapper field into a list, absence as an empty list.

    emapper separates multi-valued fields with commas and writes a bare '-' when it has
    nothing. An empty list is the honest representation of absence; a list containing '-'
    would aggregate as a real term, which is the placeholder trap in the module docstring
    applied to four more columns.
    """
    cleaned = _clean(value)
    return [v for v in (p.strip() for p in cleaned.split(",")) if v] if cleaned else []


def parse_annotations(text):
    """Parse an `.emapper.annotations` file into {query: record}.

    emapper writes '##' banner lines before the header and a '## Total time' line after the
    last row, so a parser that accepts any non-header line produces a record keyed
    '## Total time'. Both are skipped here.

    The header line itself begins '#query'; every other '#' line is commentary.
    """
    records = {}
    header = None
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("#query"):
            header = line.lstrip("#").split("\t")
            header[0] = "query"
            continue
        if line.startswith("#"):
            continue
        if header is None:
            continue
        fields = line.split("\t")
        row = dict(zip(header, fields))
        query = _clean(row.get("query"))
        if not query:
            continue
        pathways = _clean(row.get("KEGG_Pathway"))
        records[query] = {
            # Several single-letter categories are written adjacently, as "EGP". They are
            # kept whole: splitting them into characters invents three annotations, and
            # keeping only the first discards two.
            "cog_category": _clean(row.get("COG_category")),
            "kegg_pathways": [p for p in pathways.split(",") if p] if pathways else [],
            "preferred_name": _clean(row.get("Preferred_name")),
            "description": _clean(row.get("Description")),
            "eggnog_ogs": _clean(row.get("eggNOG_OGs")),
            # Four columns emapper always writes and this parser used to discard. They are
            # the controlled-vocabulary identifiers on the annotated fraction: the Pfam
            # families eggNOG assigns to the orthologous group, GO terms, EC numbers and
            # KEGG orthologs. A functional grouping has to be derived from vocabularies
            # like these rather than written by hand, so discarding them made the grouping
            # impossible to build.
            "pfams": _clean_list(row.get("PFAMs")),
            "gos": _clean_list(row.get("GOs")),
            "ec": _clean_list(row.get("EC")),
            "kegg_ko": _clean_list(row.get("KEGG_ko")),
        }
    return records


def cog_category(record):
    """The COG category of a parsed record, or "" when there is none."""
    return (record or {}).get("cog_category", "")


def kegg_pathways(record):
    """The KEGG pathways of a parsed record, or [] when there are none."""
    return list((record or {}).get("kegg_pathways", []))
