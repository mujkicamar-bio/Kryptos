"""Parser for eggNOG-mapper's `.emapper.annotations` table.

eggNOG-mapper is run on the proteins the cascade named, to give them controlled terms
(COG, KEGG, GO, EC, Pfam) where the cascade gives free text. It writes a bare '-'
(NO_VALUE, the placeholder PlasmidScope's eggNOG-mapper tables carry too) in every field it
has nothing for; that placeholder is read as absence, because taken literally it would
become a COG category or pathway named '-' shared by most proteins.
"""
from plasmidann.plasmidscope import NO_VALUE


def _clean(value):
    value = (value or "").strip()
    return "" if value == NO_VALUE else value


def _clean_list(value):
    """Parse a comma-separated emapper field into a list, absence as an empty list.

    emapper separates multi-valued fields with commas and writes a bare '-' when it has
    nothing.
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
        records[query] = {
            # Several single-letter categories are written adjacently, as "EGP". The table
            # keeps the string whole; labels.labels_from_orthology splits it into letters.
            "cog_category": _clean(row.get("COG_category")),
            "kegg_pathways": _clean_list(row.get("KEGG_Pathway")),
            "preferred_name": _clean(row.get("Preferred_name")),
            "description": _clean(row.get("Description")),
            "eggnog_ogs": _clean(row.get("eggNOG_OGs")),
            # The Pfam families eggNOG assigns to the orthologous group, GO terms, EC
            # numbers and KEGG orthologs.
            "pfams": _clean_list(row.get("PFAMs")),
            "gos": _clean_list(row.get("GOs")),
            "ec": _clean_list(row.get("EC")),
            "kegg_ko": _clean_list(row.get("KEGG_ko")),
        }
    return records

