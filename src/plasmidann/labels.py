"""The functional-label vocabulary: what each tool said, and what kind of statement it is.

Every label is copied verbatim from the tool that produced it and tagged with its KIND. No
biological role is assigned here.

The kind travels with every label because 'repA' as an eggNOG gene symbol and 'RepA_N' as
a Pfam family are different statements from different evidence; downstream code groups
on (kind, label), never on the label alone.

Not labels: organism names (the source of the reference sequence, not a statement about
the query; the full title stays in hits.tsv), uninformative descriptions such as
'hypothetical protein' (recorded by the cascade as dark evidence), and empty strings.
"""

from plasmidann import labeldb

# Every kind of label protein_labels admits; it refuses any other.
KINDS = frozenset({
    # hmmsearch against Pfam-A, tiers T1 and T2
    "pfam_family",        # RepA_N            the family name, as tier_search records it
    "pfam_description",   # 'Replication initiator protein A (RepA) N-terminus'
    "pfam_clan",          # CL0123            families too divergent to align as one
    # DIAMOND against NCBI's rendering of Swiss-Prot, and against ClusteredNR
    "swissprot_product",  # 'Toxin CcdB'      from 'RecName: Full=...'
    "nr_product",         # 'conjugal transfer protein TraG'   the ClusteredNR title
    # pharokka in protein mode: the phage families, plus CARD and VFDB from the same run
    "pharokka_annotation",  # 'terminase large subunit'   the family's own name
    "pharokka_category",    # 'head and packaging'        the curators' functional group
    "card_gene_family",     # 'TEM beta-lactamase'        CARD's AMR gene family
    "card_mechanism",       # 'antibiotic inactivation'   CARD's resistance mechanism
    "vfdb_factor",          # 'type IV pilus'             VFDB's virulence factor
    # eggNOG-mapper
    "gene_symbol",        # repA, traG, mobA  the most systematic axis available
    "cog_category",       # L, D, V           one label per letter
    "cog_id",             # COG5527
    "eggnog_pfam",        # the Pfam families eggNOG assigns to the orthologous group
    "eggnog_description",
    "go",                 # GO:0006270
    "ec",                 # 2.7.7.7
    "kegg_ko",            # ko:K02314
}) | labeldb.KINDS  # the plasmid label databases: card_amr_family, amrfinder_gene, ...

# The database a hits.tsv row came from, declared per tier in config/cascade.yaml; it
# decides the label kind.
SOURCES = frozenset({"pfam", "pharokka", "card", "vfdb", "swissprot", "nr"})

def _informative(row):
    """Whether the cascade judged a hits.tsv row informative (csv writes 'True'/'False')."""
    return row.get("informative") == "True"


def parse_ncbi_title(title):
    """Split an NCBI protein title into accession, product name and organism.

    Two shapes are handled, both verified against the databases in use:

        swissprot   P62554.1 RecName: Full=Toxin CcdB; AltName: Full=Protein LetD [E. coli]
        nr          WP_000813620.1 type II toxin-antitoxin system RelE/ParE toxin [E. coli]

    A title that matches neither returns the whole string as the product. That is the
    honest fallback: nr titles are not uniform, this parser will meet shapes it was not
    shown, and a dropped label is invisible while a strange one is not. NCBI's record
    prefixes (MULTISPECIES:, MAG:, TPA:) stay in the product as NCBI writes them.
    """
    text = (title or "").strip()
    if not text:
        return {"accession": "", "product": "", "organism": ""}

    organism = ""
    if text.endswith("]") and "[" in text:
        head, _, tail = text.rpartition("[")
        organism = tail[:-1].strip()
        text = head.strip()

    accession, product = "", text
    first, _, rest = text.partition(" ")
    # An accession is the leading token when it looks like one: NCBI and UniProt
    # accessions carry a version suffix and no spaces. Testing for the '.' rather than
    # matching a pattern per database keeps this from needing a rule per accession style.
    if rest and "." in first:
        accession, product = first, rest.strip()

    if "RecName:" in product:
        # 'RecName: Full=Toxin CcdB; AltName: Full=Protein LetD' - the recommended name is
        # the first Full=, and the AltName synonyms are dropped: they are the same protein
        # under other names, so admitting them would multiply one statement into several.
        # A qualifier before RecName ('PUTATIVE PSEUDOGENE: RecName: ...', seen in the
        # installed swissprot) is a statement about the gene and is kept in front.
        prefix, _, body = product.partition("RecName:")
        head = body.strip().split(";", 1)[0].strip()
        product = prefix + (head[len("Full="):].strip() if head.startswith("Full=") else head)

    return {"accession": accession, "product": product, "organism": organism}


def labels_from_hit(row, pfam=None):
    """Labels from one hits.tsv row.

    `pfam` is the mapping from plasmidann.pfam_meta.load; when it is None or lacks the
    family, only the family name is emitted. A family the installed release does not know
    is a database mismatch, and the name has to survive so the mismatch appears in the
    table rather than silently removing the hit.
    """
    if not _informative(row):
        return []

    source = (row.get("source") or "").strip()
    if source not in SOURCES:
        raise ValueError(
            f"hits.tsv row for {row.get('query')!r} at tier {row.get('tier')!r} has "
            f"source {source!r}, which is not one of {sorted(SOURCES)}. The label kind "
            "is decided by the source; a row without one cannot be labelled correctly.")

    label = (row.get("label") or "").strip()
    accession = (row.get("target_accession") or "").strip()
    if not label:
        return []

    out = []
    if source == "pfam":
        out.append({"kind": "pfam_family", "label": label, "accession": accession})
        meta = (pfam or {}).get(label)
        if meta:
            if meta.get("description"):
                out.append({"kind": "pfam_description", "label": meta["description"],
                            "accession": accession})
            if meta.get("clan"):
                out.append({"kind": "pfam_clan", "label": meta["clan"],
                            "accession": accession})
        return out

    category = (row.get("category") or "").strip()
    if source == "pharokka":
        out.append({"kind": "pharokka_annotation", "label": label, "accession": accession})
        if category:
            out.append({"kind": "pharokka_category", "label": category,
                        "accession": accession})
        return out
    if source == "card":
        out.append({"kind": "card_gene_family", "label": label, "accession": accession})
        if category:
            out.append({"kind": "card_mechanism", "label": category,
                        "accession": accession})
        return out
    if source == "vfdb":
        return [{"kind": "vfdb_factor", "label": label, "accession": accession}]

    parsed = parse_ncbi_title(label)
    if not parsed["product"]:
        return []
    kind = "swissprot_product" if source == "swissprot" else "nr_product"
    return [{"kind": kind, "label": parsed["product"],
             "accession": accession or parsed["accession"]}]


def _split(value):
    """Comma-separated table field to a list, with blanks removed."""
    return [v for v in (p.strip() for p in (value or "").split(",")) if v]


def labels_from_orthology(row):
    """Labels from one orthology.tsv row.

    Every controlled-vocabulary identifier eggNOG assigns, plus the gene symbol. The gene
    symbol is singled out because it is the most systematic axis available anywhere in this
    pipeline: symbols are assigned in families (rep*, tra*, trb*, mob*, par*, tnp*, ccd*,
    rel*, vap*, hig*), which free-text descriptions are not.
    """
    out = []
    symbol = (row.get("preferred_name") or "").strip()
    if symbol:
        out.append({"kind": "gene_symbol", "label": symbol, "accession": ""})

    # Adjacent letters, as 'LKV', are three separate COG categories. The orthology table
    # keeps the string whole for provenance; the vocabulary needs each category, because
    # 'LKV' is not a category any COG release defines.
    for letter in (row.get("cog_category") or "").strip():
        if letter.isalpha():
            out.append({"kind": "cog_category", "label": letter, "accession": ""})

    # eggNOG_OGs is 'COG5527@2,2QV1F@1224' - the group identifier at each taxonomic level.
    # The COG identifiers are the citable ones; the numeric eggNOG groups have no published
    # functional description, so admitting them would add labels nothing can interpret.
    for group in _split(row.get("eggnog_ogs")):
        name = group.split("@", 1)[0]
        if name.startswith("COG"):
            out.append({"kind": "cog_id", "label": name, "accession": ""})

    for family in _split(row.get("pfams")):
        out.append({"kind": "eggnog_pfam", "label": family, "accession": ""})
    for term in _split(row.get("gos")):
        out.append({"kind": "go", "label": term, "accession": ""})
    for number in _split(row.get("ec")):
        out.append({"kind": "ec", "label": number, "accession": ""})
    for ko in _split(row.get("kegg_ko")):
        out.append({"kind": "kegg_ko", "label": ko, "accession": ""})

    description = (row.get("eggnog_description") or "").strip()
    if description:
        out.append({"kind": "eggnog_description", "label": description, "accession": ""})
    return out
