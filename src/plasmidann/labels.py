"""The functional-label vocabulary: what each tool said, and what kind of statement it is.

WHY THIS REPLACES A CURATED LIST

An earlier version of this pipeline classified plasmid proteins into replication,
mobilisation, conjugation, partition, transposition, toxin-antitoxin and
restriction-modification using a hand-written list of 73 Pfam family names. It could not
work, and the measurements say why:

  * Pfam-A 38.2 holds 30,134 families, of which 67 mention replication in their
    description and 42 mention conjugation. The list named 16 and 15. It named no MobB and
    no MobD.
  * pfam2go, the published Pfam-to-GO mapping, covers 4 of the 16 listed replication
    families, 1 of 16 conjugation families and 0 of 11 mobilisation families - so the gaps
    are in the curated mapping too, not only in the list.
  * Nine of the 73 names did not exist in Pfam-A at all, so those entries had never once
    matched anything, and no output could have revealed it.
  * The role assignments had no source. They were one person's reading.

A hand list cannot reach the scope of a 30,134-family database, and its gaps are silent.
So no role is assigned here. Every label is taken verbatim from the tool that produced it,
tagged with the KIND of statement it is, and the grouping into biological categories is
derived later from the vocabulary actually observed in the data.

WHY THE KIND IS LOAD-BEARING

'repA' as a gene symbol from eggNOG and 'RepA_N' as a Pfam family from hmmsearch are
different statements about a protein, from different evidence, with different reliability.
A vocabulary that could not tell them apart could not be audited, and a category built from
both would be impossible to describe in a methods section. The kind travels with every
label for that reason, and downstream code groups on (kind, label), never on the label
alone.

WHAT IS DELIBERATELY NOT A LABEL

  * Organism names. A functional grouping is not a taxonomy, and an organism in a title is
    the source of the reference sequence, not a statement about the query. The full title
    stays in hits.tsv either way, so nothing is lost from the record.
  * Uninformative descriptions. 'hypothetical protein' names nothing; the cascade already
    records it as dark evidence. Admitted here it would become the most frequent, and
    therefore most apparently enriched, category in the collection.
  * The empty string in any field. Absence is absence, and a label of '' would group every
    protein missing that field into one enormous false category.
"""

# Every kind of label this module can emit. A kind not listed here is a statement nothing
# downstream knows how to group, so emitting one is a bug rather than a new feature.
#
# Deliberately absent: any biological role name. Roles are derived from these kinds later;
# a role appearing in this set would be the curated list growing back, and
# tests/test_labels.py asserts that it has not.
KINDS = frozenset({
    # hmmsearch against Pfam-A, tiers T1 and T2
    "pfam_family",        # RepA_N            the family name, as tier_search records it
    "pfam_description",   # 'Replication initiator protein A (RepA) N-terminus'
    "pfam_clan",          # CL0123            families too divergent to align as one
    # DIAMOND, tier T3 (NCBI rendering of Swiss-Prot) and tier T4 (nr)
    "swissprot_product",  # 'Toxin CcdB'      from 'RecName: Full=...'
    "pgap_product",       # 'conjugal transfer protein TraG'
    # eggNOG-mapper, S4b
    "gene_symbol",        # repA, traG, mobA  the most systematic axis available
    "cog_category",       # L, D, V           one label per letter
    "cog_id",             # COG5527
    "eggnog_pfam",        # the Pfam families eggNOG assigns to the orthologous group
    "eggnog_description",
    "go",                 # GO:0006270
    "ec",                 # 2.7.7.7
    "kegg_ko",            # ko:K02314
    # MacSyFinder, S8a phase 2
    "macsy_system",       # RM_Type_II        the model's own system name
    "macsy_component",    # RM_Type_II_REase
    # IntegronFinder, S8b
    "integron_element",   # intI, attC, attI
    "integron_type",      # complete, In0, CALIN
})

# Tiers whose label is a Pfam family name rather than a sequence title. Taken from the tier
# method in config/cascade.yaml: the hmmer tiers search Pfam-A and report family names.
_PFAM_TIERS = frozenset({"T1", "T2"})

# The tier whose database is NCBI's rendering of Swiss-Prot. Its titles carry
# 'RecName: Full=<name>;' rather than the UniProt 'OS=/GN=/PE=' structure, and notably no
# gene symbol - verified against the installed database, not assumed.
_SWISSPROT_TIERS = frozenset({"T3"})

# NCBI marks a title shared by several organisms with this prefix. 'MULTISPECIES: relaxase'
# and 'relaxase' are the same product, so keeping the prefix would split every widespread
# protein into two labels - exactly the proteins a grouping most needs to see as one.
_MULTISPECIES = "MULTISPECIES:"


def _informative(row):
    """Whether a hits.tsv row was judged informative by the cascade.

    The column is written by csv as the string 'True' or 'False'. Comparing the raw value
    to a boolean would make every row falsy and silently empty the vocabulary, so the
    accepted spellings are tested explicitly.
    """
    return row.get("informative") in (True, "True", "true", "1", 1)


def parse_ncbi_title(title):
    """Split an NCBI protein title into accession, product name and organism.

    Two shapes are handled, both verified against the databases in use:

        swissprot   P62554.1 RecName: Full=Toxin CcdB; AltName: Full=Protein LetD [E. coli]
        nr          WP_000813620.1 type II toxin-antitoxin system RelE/ParE toxin [E. coli]

    A title that matches neither returns the whole string as the product. That is the
    honest fallback: nr titles are not uniform, this parser will meet shapes it was not
    shown, and a dropped label is invisible while a strange one is not.
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

    if product.startswith(_MULTISPECIES):
        product = product[len(_MULTISPECIES):].strip()

    if product.startswith("RecName:"):
        # 'RecName: Full=Toxin CcdB; AltName: Full=Protein LetD' - the recommended name is
        # the first Full=, and the AltName synonyms are dropped: they are the same protein
        # under other names, so admitting them would multiply one statement into several.
        body = product[len("RecName:"):].strip()
        head = body.split(";", 1)[0].strip()
        product = head[len("Full="):].strip() if head.startswith("Full=") else head

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

    tier = row.get("tier", "")
    label = (row.get("label") or "").strip()
    accession = (row.get("target_accession") or "").strip()
    if not label:
        return []

    out = []
    if tier in _PFAM_TIERS:
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

    parsed = parse_ncbi_title(label)
    if not parsed["product"]:
        return []
    kind = "swissprot_product" if tier in _SWISSPROT_TIERS else "pgap_product"
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


def labels_from_defence(row):
    """Labels from one defence_systems.tsv row.

    `system` is MacSyFinder's model_fqn, a path such as
    'defense-finder-models/Defense/RM_Type_II'. The system name is its last element; the
    full path is kept as the accession so the model set that made the call stays visible in
    the table, which is what makes the call traceable to its publication.
    """
    out = []
    fqn = (row.get("system") or "").strip()
    if fqn:
        out.append({"kind": "macsy_system", "label": fqn.rsplit("/", 1)[-1],
                    "accession": fqn})
    component = (row.get("component") or "").strip()
    if component:
        out.append({"kind": "macsy_component", "label": component, "accession": fqn})
    return out


def labels_from_integron(row):
    """Labels from one integron row: the element type and the integron class."""
    out = []
    element = (row.get("annotation") or "").strip()
    if element:
        out.append({"kind": "integron_element", "label": element, "accession": ""})
    integron_type = (row.get("integron_type") or "").strip()
    if integron_type:
        out.append({"kind": "integron_type", "label": integron_type, "accession": ""})
    return out
