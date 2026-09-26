"""Plasmid-specific label databases: when a search hit becomes a label, and where two disagree.

WHAT THIS ADDS TO THE LABEL VOCABULARY

The cascade and eggNOG name a protein; they do not say that it is a relaxase of the MOBP
family, a type II antitoxin, a mercury-resistance gene or an anti-CRISPR. Seven curated
databases say exactly that, each within its own domain, and each is its own label kind so
that a statement from one is never mistaken for a statement from another:

    source     label kind          label                         sub_label
    tadb       tadb_ta             'type II toxin'               -
    bacmet     bacmet_compound     'Mercury', 'Acridine'         the gene (merA)
    oritdb     oritdb_role         relaxase / auxiliary / T4CP   MOB or T4CP family
    card       card_amr_family     CARD AMR gene family          the model name (TEM-1)
    mobileog   mobileog_category   major category                minor category; evidence
    dbapis     dbapis_family       APIS family (APIS030)         gene; evidence
    acrdb      acrdb_family        Acr family (AcrIF1)           CRISPR type; evidence
    amrfinder  amrfinder_gene      element symbol (blaTEM-1)     element type/subtype

PlasAnn's own database is not used, and neither are its labels (user decision, 2026-09-25):
only its tier thresholds are taken, from the paper.

WHY COVERAGE ON BOTH SEQUENCES

The tiers are PlasAnn's (Islam et al. 2026, Nucleic Acids Res., Methods, "Annotation
pipeline"). The paper states identity and coverage; this module requires the coverage on
the query AND on the subject. Query coverage alone lets a 40-residue fragment carry the
label of a 400-residue reference, and subject coverage alone lets a multidomain protein
carry the label of one of its domains. Measured on the 100-plasmid test set against the
PlasAnn database (Step 0), the stricter rule labels 1,188 proteins where query coverage
alone labels 1,267.

WHY CARD HAS ITS OWN RULE

CARD curates a bitscore cut-off for every protein homolog model, and its Perfect / Strict
paradigm is defined on it (Alcock et al. 2023, Nucleic Acids Res. 51:D690). A generic
identity tier would ignore the curation. Only protein homolog models are used: variant,
rRNA, overexpression and knockout models detect resistance from mutations or absence,
which the presence of a similar protein cannot show.

THE DISAGREEMENT TABLE CHANGES NOTHING

`disagreements` lists every protein on which two sources make incompatible statements. It
reads labels and returns new rows; no label is removed, re-ranked or rewritten because
another source disagrees. Deciding which source is right is a review question, and a
pipeline that settled it silently would hide the evidence the review needs.
"""
import csv
import json
import pathlib
import re

# Label kind per source. labels.KINDS includes these so protein_labels admits them.
KIND = {
    "card": "card_amr_family",
    "amrfinder": "amrfinder_gene",
    "tadb": "tadb_ta",
    "bacmet": "bacmet_compound",
    "oritdb": "oritdb_role",
    "mobileog": "mobileog_category",
    "dbapis": "dbapis_family",
    "acrdb": "acrdb_family",
}
KINDS = frozenset(KIND.values())

# The databases searched with the identity/coverage tiers, and every database under
# data/refs/labels/<db>/. AMRFinderPlus is a separate program with its own database.
TIERED = ("tadb", "bacmet", "oritdb", "mobileog", "dbapis", "acrdb")
DATABASES = TIERED + ("card",)

COLUMNS = ["seq_id", "source", "label_kind", "label", "sub_label", "tier", "cut_off",
           "pident", "qcov", "scov", "bitscore", "subject", "database_version"]
DISAGREEMENT_COLUMNS = ["seq_id", "source_a", "label_a", "source_b", "label_b",
                        "conflict_type"]

# Islam et al. 2026, Nucleic Acids Res. (PlasAnn), Methods "Annotation pipeline": primary
# tier >=80% identity and >=90% coverage, secondary tier >60% identity and >70% coverage.
# Coverage is required on query and subject alike (module docstring).
TIER1_IDENTITY = 80.0
TIER1_COVERAGE = 90.0
TIER2_IDENTITY = 60.0
TIER2_COVERAGE = 70.0

# DIAMOND for the tiered databases. e-value 1e-5 is PlasAnn's search threshold (Islam et
# al. 2026, released code essential_annotation.py); --more-sensitive is the mode RGI passes
# to DIAMOND (arpcard/rgi app/Diamond.py), used for every database so the sensitivity is
# the same across sources.
# Every target is reported (--max-target-seqs 0), pre-filtered at the tier-2 minimum, which
# no tier-qualifying hit falls below. A target limit is ranked by bitscore, not by tier, so
# a partial high-scoring hit can push the tier-1 hit out: measured on the 100-plasmid test
# set (leaf 1.2, G2 and G4), 50 targets without the pre-filter (Step 0's search) lost 1
# oriTDB and 3 mobileOG-db labels and gave 2 + 2 more the wrong tier, and even 1,000 targets
# truncated oriTDB queries that have more qualifying relaxases than that. The pre-filter
# keeps the output small: 17,856 oriTDB and 83,770 mobileOG-db (all 775,257 entries) rows
# for 5,304 proteins.
DIAMOND_TIERED_ARGS = (f"--more-sensitive --evalue 1e-5 --max-target-seqs 0 "
                       f"--id {TIER2_IDENTITY:g} --query-cover {TIER2_COVERAGE:g} "
                       f"--subject-cover {TIER2_COVERAGE:g}")
# DIAMOND for CARD exactly as RGI runs it (app/Diamond.py): --more-sensitive with DIAMOND's
# default e-value and target count. The curated bitscore cut-off, not the e-value, decides
# a Strict call; Step 0 reproduced the July RGI calls on all 16 positive test plasmids.
DIAMOND_CARD_ARGS = "--more-sensitive"

# Term prefixes for context terms (contract of the labels-build plan). AMRFinderPlus is
# split by element type below; other STRESS subtypes and VIRULENCE get no prefix.
_PREFIX = {"card": "amr", "bacmet": "metal", "tadb": "ta", "oritdb": "conj_role",
           "mobileog": "mge", "dbapis": "antidefence", "acrdb": "antidefence"}
# AMRFinderPlus Type/Subtype values (NCBI AMRFinderPlus documentation, "Output format").
_AMRFINDER_PREFIX = {("AMR", "AMR"): "amr", ("AMR", "POINT"): "amr",
                     ("STRESS", "METAL"): "metal", ("STRESS", "BIOCIDE"): "metal",
                     ("STRESS", "ACID"): "", ("STRESS", "HEAT"): "",
                     ("VIRULENCE", "VIRULENCE"): "", ("VIRULENCE", "ANTIGEN"): ""}


def tier(pident, qcov, scov):
    """1, 2 or None for one hit: identity in percent, coverages in percent of each sequence."""
    coverage = min(qcov, scov)
    if pident >= TIER1_IDENTITY and coverage >= TIER1_COVERAGE:
        return 1
    if pident > TIER2_IDENTITY and coverage > TIER2_COVERAGE:
        return 2
    return None


def best_tiered_hits(hits):
    """{query: (tier, hit)}: the best tier, then the highest bitscore, then the subject key.

    The better tier wins over a higher bitscore, because the tier is the statement a label
    is cited with; the key breaks exact ties so the result does not depend on hit order.
    """
    best = {}
    for hit in hits:
        t = tier(hit["pident"], hit["qcov"], hit["scov"])
        if t is not None:
            _keep_best(best, hit, t, (-t, hit["bitscore"]))
    return {q: (t, hit) for q, (_, t, hit) in best.items()}


def _keep_best(best, hit, grade, rank):
    """Keep `hit` for its query if its rank is higher, or equal with a smaller subject key."""
    current = best.get(hit["query"])
    if (current is None or rank > current[0]
            or (rank == current[0] and int(hit["key"]) < int(current[2]["key"]))):
        best[hit["query"]] = (rank, grade, hit)


def card_call(pident, sstart, send, slen, bitscore, cut_off):
    """'Perfect', 'Strict' or None (Loose, discarded) for one CARD protein homolog hit.

    Alcock et al. 2023 and CARD's model description: "A Perfect RGI match is 100% identical
    to the reference protein sequence along its entire length, a Strict RGI match is not
    identical but the bit-score of the matched sequence is greater than the curated BLASTP
    bit-score cutoff". RGI's code (app/HomologModel.py) tests `>=`, which is followed here.
    """
    if pident == 100.0 and sstart == 1 and send == slen:
        return "Perfect"
    if bitscore >= cut_off:
        return "Strict"
    return None


def best_card_hits(hits, models):
    """{query: (call, hit)}: a Perfect hit if any, else the Strict hit with the top bitscore.

    Hits to a model not in `models` (not a protein homolog model) are ignored.
    """
    best = {}
    for hit in hits:
        model = models.get(hit["key"])
        if model is None:
            continue
        call = card_call(hit["pident"], hit["sstart"], hit["send"], hit["slen"],
                         hit["bitscore"], model["cut_off"])
        if call is not None:
            _keep_best(best, hit, call, (call == "Perfect", hit["bitscore"]))
    return {q: (call, hit) for q, (_, call, hit) in best.items()}


def card_models(path):
    """{model_id: entry} for every CARD protein homolog model in card.json.

    The cut-off is the model's curated `model_param.blastp_bit_score`. There is one label
    per AMR gene family the model belongs to, with the model name as sub_label.
    """
    data = json.loads(pathlib.Path(path).read_text())
    models = {}
    for model_id, m in data.items():
        if model_id.startswith("_") or m["model_type"] != "protein homolog model":
            continue
        (seq,) = m["model_sequences"]["sequence"].values()
        families = sorted({c["category_aro_name"] for c in m["ARO_category"].values()
                           if c["category_aro_class_name"] == "AMR Gene Family"})
        models[model_id] = {
            "key": model_id, "subject": "ARO:" + m["ARO_accession"],
            "sequence": seq["protein_sequence"]["sequence"].upper(),
            "cut_off": float(m["model_param"]["blastp_bit_score"]["param_value"]),
            "labels": [(f, m["model_name"]) for f in families],
        }
    return models


# ------------------------------------------------------------------------------------
# Reference FASTA files and their headers
# ------------------------------------------------------------------------------------
def fasta_records(path):
    """Yield (header, sequence): the full header line and the upper-case sequence without
    whitespace or a terminal stop '*'. plasmidann.fasta.iter_fasta keeps only the
    identifier, and here the labels are in the rest of the header.
    """
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper().rstrip("*")
                header, chunks = line[1:].strip(), []
            else:
                # Whitespace is never sequence. Anti-CRISPRdb's core dataset carries
                # spaces inside some sequences, and DIAMOND rejects a FASTA holding one.
                chunks.append("".join(line.split()))
    if header is not None:
        yield header, "".join(chunks).upper().rstrip("*")


_TADB_STEM = re.compile(r"type_([IVX]+)_(AT|T)")


def parse_tadb(header):
    """TADB 3.0 as installed: '<TADB id>|<file stem> <original header>'.

    TADB distributes one file per type and role (type_II_T_exp.fas, regulator_exp.fas) and
    its own headers carry neither, so the installer keeps the file stem in the id.
    """
    subject, stem = header.split()[0].split("|")
    if stem == "regulator":
        return subject, [("regulator", "")]
    m = _TADB_STEM.fullmatch(stem)
    if not m:
        raise ValueError(f"TADB header {header!r} does not name a type and role "
                         "(type_<N>_T, type_<N>_AT or regulator)")
    role = "antitoxin" if m.group(2) == "AT" else "toxin"
    return subject, [(f"type {m.group(1)} {role}", "")]


_METAL = re.compile(r"([A-Z][a-z]+) \([A-Za-z]+\)")


def read_bacmet_compounds(path):
    """{BacMet_ID: Compound} from BacMet's own mapping file (BacMet2_EXP.753.mapping.txt)."""
    with open(path, newline="") as fh:
        return {r["BacMet_ID"]: r["Compound"] for r in csv.DictReader(fh, delimiter="\t")}


def parse_bacmet(header, compounds):
    """BacMet 2.0: one label per compound, with the gene as sub_label.

    A metal is labelled by its element ('Mercury (Hg)' -> 'Mercury'), a biocide by the
    class BacMet assigns ('Triclosan [class: Phenolic compounds]' -> 'Phenolic compounds'):
    the class is the level at which neighbours can share a term.
    """
    subject, gene = header.split("|")[:2]
    out = set()
    for part in compounds[subject].split(", "):
        part = part.strip().rstrip(".")
        if "[class:" in part:
            out.add(part.split("[class:", 1)[1].strip(" ]"))
        elif _METAL.fullmatch(part):
            out.add(_METAL.fullmatch(part).group(1))
        elif part:
            out.add(part)
    return subject, [(label, gene) for label in sorted(out)]


_ORITDB_ROLE = {"relaxase": "relaxase", "auxiliary": "auxiliary protein", "t4cp": "T4CP"}


def parse_oritdb(header):
    """oriTDB 2.0 as installed: '<role>_<n> <original header>', role relaxase / auxiliary /
    t4cp from the file it came from. The subject is oriTDB's entry name (TraI_RP4), and the
    family (MOBP, VirD4/TraG) is the original header's token before 'id=', '' when absent.
    """
    tokens = header.split()
    role = _ORITDB_ROLE.get(tokens[0].split("_")[0])
    if role is None:
        raise ValueError(f"oriTDB header {header!r} does not start with a role "
                         f"({', '.join(sorted(_ORITDB_ROLE))})")
    original = tokens[1:]
    before_id = next((i for i, t in enumerate(original) if t.startswith("id=")),
                     len(original))
    family = original[before_id - 1] if before_id > 2 else ""
    return original[0], [(role, "" if family == "_" else family)]


def _with_evidence(detail, evidence):
    """sub_label '<detail>; evidence=<class>', or 'evidence=<class>' when there is no detail.

    Every evidence class of a database is kept, and the class of the reference entry that
    gave the label travels with it, so a label from a curated entry and one from a predicted
    or keyword-recovered entry can be told apart downstream.
    """
    return "; ".join(p for p in (detail, f"evidence={evidence}") if p)


def _fields(header):
    """{key: value} of the 'key=value' tokens of an installed header."""
    return dict(t.split("=", 1) for t in header.split()[1:] if "=" in t)


def parse_mobileog(header):
    """mobileOG-db beatrix-1.6: 'id|name|uniprot|major|minor|element|evidence', evidence
    Manual, Homology or Keyword Search. The element field can hold a space ('Plasmid
    RefSeq'), so the header is split on '|' alone.
    """
    fields = header.split("|")
    subject, major, minor, evidence = fields[0], fields[3], fields[4], fields[6]
    return subject, [(major, _with_evidence("" if minor == "NA" else minor, evidence))]


def parse_dbapis(header):
    """dbAPIS as installed: '<accession> gene=<gene> family=<family> evidence=verified|homolog
    <description>'. A verified seed that formed no family has its gene as family (gp54).
    """
    fields = _fields(header)
    if not fields.get("family") or not fields.get("evidence"):
        raise ValueError(f"dbAPIS header lacks family= or evidence=: {header!r}")
    return header.split()[0], [(fields["family"],
                                _with_evidence(fields.get("gene", ""), fields["evidence"]))]


def parse_acrdb(header):
    """Anti-CRISPRdb v2.2 as installed from its core dataset: '<anti_CRISPR_id>
    family=<Family> type=<Anti_type> acc=<Accession> evidence=<Verified|PLiterature|
    Putative>'. The family is the label; the CRISPR type it inhibits and the evidence class
    are the sub_label.
    """
    fields = _fields(header)
    if not fields.get("family") or not fields.get("evidence"):
        raise ValueError(f"Anti-CRISPRdb header lacks family= or evidence=: {header!r}")
    return header.split()[0], [(fields["family"],
                                _with_evidence(fields.get("type", ""), fields["evidence"]))]


def load_reference(db, directory):
    """Every entry of one tiered database, from <directory>/<db>.faa as the installer
    (tools/download_label_dbs.py) writes it: [{key, subject, sequence, labels, cut_off}].

    Keys are sequential numbers in file order, used as the DIAMOND subject id: oriTDB entry
    names repeat across roles, and a number maps back unambiguously.
    """
    directory = pathlib.Path(directory)
    if db == "bacmet":
        # BacMet's own mapping file, which the installer keeps under raw/.
        (mapping,) = directory.rglob("*mapping*.txt")
        compounds = read_bacmet_compounds(mapping)
    parse = {"tadb": parse_tadb, "oritdb": parse_oritdb, "mobileog": parse_mobileog,
             "dbapis": parse_dbapis, "acrdb": parse_acrdb,
             "bacmet": lambda header: parse_bacmet(header, compounds)}[db]
    entries = []
    for header, seq in fasta_records(directory / f"{db}.faa"):
        subject, labels = parse(header)
        entries.append({"key": str(len(entries)), "subject": subject, "sequence": seq,
                        "labels": labels, "cut_off": ""})
    return entries


def label_rows(source, best, entries, version):
    """Rows in COLUMNS order for {query: (tier or call, hit)} against {key: entry}."""
    rows = []
    for query in sorted(best):
        grade, hit = best[query]
        entry = entries[hit["key"]]
        for label, sub_label in entry["labels"]:
            rows.append({
                "seq_id": query, "source": source, "label_kind": KIND[source],
                "label": label, "sub_label": sub_label, "tier": str(grade),
                "cut_off": str(entry["cut_off"]), "pident": str(hit["pident"]),
                "qcov": str(hit["qcov"]), "scov": str(hit["scov"]),
                "bitscore": str(hit["bitscore"]), "subject": entry["subject"],
                "database_version": version})
    return rows


# ------------------------------------------------------------------------------------
# AMRFinderPlus (Feldgarden et al. 2021, Sci. Rep. 11:12728), protein mode with --plus
# ------------------------------------------------------------------------------------
# Output column names in AMRFinderPlus 4, then the name the same column has in 3.x.
_AMR_COLUMNS = {
    "seq_id": ("Protein id", "Protein identifier"),
    "symbol": ("Element symbol", "Gene symbol"),
    "type": ("Type", "Element type"),
    "subtype": ("Subtype", "Element subtype"),
    "method": ("Method",),
    "scov": ("% Coverage of reference", "% Coverage of reference sequence"),
    "pident": ("% Identity to reference", "% Identity to reference sequence"),
    "subject": ("Closest reference accession", "Accession of closest sequence"),
}


def parse_amrfinder(path, version):
    """Label rows from an `amrfinder -p <faa> --plus` table.

    Every element becomes an amrfinder_gene label with Type/Subtype as sub_label; its
    Method (EXACTP, BLASTP, PARTIALP, HMM...) is the tier. AMRFinderPlus reports identity
    and coverage of the reference only, so qcov and bitscore stay empty, as does every 'NA'.
    An element type this module does not know is refused, because its term prefix would
    otherwise be guessed.
    """
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        col = {k: next(n for n in names if n in reader.fieldnames)
               for k, names in _AMR_COLUMNS.items()}
        for r in reader:
            value = {k: ("" if r[c] == "NA" else r[c]) for k, c in col.items()}
            if (value["type"], value["subtype"]) not in _AMRFINDER_PREFIX:
                raise ValueError(f"AMRFinderPlus element type {value['type']}/"
                                 f"{value['subtype']} is not known to plasmidann.labeldb")
            rows.append({
                "seq_id": value["seq_id"], "source": "amrfinder",
                "label_kind": KIND["amrfinder"], "label": value["symbol"],
                "sub_label": f"{value['type']}/{value['subtype']}",
                "tier": value["method"], "cut_off": "", "pident": value["pident"],
                "qcov": "", "scov": value["scov"], "bitscore": "",
                "subject": value["subject"], "database_version": version})
    return rows


def term_prefix(row):
    """The context-term prefix of a label row, '' when it has none."""
    if row["source"] == "amrfinder":
        return _AMRFINDER_PREFIX[tuple(row["sub_label"].split("/", 1))]
    return _PREFIX[row["source"]]


# ------------------------------------------------------------------------------------
# Disagreements
# ------------------------------------------------------------------------------------
def read_protein_map(path):
    """{orf_id: seq_id} from 03_dereplication/protein_map.tsv (seq_id, comma-joined ORFs)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            seq_id, orfs = line.rstrip("\n").split("\t")
            for orf in orfs.split(","):
                out[orf] = seq_id
    return out


def by_protein(rows, orf_to_protein):
    """{seq_id: {(system, component)}} from a per-ORF system table (DefenseFinder, CONJScan).

    The system is the last element of MacSyFinder's model path. ORFs absent from the map
    (not in the dereplicated protein set) are skipped.
    """
    out = {}
    for r in rows:
        seq_id = orf_to_protein.get(r["orf_id"])
        if seq_id is not None:
            out.setdefault(seq_id, set()).add((r["system"].rsplit("/", 1)[-1],
                                               r["component"]))
    return out


def read_ko_symbols(path):
    """{KO: [gene symbols]} from KEGG's KO list (rest.kegg.jp/list/ko), whose lines read
    'K00002<TAB>AKR1A1, adh; alcohol dehydrogenase (NADP+) [EC:1.1.1.2]'. A KO whose
    description has no ';' names no symbol, and neither does one whose symbol is the KO
    itself ('K00243; uncharacterized protein', 1,185 KOs in the 2026-09-25 list).
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            ko, desc = line.rstrip("\n").split("\t")
            symbols = desc.split(";")[0] if ";" in desc else ""
            out[ko] = [s.strip() for s in symbols.split(",")
                       if s.strip() and not re.fullmatch(r"K\d{5}", s.strip())]
    return out


def tier0_symbols(orthology_rows, ko_symbols):
    """{seq_id: gene symbols} for the Tier 0 proteins of orthology.tsv.

    PlasmidScope's eggNOG annotation carries KOs but, measured on the test set, not one
    eggNOG preferred name (0 of 2,724 Tier 0 proteins), so the symbols are the KEGG gene
    symbols of its KOs, plus the preferred name where one exists.
    """
    out = {}
    for r in orthology_rows:
        if r["orthology_source"] != "plasmidscope":
            continue
        symbols = set()
        if r["preferred_name"] not in ("", "-"):
            symbols.add(r["preferred_name"])
        for ko in r["kegg_ko"].split(","):
            symbols.update(ko_symbols.get(ko.strip().removeprefix("ko:"), ()))
        if symbols:
            out[r["seq_id"]] = symbols
    return out


def _names(name):
    """Comparable forms of a gene name: lower case, alphanumerics only, with and without
    a leading 'bla' (CARD names beta-lactamases without it, AMRFinderPlus and eggNOG with).
    """
    n = re.sub(r"[^a-z0-9]", "", name.lower())
    forms = {n}
    if n.startswith("bla") and len(n) > 3:
        forms.add(n[3:])
    return {f for f in forms if f}


def same_gene(a, b):
    """True when two gene names are the same gene up to punctuation and allele suffix.

    Equal after normalisation, or one a prefix of the other of at least three characters:
    'tetA' = 'tet(A)', 'bla' ~ 'blaTEM-1', 'aac(6')-Ib' ~ 'aac(6')-Ib-cr', 'TEM-1' =
    'blaTEM-1'. 'sul1' and 'sul2' differ. Three characters is the gene-symbol stem of the
    bacterial nomenclature (Demerec et al. 1966, Genetics 54:61): a shorter prefix names no
    gene, and a stem shared with an allele or locus suffix names the same gene family.
    Synonyms of different stems (aad/ant, virD4/traD) are not resolved and appear as
    conflicts for review.
    """
    for x in _names(a):
        for y in _names(b):
            short, long_ = sorted((x, y), key=len)
            if short == long_ or (len(short) >= 3 and long_.startswith(short)):
                return True
    return False


_GENE_SYMBOL = re.compile(r"[A-Za-z][a-z]{2}[A-Z][A-Za-z0-9]*")


def card_gene(model_name):
    """The gene a CARD model name names, or '' when it names none.

    5,914 of the 6,059 protein homolog model names in CARD 4.0.2 are one token (TEM-1,
    sul1). For the other 145: a gene in brackets at the end is the gene ('PC1
    beta-lactamase (blaZ)'); a leading binomial is the organism ('Streptomyces lividans
    cmlR'); then the first token that is not an English word ('beta-lactamase',
    'intrinsic', 'Trimethoprim-resistant'), a number-led token ('23S') or a '-type'
    qualifier ('mecA-type mecI') is the gene. Checked against all 145, this is right for
    144 - the gene, or '' for the 9 that name none ('Enterococcus faecium chloramphenicol
    acetyltransferase') - and '23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(A)'
    gives 'rRNA'.
    """
    tokens = model_name.split()
    if len(tokens) == 1:
        return model_name
    bracketed = re.search(r"\s\(([^()\s]+)\)$", model_name)
    if bracketed:
        return bracketed.group(1)
    if (len(tokens) > 2 and re.fullmatch(r"[A-Z][a-z]+", tokens[0])
            and re.fullmatch(r"[a-z]+", tokens[1])):
        tokens = tokens[2:]
    for t in tokens:
        # A word is lower case after its first letter and longer than four characters:
        # the only all-lower-case gene tokens among the 145 names are 'cmr' and 'rox'.
        word = len(t) > 4 and re.fullmatch(r"[A-Za-z][a-z]+([-/][a-z]+)*", t)
        if not (word or t[0].isdigit() or t.endswith("-type")):
            return t
    return ""


def _gene(row):
    """The gene a label row names, or '' when its source names none."""
    source = row["source"]
    if source == "bacmet":
        return row["sub_label"]
    if source == "card":
        return card_gene(row["sub_label"])
    if source == "amrfinder":
        return row["label"]
    if source == "oritdb":
        # oriTDB entry names start with the gene (traD_pHCM1, TraI_RP4) - or, for many
        # predicted entries, with a family or a word (t4cp2_..., Relaxase_..., Mob_...).
        # Only a name of the bacterial gene-symbol form, three letters and an upper-case
        # locus letter (Demerec et al. 1966, Genetics 54:61), is taken as a gene.
        name = row["subject"].split("_")[0]
        return name if _GENE_SYMBOL.fullmatch(name) else ""
    return ""


def _conj_role(component):
    """(role, MOB family) of a CONJScan component: 'T4SS_MOBP1' -> ('relaxase', 'MOBP')."""
    name = component.split("T4SS_", 1)[-1]
    if name.upper().startswith("MOB"):
        return "relaxase", re.match(r"MOB[A-Z]", name.upper()).group(0)
    if name.lower() in ("t4cp1", "t4cp2", "tcpa"):
        return "T4CP", ""
    return "MPF", ""


def disagreements(labels, tier0=None, defence=None, conj=None):
    """Rows in DISAGREEMENT_COLUMNS for every cross-source conflict. Changes no label.

    labels   label rows (COLUMNS) from every source of this module
    tier0    {seq_id: gene symbols} of the Tier 0 (PlasmidScope eggNOG) annotation, from
             tier0_symbols
    defence  {seq_id: {(system, component)}} from DefenseFinder, via by_protein
    conj     {seq_id: {(system, component)}} from CONJScan, via by_protein

    Conflict types:
      tier0_vs_<source>     the Tier 0 gene symbol and a source that names a gene (bacmet,
                            card, amrfinder, oritdb) name different genes
      card_vs_amrfinder     both call an AMR gene, and no gene of one is a gene of the other
      bacmet_vs_amrfinder   BacMet and an AMRFinderPlus STRESS METAL/BIOCIDE element name
                            no gene in common
      card_vs_bacmet        both label the protein and name no gene in common
      tadb_vs_defencefinder TADB calls it a toxin-antitoxin gene and DefenseFinder a
                            defence-system component. Recorded for every such protein:
                            no cited source says which defence systems are TA-derived
                            (MazEF, AbiE) and which contradict a TA call, so the pair is
                            listed for review rather than judged by a hand-written list.
      oritdb_vs_conjscan    the roles differ (relaxase / T4CP against CONJScan's MOB* /
                            t4cp / MPF components; an auxiliary protein that CONJScan calls
                            a relaxase or T4CP), or both say relaxase and the MOB families
                            differ. oriTDB 'Other' carries no family to compare.
    A one-sided call (one source labels the protein, the other does not) is not a
    conflict here: it is read directly from protein_labels_plasmid.tsv.
    """
    by_seq = {}
    for row in labels:
        by_seq.setdefault(row["seq_id"], []).append(row)

    out = set()

    def add(seq_id, source_a, label_a, source_b, label_b, conflict):
        out.add((seq_id, source_a, label_a, source_b, label_b, conflict))

    for seq_id, rows in by_seq.items():
        genes = {}
        for row in rows:
            gene = _gene(row)
            if not gene:
                continue
            if row["source"] == "amrfinder":
                kind = row["sub_label"]
                if kind.startswith("AMR/"):
                    genes.setdefault("amrfinder_amr", set()).add(gene)
                elif kind in ("STRESS/METAL", "STRESS/BIOCIDE"):
                    genes.setdefault("amrfinder_metal", set()).add(gene)
                genes.setdefault("amrfinder", set()).add(gene)
            else:
                genes.setdefault(row["source"], set()).add(gene)

        symbols = sorted((tier0 or {}).get(seq_id, ()))
        if symbols:
            for source in ("bacmet", "card", "amrfinder", "oritdb"):
                for gene in sorted(genes.get(source, ())):
                    if not any(same_gene(s, gene) for s in symbols):
                        add(seq_id, "tier0", ",".join(symbols), source, gene,
                            f"tier0_vs_{source}")

        for a, b, key_b, conflict in (("card", "amrfinder", "amrfinder_amr",
                                       "card_vs_amrfinder"),
                                      ("bacmet", "amrfinder", "amrfinder_metal",
                                       "bacmet_vs_amrfinder"),
                                      ("card", "bacmet", "bacmet", "card_vs_bacmet")):
            ga, gb = sorted(genes.get(a, ())), sorted(genes.get(key_b, ()))
            if ga and gb and not any(same_gene(x, y) for x in ga for y in gb):
                add(seq_id, a, ",".join(ga), b, ",".join(gb), conflict)

        components = sorted(f"{s}/{c}" for s, c in (defence or {}).get(seq_id, ()))
        if components:
            for row in rows:
                if row["source"] == "tadb":
                    add(seq_id, "tadb", row["label"], "defencefinder",
                        ",".join(components), "tadb_vs_defencefinder")

        conj_components = sorted(c for _, c in (conj or {}).get(seq_id, ()))
        if conj_components:
            called = [_conj_role(c) for c in conj_components]
            for row in rows:
                if row["source"] != "oritdb":
                    continue
                role, family = row["label"], row["sub_label"]
                if role == "auxiliary protein":
                    agrees = not any(r in ("relaxase", "T4CP") for r, _ in called)
                elif role == "relaxase" and re.fullmatch(r"MOB[A-Z]", family):
                    agrees = ("relaxase", family) in called
                else:
                    agrees = any(r == role for r, _ in called)
                if not agrees:
                    add(seq_id, "oritdb", f"{role} {family}".strip(), "conjscan",
                        ",".join(conj_components), "oritdb_vs_conjscan")

    return [dict(zip(DISAGREEMENT_COLUMNS, r)) for r in sorted(out)]
