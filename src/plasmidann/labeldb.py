"""Plasmid-specific label databases: when a search hit becomes a label.

WHAT THIS ADDS TO THE LABEL VOCABULARY

The cascade and eggNOG name a protein; they do not say that it is a relaxase of the MOBP
family, a type II antitoxin, a mercury-resistance gene or an anti-CRISPR. Seven curated
databases and AMRFinderPlus say exactly that, each within its own domain, and each is its own label kind so
that a statement from one is never mistaken for a statement from another:

    source     label kind          label                         sub_label
    tadb       tadb_ta             'type II toxin'               -
    bacmet     bacmet_compound     compound ('Mercury (Hg)')     gene; BacMet's class
    oritdb     oritdb_role         relaxase / auxiliary / T4CP   MOB or T4CP family
    card       card_amr_family     CARD AMR gene family          the model name (TEM-1)
    mobileog   mobileog_category   major category                minor category; evidence
    dbapis     dbapis_family       APIS family (APIS030)         gene; evidence
    acrdb      acrdb_family        Acr family (AcrIF1)           CRISPR type; evidence
    amrfinder  amrfinder_gene      element symbol (blaTEM-1)     element type/subtype

PlasAnn's database and labels are not used; only its tier thresholds are, from the paper.

WHY COVERAGE ON BOTH SEQUENCES

The tiers are PlasAnn's (Islam et al. 2026, Nucleic Acids Res., Methods, "Annotation
pipeline"). The paper states identity and coverage; this module requires the coverage on
the query AND on the subject. Query coverage alone lets a 40-residue fragment carry the
label of a 400-residue reference, and subject coverage alone lets a multidomain protein
carry the label of one of its domains. Measured on the 100-plasmid test set against the
PlasAnn database, the stricter rule labels 1,188 proteins where query coverage alone
labels 1,267.

WHY CARD HAS ITS OWN RULE

CARD curates a bitscore cut-off for every protein homolog model, and its Perfect / Strict
paradigm is defined on it (Alcock et al. 2023, Nucleic Acids Res. 51:D690). A generic
identity tier would ignore the curation. Only protein homolog models are used: variant,
rRNA, overexpression and knockout models detect resistance from mutations or absence,
which the presence of a similar protein cannot show.
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
# Every target is reported (--max-target-seqs 0), pre-filtered at the tier-2 minimum: a
# target limit ranks by bitscore, not tier, and on the 100-plasmid test set it dropped
# tier-1 oriTDB and mobileOG-db hits.
DIAMOND_TIERED_ARGS = (f"--more-sensitive --evalue 1e-5 --max-target-seqs 0 "
                       f"--id {TIER2_IDENTITY:g} --query-cover {TIER2_COVERAGE:g} "
                       f"--subject-cover {TIER2_COVERAGE:g}")
# DIAMOND for CARD exactly as RGI runs it (app/Diamond.py): --more-sensitive with DIAMOND's
# default e-value and target count. The curated bitscore cut-off, not the e-value, decides
# a Strict call. These settings reproduced the RGI calls on all 16 AMR-positive plasmids of
# the 100-plasmid test set.
DIAMOND_CARD_ARGS = "--more-sensitive"

# AMRFinderPlus Type/Subtype values (NCBI AMRFinderPlus documentation, "Output format").
_AMRFINDER_TYPES = {("AMR", "AMR"), ("AMR", "POINT"), ("STRESS", "METAL"),
                    ("STRESS", "BIOCIDE"), ("STRESS", "ACID"), ("STRESS", "HEAT"),
                    ("VIRULENCE", "VIRULENCE"), ("VIRULENCE", "ANTIGEN")}


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


def read_bacmet_compounds(path):
    """{BacMet_ID: Compound} from BacMet's own mapping file (BacMet2_EXP.753.mapping.txt)."""
    with open(path, newline="") as fh:
        return {r["BacMet_ID"]: r["Compound"] for r in csv.DictReader(fh, delimiter="\t")}


def parse_bacmet(header, compounds):
    """BacMet 2.0: one label per compound, as BacMet's mapping file names it, with the gene
    and, when BacMet gives one, the compound class as sub_label:
    'Triclosan [class: Phenolic compounds]' -> ('Triclosan', 'abeM; class=Phenolic compounds').
    The final '.' that ends some compound lists is list punctuation and is dropped.
    """
    subject, gene = header.split("|")[:2]
    out = set()
    for part in compounds[subject].split(", "):
        compound, _, cls = part.strip().rstrip(".").partition("[class:")
        if compound.strip():
            out.add((compound.strip(), f"{gene}; class={cls.strip(' ]')}" if cls else gene))
    return subject, sorted(out)


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
    gene=NA marks a homologue whose family has no verified seed, so it names no gene. An
    accession dbAPIS assigns to several families carries them comma-joined in family=, and
    gives one label per family.
    """
    fields = _fields(header)
    if not fields.get("family") or not fields.get("evidence"):
        raise ValueError(f"dbAPIS header lacks family= or evidence=: {header!r}")
    gene = fields.get("gene", "")
    detail = _with_evidence("" if gene == "NA" else gene, fields["evidence"])
    return header.split()[0], [(f, detail) for f in fields["family"].split(",")]


def parse_acrdb(header):
    """Anti-CRISPRdb version 2.2 as installed from its core dataset: '<anti_CRISPR_id>
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
# Output column names of AMRFinderPlus 4.
_AMR_COLUMNS = {
    "seq_id": "Protein id",
    "symbol": "Element symbol",
    "type": "Type",
    "subtype": "Subtype",
    "method": "Method",
    "scov": "% Coverage of reference",
    "pident": "% Identity to reference",
    "subject": "Closest reference accession",
}


def parse_amrfinder(path, version):
    """Label rows from an `amrfinder -p <faa> --plus` table.

    Every element becomes an amrfinder_gene label with Type/Subtype as sub_label; its
    Method (EXACTP, BLASTP, PARTIALP, HMM...) is the tier. AMRFinderPlus reports identity
    and coverage of the reference only, so qcov and bitscore stay empty, as does every 'NA'.
    An element type this module does not know is refused, because its context-term type
    (plasmidann.context_terms.label_term) would otherwise be guessed.
    """
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            value = {k: ("" if r[c] == "NA" else r[c]) for k, c in _AMR_COLUMNS.items()}
            if (value["type"], value["subtype"]) not in _AMRFINDER_TYPES:
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
