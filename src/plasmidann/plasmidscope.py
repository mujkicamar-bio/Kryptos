"""PlasmidScope's eggNOG annotations, transferred to our proteins by identical sequence.

PlasmidScope (Li et al., Nucleic Acids Res. 2025, 53:D179) ran eggNOG-mapper 2.1.12 on
every protein of every plasmid in its `ALL` set. A protein of ours whose sequence is
identical to one of theirs has, by construction, the same eggNOG result, so re-searching
it through our cascade repeats published work. Proteins it annotates skip the cascade;
everything else - proteins it leaves dark, and proteins it does not contain - is searched
by all five tiers.

WHY THE CLASS IGNORES THE PRODUCT NAME

PlasmidScope's ORFs are not one method. For COMPASS, IMG-PR, PLSDB, mMGE and Kraken2 they
are Prodigal calls with Prokka products; for GenBank, RefSeq, DDBJ and EMBL they are the
deposited CDS with the submitters' products (mostly NCBI PGAP). A named product therefore
means different things for different sources. eggNOG-mapper was run on all of them, so
the eggNOG fields are the only evidence that means the same thing on every row.

WHY THE JOIN IS IDENTITY, NOT COORDINATES

The key is plasmidann.dereplicate's seq_id, a hash of the protein sequence. A protein that
shares a stop codon with a PlasmidScope ORF but starts elsewhere is a different sequence,
and eggNOG's result for one does not transfer to the other.
"""
from plasmidann.cascade import is_informative
from plasmidann.dereplicate import _seq_id

# eggNOG-mapper's placeholder for "no value", which PlasmidScope keeps.
NO_VALUE = "-"

# Most informative first. A sequence that occurs on several plasmids takes the most
# informative class any occurrence received.
CLASSES = ("ANNOTATED", "UNKNOWN_ORTHOLOG", "NONE")
RANK = {c: i for i, c in enumerate(CLASSES)}

# PlasmidScope column -> our column. KEGG_Pathway and GOs are carried for S4b orthology.
FIELDS = {"COG_category": "cog_category", "COG_id": "cog_id", "KEGG_ko": "kegg_ko",
          "KEGG_Pathway": "kegg_pathways", "PFAMs": "pfams", "GOs": "gos",
          "EC_number": "ec"}

# The tier name PlasmidScope-resolved proteins carry in protein_annotation.tsv.
TIER = "PS"


def ps_class(row):
    """Class one PlasmidScope protein row from its eggNOG fields alone.

    ANNOTATED         an informative Pfam name, a KEGG KO or an EC number
    UNKNOWN_ORTHOLOG  eggNOG placed it (a COG/OG, or only DUF/UPF Pfams) but named nothing
    NONE              no eggNOG identifier of any kind

    The bar is the cascade's own: a label must NAME a function (cascade.is_informative).
    A DUF or UPF family is a domain of unknown function, and a COG category letter is a
    broad class, not a function; a protein resting on either is searched by the cascade
    like any other unnamed protein. On the 100-plasmid test set these were 153 and 81 of
    the 2,956 proteins a looser rule had called annotated.

    PlasmidScope writes category S on rows with no identifier as well, so NONE cannot tell
    "no eggNOG hit" from "hit without an identifier". Both are dark here.
    """
    pfams = [] if row["PFAMs"] == NO_VALUE else row["PFAMs"].split(",")
    if any(is_informative(p) for p in pfams) or \
            any(row[k] != NO_VALUE for k in ("KEGG_ko", "EC_number")):
        return "ANNOTATED"
    if pfams or row["COG_id"] != NO_VALUE:
        return "UNKNOWN_ORTHOLOG"
    return "NONE"


def protein_seq_id(row):
    """The seq_id our pipeline gives the same protein. Both sides drop the stop '*'."""
    return _seq_id(row["Sequence"].rstrip("*"))


def reduce_rows(rows, wanted):
    """{seq_id: record} for the PlasmidScope rows whose protein is in `wanted`.

    Each record holds the class, the eggNOG fields (placeholders as '') and the ORF
    source, taken from the most informative occurrence of the sequence.
    """
    out = {}
    for row in rows:
        sid = protein_seq_id(row)
        if sid not in wanted:
            continue
        cls = ps_class(row)
        known = out.get(sid)
        if known is not None:
            known["n_orfs"] += 1
            if RANK[cls] >= RANK[known["ps_class"]]:
                continue
        rec = {"seq_id": sid, "ps_class": cls,
               "orf_source": row["Orf Prediction Source"],
               "n_orfs": known["n_orfs"] if known else 1}
        for src, dst in FIELDS.items():
            rec[dst] = "" if row[src] == NO_VALUE else row[src]
        out[sid] = rec
    return out


def annot_label(rec):
    """The label a PlasmidScope-resolved protein carries: first informative Pfam, else KO,
    else EC. A DUF listed before a named family must not become the protein's label."""
    named = [p for p in rec.get("pfams", "").split(",") if is_informative(p)]
    if named:
        return named[0]
    for key in ("kegg_ko", "ec"):
        if rec.get(key):
            return rec[key].split(",")[0]
    return ""
