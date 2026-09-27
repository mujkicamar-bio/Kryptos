"""The phage tier: parsing pharokka's protein-mode output.

pharokka distributes the prokaryotic virus protein families (PHROGs) as an MMseqs2
profile database and as HMMER3 profiles, and names and categorises every family. Its
`*_full_merged_output.tsv` has one row per input protein, whether or not anything hit it,
with the family, its annotation and category and the statistics of both searches, plus
hits against CARD and VFDB, which it searches in the same run.

The output carries no alignment coordinates: pharokka deletes its raw alignment tables on
exit. A hit therefore has no start and no end, and the cascade reports its completeness as
NOT_MEASURED. Read from pharokka 1.10.1 output on 200 test-set proteins:

  * A protein with no hit carries 'No_MMseqs' in the `phrog` column, not 'No_PHROG', so a
    family hit is recognised by an integer in that column, not by a sentinel string.
  * 90 of 92 hits were found by both searches, each with its own E-value; the stronger
    one is reported, because that is the evidence the assignment rests on.
"""
import csv
import io

# The sentinels pharokka writes into a field with no value. Not relied on to detect a hit
# (see the module docstring); listed so that CARD and VFDB fields can be read as absent.
_ABSENT = frozenset({"None", "No_PHROG", "No_MMseqs", "No_PHROGs_HMM", "NA", ""})


def _present(value):
    return (value or "").strip() not in _ABSENT


def _strongest(*evalues):
    """The smallest of the E-values that are present, as the string pharokka wrote."""
    return min((v.strip() for v in evalues if _present(v)), key=float, default="")


def parse_merged(text):
    """Hits from pharokka's `*_full_merged_output.tsv`, as a list of dicts.

    Each dict has: query, source ('pharokka', 'card' or 'vfdb'), family_id, label,
    category, evalue, query_length. No coordinates - see the module docstring.
    """
    hits = []
    for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
        query = row["ID"]
        length = int(row["length"])

        # --- the phage family --------------------------------------------------------
        phrog = (row.get("phrog") or "").strip()
        if phrog.isdigit():
            hits.append({
                "query": query,
                "source": "pharokka",
                "family_id": f"phrog_{phrog}",
                "label": (row.get("annot") or "").strip(),
                "category": (row.get("category") or "").strip(),
                "evalue": _strongest(row.get("mmseqs_eVal"), row.get("pyhmmer_evalue")),
                "query_length": length,
            })

        # --- CARD: antimicrobial resistance ------------------------------------------
        if _present(row.get("CARD_hit")):
            hits.append({
                "query": query,
                "source": "card",
                "family_id": (row.get("ARO_Accession") or "").strip(),
                "label": (row.get("AMR_Gene_Family") or row.get("CARD_short_name")
                          or "").strip(),
                "category": (row.get("Resistance_Mechanism") or "").strip(),
                "evalue": (row.get("CARD_eVal") or "").strip(),
                "query_length": length,
            })

        # --- VFDB: virulence factors -------------------------------------------------
        if _present(row.get("vfdb_hit")):
            hits.append({
                "query": query,
                "source": "vfdb",
                "family_id": (row.get("vfdb_hit") or "").strip(),
                "label": (row.get("vfdb_description") or row.get("vfdb_short_name")
                          or "").strip(),
                "category": "",
                "evalue": (row.get("vfdb_eVal") or "").strip(),
                "query_length": length,
            })
    return hits
