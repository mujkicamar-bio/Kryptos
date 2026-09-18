"""The phage tier: pharokka in protein mode (spec section 18).

WHY PHAROKKA

The prokaryotic virus protein families are distributed as HH-suite HHM profiles, which
nothing in this pipeline can search without the HMM-HMM tooling spec section 79 excludes.
Two conversions of those profiles were implemented and abandoned. pharokka redistributes
the same families as a MMseqs2 profile database - searched with the PROFILES AS QUERY,
which is the orientation that works - and as HMMER3 profiles searched with pyhmmer, and
carries a functional category for every family. It is a tool used as its authors ship it.

WHAT THE OUTPUT IS, AND IS NOT

One row per input protein, whether or not anything hit it. A row carries the family, its
annotation and category, and the statistics of whichever of the two searches found it.
Alongside it, hits against CARD (antimicrobial resistance) and VFDB (virulence factors),
which pharokka searches in the same run.

It does NOT carry alignment coordinates. pharokka deletes its raw alignment tables on exit,
unconditionally, so no span survives. The parser reports that honestly: a hit has no start
and no end, and the cascade treats it as a family-level assignment - the families are
whole-protein clusters, not domains - whose completeness is NOT_MEASURED rather than
invented.

TWO THINGS THE REAL OUTPUT DOES THAT THE DOCUMENTATION DOES NOT SAY

Built against pharokka 1.10.1 on 200 proteins from the test set.

  * A protein with no hit carries 'No_MMseqs' in the `phrog` column, not 'No_PHROG': the
    null-fill runs in column order and mmseqs_phrog is filled first. So "has a family" is
    tested as "the column is an integer", which is what a family identifier is, rather
    than against any sentinel string.
  * 90 of 92 hits were found by BOTH searches, each with its own E-value. The tier reports
    the stronger, because that is the evidence the assignment rests on, and the cascade
    ranks hits by E-value.
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
    candidates = []
    for value in evalues:
        if _present(value):
            try:
                candidates.append((float(value), value.strip()))
            except ValueError:
                continue
    return min(candidates)[1] if candidates else ""


def parse_merged(text):
    """Hits from pharokka's `*_full_merged_output.tsv`, as a list of dicts.

    Each dict has: query, source ('pharokka', 'card' or 'vfdb'), family_id, label,
    category, evalue, found_by, query_length. No coordinates - see the module docstring.
    """
    hits = []
    for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
        query = row["ID"]
        try:
            length = int(row.get("length") or 0)
        except ValueError:
            length = 0

        # --- the phage family --------------------------------------------------------
        phrog = (row.get("phrog") or "").strip()
        if phrog.isdigit():
            mm = _present(row.get("mmseqs_eVal"))
            hm = _present(row.get("pyhmmer_evalue"))
            hits.append({
                "query": query,
                "source": "pharokka",
                "family_id": f"phrog_{phrog}",
                "label": (row.get("annot") or "").strip(),
                "category": (row.get("category") or "").strip(),
                "evalue": _strongest(row.get("mmseqs_eVal"), row.get("pyhmmer_evalue")),
                "found_by": "+".join(n for n, p in (("mmseqs", mm), ("pyhmmer", hm)) if p),
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
                "found_by": "mmseqs",
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
                "found_by": "mmseqs",
                "query_length": length,
            })
    return hits
