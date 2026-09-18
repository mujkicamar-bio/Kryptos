"""PHROGs: HH-suite profile search over the prokaryotic virus protein families.

WHY HH-SUITE, AND WHY THAT IS A RECORDED DEVIATION

PHROGs distributes 38,880 HH-suite HHM profiles and documents HH-suite as the way to search
them. Spec section 79 excludes HMM-HMM search; using PHROGs as its authors distribute it is
an explicit, narrow exception to that exclusion, for this database alone, and it is recorded
in PLASMID_ANALYSIS.md rather than left as a contradiction between the code and the spec.

The alternative was tried first and abandoned. Converting the profiles to an MMseqs2 profile
database lost the PHROG identifier outright: `mmseqs convertprofiledb` takes each profile's
header from its NAME line, which names the SEED PROTEIN rather than the PHROG, and it does
not assign database keys in ffindex order, so the identifier could not be recovered from the
key either. Section 18 lists phrog_id as required, so that database could not produce the
output the spec asks for. It was also slow enough that a search of 38,880 converted profiles
against 5,304 proteins did not finish.

With HH-suite the identifier is native: it is the ffindex key, which ffindex_build takes from
the profile's file name.

THE SCORE IS THE PROBABILITY

Section 18 requires phrog_score. HH-suite reports Prob, E-value, P-value and Score, and
PROBABILITY is the statistic PHROGs' own annotation work reports and filters on - it is
HH-suite's calibrated estimate that the match is a true homology, which an E-value over a
profile database is not. Both are carried; the probability is what the threshold applies to.

This is why the native .hhr hit table is parsed rather than -blasttab, which is easier to
read but emits no probability at all.
"""
import re

# The hit table's column header. Everything before it is the run's own metadata, and the
# alignments after the table repeat each hit in a different shape - parsing those as hits
# would double-count every match.
_HEADER = re.compile(r"^\s*No\s+Hit\s+Prob\s+E-value")

# One hit row. The fields are column-aligned rather than delimited, and the description is
# free text that may be empty or contain spaces, so the row is matched from BOTH ends:
# the rank and name from the left, the seven numeric fields and two ranges from the right.
# Splitting on whitespace instead shifts every field whenever a description is absent.
_HIT = re.compile(
    r"^\s*(?P<rank>\d+)\s+(?P<name>\S+)\s*(?P<description>.*?)\s+"
    r"(?P<prob>\d+\.\d+)\s+(?P<evalue>\S+)\s+(?P<pvalue>\S+)\s+"
    r"(?P<score>-?\d+\.\d+)\s+(?P<ss>-?\d+\.\d+)\s+(?P<cols>\d+)\s+"
    r"(?P<qstart>\d+)-(?P<qend>\d+)\s+"
    r"(?P<tstart>\d+)-(?P<tend>\d+)\s*\((?P<tlen>\d+)\)\s*$")


def parse_hhr(text, min_probability=None):
    """Hits from one HH-suite .hhr result, as a list of dicts.

    `min_probability` drops hits below an HH-suite probability. It is a parameter and not a
    constant because it governs what enters the annotated set, and the dark set is defined
    as the complement of that - a threshold that decides the deliverable belongs in
    configuration where it can be audited and swept.

    A query that matched nothing gives an empty list. That is the commonest outcome for a
    dark protein and is not a parse failure.
    """
    query = ""
    match_columns = 0
    hits = []
    in_table = False

    for line in text.splitlines():
        if line.startswith("Query "):
            parts = line.split(None, 1)
            query = parts[1].strip() if len(parts) > 1 else ""
            continue
        if line.startswith("Match_columns"):
            try:
                match_columns = int(line.split()[1])
            except (IndexError, ValueError):
                match_columns = 0
            continue
        if _HEADER.match(line):
            in_table = True
            continue
        if not in_table:
            continue
        if not line.strip():
            # The blank line after the table ends it. Without this the alignment blocks
            # below would be scanned, and 'No 1' lines there look enough like hit rows to
            # be worth refusing explicitly.
            in_table = False
            continue

        m = _HIT.match(line)
        if not m:
            continue

        probability = float(m.group("prob"))
        if min_probability is not None and probability < min_probability:
            continue

        qstart, qend = int(m.group("qstart")), int(m.group("qend"))
        hits.append({
            "query": query,
            # ffindex_build keys entries by file name, so the database reports
            # 'phrog_1234.hhm'. The stem is the identifier the spec asks for and the key
            # that joins to the PHROGs annotation table.
            "phrog_id": m.group("name")[:-4] if m.group("name").endswith(".hhm")
                        else m.group("name"),
            "phrog_description": m.group("description").strip(),
            "phrog_probability": probability,
            "phrog_evalue": m.group("evalue"),
            "phrog_score": float(m.group("score")),
            # Against the QUERY's length. The Template HMM range and its parenthesised
            # length belong to the PHROG; using those would give a coverage above 1
            # whenever the profile is shorter than the protein.
            "phrog_query_coverage": (round((qend - qstart + 1) / match_columns, 4)
                                     if match_columns else 0.0),
            "query_start": qstart,
            "query_end": qend,
            "query_length": match_columns,
            "profile_length": int(m.group("tlen")),
        })

    return hits


def read_ffindex(data_path, index_path):
    """Yield (name, text) for every entry of an ffindex pair.

    hhblits_omp writes its .hhr results as an ffindex rather than as one file per query,
    because a per-file layout does not survive a few million queries on a shared
    filesystem. Reading it back is therefore part of the stage rather than an option.

    The format is `<name>\\t<offset>\\t<length>` per index line, and each entry is
    NUL-terminated with the recorded length INCLUDING that byte. Taking the length
    verbatim leaves a stray NUL on the end of every entry - for a .hhr that means the last
    hit row never matches its pattern, which drops one hit per query silently.
    """
    with open(index_path) as fh:
        entries = [l.split("\t") for l in fh.read().splitlines() if l.strip()]
    if not entries:
        return

    with open(data_path, "rb") as fh:
        for fields in entries:
            name, offset, length = fields[0], int(fields[1]), int(fields[2])
            fh.seek(offset)
            yield name, fh.read(length).rstrip(b"\0").decode("utf-8", "replace")
