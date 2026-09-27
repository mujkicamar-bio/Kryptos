"""Pfam family metadata, from the release file rather than from the search output.

hmmsearch --domtblout carries only the family name and accession. The description, type
and clan come from Pfam-A.hmm.dat, which ships with the release: the description gives a
family name a matchable meaning, and the clan groups families that are homologous but too
divergent to align as one.
"""
import gzip
import pathlib

# Stockholm '#=GF <tag> <value>' tags this module keeps, mapped to output key.
_TAGS = {
    "ID": "name",
    "AC": "accession",
    "DE": "description",
    "TP": "type",
    "CL": "clan",
}


def parse_pfam_dat(text):
    """Parse Pfam-A.hmm.dat contents into {family name: metadata}.

    Records are separated by '//'. A family with no '#=GF CL' line belongs to no clan and
    gets an empty string, not a missing key: most families have no clan, so a missing key
    would push a guard into every caller and fail on the common case.
    """
    families = {}
    record = {}
    for line in text.splitlines():
        if line.startswith("//"):
            name = record.pop("name", "")
            if name:
                families[name] = {
                    "accession": record.get("accession", ""),
                    "description": record.get("description", ""),
                    "type": record.get("type", ""),
                    "clan": record.get("clan", ""),
                }
            record = {}
            continue
        if not line.startswith("#=GF "):
            continue
        parts = line[len("#=GF "):].split(None, 1)
        if len(parts) != 2:
            continue
        tag, value = parts
        key = _TAGS.get(tag)
        if key:
            record[key] = value.strip()
    return families


def load(path):
    """Parse Pfam-A.hmm.dat from a path; a `.gz` suffix means gzipped."""
    path = pathlib.Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            return parse_pfam_dat(fh.read())
    return parse_pfam_dat(path.read_text())
