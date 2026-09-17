"""Pfam family metadata, from the release file rather than from the search output.

hmmsearch --domtblout carries the family name (field 3) and the family accession
(field 4), and nothing else about the family. Its "description of target" column describes
the PROTEIN, not the profile, and is empty for our FASTA.

The description, the type and the clan come from Pfam-A.hmm.dat, which ships with the
release and is already on disk at data/refs/pfam/Pfam-A.hmm.dat.gz for Pfam 38.2. They are
the difference between a hit reading 'RepA_N' and a hit reading 'Replication initiator
protein A (RepA) N-terminus, Family, clan CL0123'.

The description is what makes a functional grouping possible at all: a family NAME carries
no matchable word, and a hand-written list of names cannot reach the scope of a
30,134-family database. Measured on this release, 67 family descriptions mention
replication where a curated list named 16, and 42 mention conjugation where it named 15.

The clan matters beyond readability: clans group families that are homologous but too
divergent to align as one, so two proteins hitting different families of one clan share an
origin. A grouping built on family names alone would treat them as unrelated.
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
    """Parse Pfam-A.hmm.dat from a path, gzipped or not.

    Decided by the suffix rather than by sniffing the magic bytes: the release file name is
    fixed, and a wrong guess here would be reported as a parse failure on a 30,134-entry
    file rather than as the wrong file being passed.
    """
    path = pathlib.Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            return parse_pfam_dat(fh.read())
    return parse_pfam_dat(path.read_text())
