"""Pfam family metadata, which the search output does not carry.

hmmsearch --domtblout gives the family NAME and accession. The description and the clan
are in Pfam-A.hmm.dat. Pfam 38.2 ships 30,134
families; the parser must handle every one without a special case.
"""
import gzip

from plasmidann import pfam_meta

SAMPLE = """# STOCKHOLM 1.0
#=GF ID   RepA_N
#=GF AC   PF06970.16
#=GF DE   Replication initiator protein A (RepA) N-terminus
#=GF GA   25.5; 25.5;
#=GF TP   Family
#=GF CL   CL0123
//
# STOCKHOLM 1.0
#=GF ID   MobA_MobL
#=GF AC   PF03389.20
#=GF DE   MobA/MobL family protein
#=GF GA   21.7; 21.7;
#=GF TP   Family
//
"""


def test_a_family_with_a_clan_is_parsed_completely():
    meta = pfam_meta.parse_pfam_dat(SAMPLE)

    assert meta["RepA_N"] == {
        "description": "Replication initiator protein A (RepA) N-terminus",
        "clan": "CL0123",
    }


def test_a_family_without_a_clan_gets_an_empty_clan_not_a_missing_key():
    """Most families belong to no clan. A missing key would make every consumer guard for
    it, and one that forgot would raise on the common case rather than the rare one."""
    meta = pfam_meta.parse_pfam_dat(SAMPLE)

    assert meta["MobA_MobL"]["clan"] == ""
    assert meta["MobA_MobL"]["description"] == "MobA/MobL family protein"


def test_every_family_in_the_sample_is_present():
    assert set(pfam_meta.parse_pfam_dat(SAMPLE)) == {"RepA_N", "MobA_MobL"}


def test_a_gzipped_file_is_read_transparently(tmp_path):
    """Pfam ships the file gzipped and it is stored gzipped, so requiring the caller to
    decompress it first would mean every caller decompressing a 30,134-entry file."""
    path = tmp_path / "Pfam-A.hmm.dat.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(SAMPLE)

    assert pfam_meta.load(path)["RepA_N"]["clan"] == "CL0123"


def test_an_uncompressed_file_is_also_read(tmp_path):
    path = tmp_path / "Pfam-A.hmm.dat"
    path.write_text(SAMPLE)

    assert pfam_meta.load(path)["MobA_MobL"]["description"] == "MobA/MobL family protein"
