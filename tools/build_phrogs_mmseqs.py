"""Build an MMseqs2 profile database from the PHROGs HH-suite profiles.

WHY MMseqs2 AND NOT hmmsearch

PHROGs ships as HH-suite `.hhm` profiles - 38,880 of them. hmmsearch cannot read that
format at all, and spec section 79 excludes HMM-HMM search ("HMM-HMM and InterPro are
absent unless the pipeline specification is formally revised"), which rules out running
HHsearch over them.

MMseqs2 is the appropriate method and is already in the environment. `mmseqs
convertprofiledb` reads HH-suite HHM format directly and produces a profile database that
`mmseqs search` queries with a sequence set - profile-versus-sequence, not HMM-HMM. Spec
section 18 asks for "a profile-based search against the configured PHROGs database", which
is exactly this.

WHAT THIS SCRIPT DOES

convertprofiledb takes an ffindex database, not a directory of files, so the profiles are
first concatenated into one data file with an index beside it. The ffindex format is:

    <name>.ffdata    the entries, each terminated by a NUL byte
    <name>.ffindex   one line per entry: name<TAB>offset<TAB>length

The PHROG identifier is the file stem (phrog_10000), which becomes the entry name and
therefore the target identifier a hit reports.

THE ANNOTATION TABLE IS SEPARATE

PHROGs distributes a table mapping each PHROG to one of nine functional categories. It is
not in this tree and its host was unreachable when this was written, so the category is not
available. What IS available is the description on each profile's NAME line, which names
the representative protein, and that is extracted here into a lookup the pipeline reads for
`phrog_description`. When the annotation table can be downloaded, add it beside this file
and the category becomes a join rather than a re-run.
"""
import argparse
import pathlib
import subprocess
import sys


def parse_name_line(text):
    """The description from an HH-suite NAME line.

    Format, from the PHROGs profiles:

        NAME  JN699628_p59 !!  !! BONGO_59 !! hypothetical protein !! AER26102.1 !! ...

    The fields are '!!'-separated and the protein description is the fourth. A profile
    whose line does not have that shape returns the line as-is rather than nothing: a
    missing description is invisible in the output, an odd one is not.
    """
    for line in text.splitlines():
        if not line.startswith("NAME"):
            continue
        body = line[4:].strip()
        parts = [p.strip() for p in body.split("!!")]
        if len(parts) >= 4 and parts[3]:
            return parts[3]
        return body
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hhm-dir", required=True, help="directory of .hhm profiles")
    ap.add_argument("--out-prefix", required=True,
                    help="output prefix, e.g. data/refs/phrogs/phrogs_profile")
    ap.add_argument("--mmseqs", default="mmseqs")
    args = ap.parse_args()

    hhm_dir = pathlib.Path(args.hhm_dir)
    profiles = sorted(hhm_dir.glob("*.hhm"))
    if not profiles:
        sys.exit(f"no .hhm profiles in {hhm_dir}")

    out = pathlib.Path(args.out_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    ffdata = out.with_suffix(".ffdata")
    ffindex = out.with_suffix(".ffindex")
    descriptions = out.with_suffix(".descriptions.tsv")

    offset = 0
    with open(ffdata, "wb") as data, open(ffindex, "w") as index, \
            open(descriptions, "w") as desc:
        desc.write("phrog_id\tphrog_description\n")
        for path in profiles:
            name = path.stem
            text = path.read_text()
            # ffindex entries are NUL-terminated; the length recorded INCLUDES that byte,
            # which is what ffindex readers expect. Getting this wrong truncates the last
            # residue of every profile, silently.
            payload = text.encode() + b"\0"
            data.write(payload)
            index.write(f"{name}\t{offset}\t{len(payload)}\n")
            offset += len(payload)
            desc.write(f"{name}\t{parse_name_line(text)}\n")

    print(f"ffindex: {len(profiles)} profiles, {offset / 1e6:.0f} MB -> {ffdata}")

    # convertprofiledb wants the prefix; it reads <prefix> and <prefix>.index, so the
    # ffindex pair is presented under the names it expects.
    staged = out.with_name(out.name + "_hhm")
    staged.write_bytes(ffdata.read_bytes())
    pathlib.Path(str(staged) + ".index").write_text(ffindex.read_text())

    profile_db = out.with_name(out.name + "_db")
    subprocess.run([args.mmseqs, "convertprofiledb", str(staged), str(profile_db)],
                   check=True)
    print(f"profile database: {profile_db}")
    print(f"descriptions:     {descriptions}")


if __name__ == "__main__":
    main()
