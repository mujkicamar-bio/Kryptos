"""Build the PHROGs HH-suite database, as HH-suite documents building a custom one.

WHY HH-SUITE

PHROGs distributes 38,880 HH-suite HHM profiles and documents HH-suite as the way to search
them. Spec section 79 excluded HMM-HMM search; that exclusion was formally revised for this
database alone on 2026-09-18, and the revision is recorded in PLASMID_ANALYSIS.md
section 22.

The MMseqs2 route was implemented first and abandoned. `mmseqs convertprofiledb` takes each
profile's header from its NAME line, which names the SEED PROTEIN rather than the PHROG, and
it does not assign database keys in ffindex order - so the PHROG identifier could be
recovered neither from the header nor from the key, and spec section 18 lists phrog_id as
required. With HH-suite the identifier is native: ffindex_build keys each entry by its file
name, so phrog_1.hhm is reported as phrog_1.

WHAT THE THREE STEPS PRODUCE

    <prefix>_hhm.ff{data,index}     the profiles themselves
    <prefix>_a3m.ff{data,index}     each profile's alignment, from its SEQ block
    <prefix>_cs219.ff{data,index}   the column-state prefilter hhblits uses

Without cs219, hhblits falls back to scoring every query against all 38,880 profiles. That
is the same search and the same answer, orders of magnitude slower - which is why pre-flight
reports a missing cs219 as an incomplete database rather than as a warning. Nothing in the
output would reveal it.

WHY THE a3m STEP EXISTS

cstranslate builds cs219 from an ALIGNMENT, and HH-suite 3.3 accepts prf, seq, fas, a2m,
a3m or ca3m - not HHM. PHROGs ships HHM only. The first attempt at this build passed
`-I hhm` and cstranslate rejected all 38,880 entries. But hhmake writes the (filtered)
alignment into each profile's SEQ block as a3m, so the documented input is inside the file:
it is extracted, indexed, and cstranslate runs on that with the documented a3m flags.

The cstranslate flags -x 0.3 -c 4 are HH-suite's own documented values for building a
database. They are not values chosen here.
"""
import argparse
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from plasmidann.phrogs import a3m_from_hhm  # noqa: E402


def run(command, cwd=None):
    print("+ " + " ".join(str(c) for c in command), flush=True)
    subprocess.run([str(c) for c in command], check=True, cwd=cwd)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hhm-dir", required=True, help="directory of PHROGs .hhm profiles")
    ap.add_argument("--out-prefix", required=True,
                    help="output prefix, e.g. data/refs/phrogs/phrogs")
    ap.add_argument("--bin-dir", default="",
                    help="directory holding ffindex_build and cstranslate")
    args = ap.parse_args()

    hhm_dir = pathlib.Path(args.hhm_dir).resolve()
    profiles = sorted(hhm_dir.glob("*.hhm"))
    if not profiles:
        sys.exit(f"no .hhm profiles in {hhm_dir}")
    print(f"profiles: {len(profiles)}")

    prefix = pathlib.Path(args.out_prefix).resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    binary = (lambda name: str(pathlib.Path(args.bin_dir).resolve() / name)) \
        if args.bin_dir else (lambda name: name)

    # ffindex_build is run from INSIDE the profile directory so entry names are bare file
    # names rather than paths. The entry name becomes the target identifier a hit reports,
    # and a path there would put the whole directory into every PHROG id.
    run([binary("ffindex_build"), "-s",
         f"{prefix}_hhm.ffdata", f"{prefix}_hhm.ffindex", "."], cwd=hhm_dir)

    # ffindex_build keys each entry by its file name, suffix included, so the hhm entries
    # come out as 'phrog_1.hhm'. hhblits joins the hhm, a3m and cs219 indexes on the bare
    # entry name, and a hit is reported under it, so the suffix is stripped here: the
    # three indexes then agree and the reported identifier is the PHROG id itself.
    index_path = pathlib.Path(f"{prefix}_hhm.ffindex")
    stripped = []
    for line in index_path.read_text().splitlines():
        name, rest = line.split("\t", 1)
        stripped.append(f"{name[:-4] if name.endswith('.hhm') else name}\t{rest}")
    index_path.write_text("\n".join(stripped) + "\n")

    # The alignments, one a3m per profile, written as an ffindex directly. Entry names
    # must match the hhm entries exactly (phrog_1.hhm -> phrog_1.a3m is NOT a match:
    # hhblits joins the two ffindexes on the bare name), so the stem is used for both.
    print("extracting alignments from the SEQ blocks", flush=True)
    offset = 0
    with open(f"{prefix}_a3m.ffdata", "wb") as data, \
            open(f"{prefix}_a3m.ffindex", "w") as index:
        for path in profiles:
            payload = a3m_from_hhm(path.read_text()).encode() + b"\0"
            data.write(payload)
            index.write(f"{path.stem}\t{offset}\t{len(payload)}\n")
            offset += len(payload)

    run([binary("cstranslate"), "-f", "-x", "0.3", "-c", "4", "-I", "a3m",
         "-i", f"{prefix}_a3m", "-o", f"{prefix}_cs219"])

    n_profiles = sum(1 for _ in open(f"{prefix}_hhm.ffindex"))
    n_prefilter = sum(1 for _ in open(f"{prefix}_cs219.ffindex"))
    print(f"built {prefix}: {n_profiles} profiles, {n_prefilter} prefilter entries")
    if n_profiles != n_prefilter:
        sys.exit(
            f"the prefilter covers {n_prefilter} of {n_profiles} profiles. hhblits scores "
            "a query only against profiles the prefilter knows, so the rest would be "
            "silently unsearchable.")


if __name__ == "__main__":
    main()
