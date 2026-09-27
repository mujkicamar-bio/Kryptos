"""The per-family work of rule family_evolution (workflow/scripts/family_evolution.py), in a module.

A process pool can only run a function it can import by name, and a Snakemake script is
not importable - hence this module, with the same configure()-per-worker pattern as
darkorf.genecall. See the script's docstring for what is measured and why.
"""
import math
import shutil
import statistics
import subprocess

from plasmidann.evolution import back_translate, consensus, yn00

COLS = ["family_id", "n_aligned", "dnds_median", "dnds_min", "n_pairs",
        "under_purifying_selection", "dnds_status", "rnacode_p", "rnacode_p_antisense",
        "rnacode_status", "coding_signal", "evidence_note"]

_proteins = _cds = _cfg = _tmpdir = None


def configure(proteins, cds, cfg, tmpdir):
    """Give this worker the sequences, the thresholds and the scratch directory."""
    global _proteins, _cds, _cfg, _tmpdir
    _proteins, _cds, _cfg, _tmpdir = proteins, cds, cfg, tmpdir


def clustal(alignment, path):
    """Write a codon alignment as Clustal W.

    RNAcode reads Clustal W or MAF and nothing else - handed FASTA it prints
    "ERROR: Unknown alignment file format" and EXITS 0, so the caller has to look at what
    it wrote rather than at its status.
    """
    names = list(alignment)
    width = 60
    with open(path, "w") as fh:
        fh.write("CLUSTAL W (1.81) multiple sequence alignment\n\n\n")
        length = len(alignment[names[0]])
        for i in range(0, length, width):
            for n in names:
                fh.write(f"{n[:15]:<16}{alignment[n][i:i + width]}\n")
            fh.write("\n\n")


def rnacode(alignment, path):
    """Best sense and antisense P from RNAcode over a codon alignment.

    Returns (p_sense, p_antisense, status). Tabular columns are
    HSS, strand, frame, length, from, to, name, start, end, score, P.

    NO_OUTPUT when RNAcode did not run: a non-zero exit (including a missing executable),
    or an ERROR line on stdout, which RNAcode prints for an unreadable alignment while
    exiting 0. NO_SIGNAL when it ran and reported no segment. The two differ because no
    signal is evidence against a family and a tool that did not run is not.
    """
    clustal(alignment, path)
    proc = subprocess.run(f"RNAcode -t {path}", shell=True, capture_output=True, text=True)
    best = {"+": None, "-": None}
    for line in proc.stdout.splitlines():
        f = line.split()
        if len(f) < 11 or f[1] not in best:
            continue
        try:
            p = float(f[10])
        except ValueError:
            continue
        if best[f[1]] is None or p < best[f[1]]:
            best[f[1]] = p
    failed = proc.returncode != 0 or any(
        line.startswith("ERROR") for line in proc.stdout.splitlines())
    if failed:
        return None, None, "NO_OUTPUT"
    if best["+"] is None and best["-"] is None:
        return None, None, "NO_SIGNAL"
    return best["+"], best["-"], "MEASURED"


def measure(job):
    """(family_id, member_set, members) -> (row, consensus or None, yn00 pairs).

    One family, one member set. The pairs are yn00's rows, each with family_id and
    member_set added.
    """
    fid, label, members = job
    row = dict.fromkeys(COLS, "")
    row.update(family_id=fid, n_aligned=0, n_pairs=0)

    usable = [m for m in members if m in _cds and m in _proteins]
    if len(usable) < _cfg["min_members_for_dnds"]:
        # Not an error: most families are small. Reported so the shortfall is visible
        # rather than looking like a failed test.
        row["evidence_note"] = "too_few_members"
        row["dnds_status"] = "TOO_FEW_MEMBERS"
        row["rnacode_status"] = "TOO_FEW_MEMBERS"
        return row, None, []

    # Alignment cost grows faster than linearly, so at most evolution.max_members_aligned
    # members are aligned, in file order.
    usable = usable[:_cfg["max_members_aligned"]]
    row["n_aligned"] = len(usable)

    pf = f"{_tmpdir}/{fid}.{label}.faa"
    with open(pf, "w") as fh:
        for m in usable:
            fh.write(f">{m}\n{_proteins[m]}\n")

    # Proteins are aligned, not nucleotides: protein alignment is far more reliable at
    # the 30% identities these families show. The alignment is then projected onto
    # codons, so every protein gap becomes exactly three nucleotide gaps.
    aln = subprocess.run(f"mafft --auto --quiet --thread 1 {pf}",
                         shell=True, capture_output=True, text=True)
    if aln.returncode != 0:
        row["evidence_note"] = "alignment_failed"
        row["dnds_status"] = "ALIGNMENT_FAILED"
        row["rnacode_status"] = "ALIGNMENT_FAILED"
        return row, None, []

    aligned = {}
    name, buf = None, []
    for line in aln.stdout.splitlines():
        if line.startswith(">"):
            if name:
                aligned[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        aligned[name] = "".join(buf)

    # One consensus per family, from the PROTEIN alignment. Written for every family
    # that could be aligned at all, including those with no usable codon pairs: the
    # re-check asks whether the family is collectively recognisable, which is a
    # different question from whether its divergence can be measured.
    family_consensus = consensus(aligned)

    # Every aligned member has the CDS of an ORF encoding its protein, so back_translate
    # raises only when the CDS and protein inputs disagree.
    codon_aln = {m: back_translate(ap, _cds[m]) for m, ap in aligned.items()}

    # yn00 writes several files into its working directory, so each call gets its own
    # directory, removed once read.
    ydir = f"{_tmpdir}/{fid}.{label}.yn00"
    result = yn00(codon_aln, ydir)
    shutil.rmtree(ydir, ignore_errors=True)
    pairs = [{"family_id": fid, "member_set": label, **p} for p in result[1]] if result else []

    # Coding potential, independent of the gene caller and of dN/dS. Both strands: for
    # a shadow ORF the antisense signal is expected to be the stronger one.
    p_sense, p_anti, rc_status = rnacode(codon_aln, f"{_tmpdir}/{fid}.{label}.aln")
    row["rnacode_status"] = rc_status
    if p_sense is not None:
        row["rnacode_p"] = p_sense
        row["coding_signal"] = int(p_sense < _cfg["rnacode_max_p"])
    if p_anti is not None:
        row["rnacode_p_antisense"] = p_anti

    # The family summary is the median and minimum of the omega values yn00 reports,
    # including its 99.0000 for a pair with dS = 0; a pair yn00 reports as nan has no
    # value. evolution.min_codons applies to the alignment length yn00 used.
    omegas = [float(p["omega"]) for p in pairs if not math.isnan(float(p["omega"]))]
    if result is None:
        row["dnds_status"] = "YN00_FAILED"
    elif result[0] < _cfg["min_codons"]:
        row["dnds_status"] = "TOO_SHORT"
    elif not omegas:
        row["dnds_status"] = "NO_ESTIMATE"
    else:
        row["n_pairs"] = len(omegas)
        row["dnds_median"] = round(statistics.median(omegas), 4)
        row["dnds_min"] = round(min(omegas), 4)
        row["under_purifying_selection"] = int(row["dnds_median"] < _cfg["dnds_purifying_max"])
        row["dnds_status"] = "MEASURED"
    return row, family_consensus, pairs
