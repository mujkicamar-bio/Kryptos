"""The per-family work of S7b (workflow/scripts/family_evolution.py), in a module.

A process pool can only run a function it can import by name, and a Snakemake script is
not importable - hence this module, with the same configure()-per-worker pattern as
darkorf.genecall. See the script's docstring for what is measured and why.
"""
import collections
import itertools
import statistics
import subprocess

from plasmidann.evolution import back_translate, consensus, dnds_detail

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
    """(family_id, member_set, members) -> (row, consensus or None). One family, one set."""
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
        return row, None

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
        return row, None

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

    # The status of every pair is counted as well as the values: a family of identical
    # sequences (NO_DIVERGENCE) is not a family measured and found neutral. A pair with
    # dN/dS = inf (dS = 0, dN > 0) is counted in the statuses but not in the median.
    ratios, statuses = [], collections.Counter()
    for a, b in itertools.combinations(sorted(codon_aln), 2):
        r, status = dnds_detail(codon_aln[a], codon_aln[b],
                                min_codons=_cfg["min_codons"])
        statuses[status] += 1
        if r is not None and r != float("inf"):
            ratios.append(r)

    # Coding potential, independent of the gene caller and of dN/dS. Both strands: for
    # a shadow ORF the antisense signal is expected to be the stronger one.
    p_sense, p_anti, rc_status = rnacode(codon_aln, f"{_tmpdir}/{fid}.{label}.aln")
    row["rnacode_status"] = rc_status
    if p_sense is not None:
        row["rnacode_p"] = p_sense
        row["coding_signal"] = int(p_sense < _cfg["rnacode_max_p"])
    if p_anti is not None:
        row["rnacode_p_antisense"] = p_anti

    row["n_pairs"] = len(ratios)
    if ratios:
        row["dnds_median"] = round(statistics.median(ratios), 4)
        row["dnds_min"] = round(min(ratios), 4)
        row["under_purifying_selection"] = int(row["dnds_median"] < _cfg["dnds_purifying_max"])
        row["dnds_status"] = "MEASURED"
    else:
        row["evidence_note"] = "no_informative_pairs"
        # The commonest status across the pairs says which kind of absence this is.
        row["dnds_status"] = (statuses.most_common(1)[0][0] if statuses
                              else "NO_INFORMATIVE_PAIRS")
    return row, family_consensus
