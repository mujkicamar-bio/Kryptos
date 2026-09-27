"""S7: evolutionary evidence that a dark ORF is a real protein.

Pairwise dN, dS and omega = dN/dS by the Yang and Nielsen (2000, Mol Biol Evol 17:32)
method, computed by yn00 of PAML 4.10.7 (Yang 2007, Mol Biol Evol 24:1586) run as is on a
codon alignment projected from the protein alignment; and the consensus of a protein
alignment. One yn00 call per family computes every pair of its aligned members.
"""
import pathlib
import subprocess

YN00_COLS = ["S", "N", "t", "kappa", "omega", "dN", "dN_SE", "dS", "dS_SE"]


def yn00(codon_aln, workdir):
    """Run yn00 on a codon alignment; return (codons, pairs), or None if yn00 failed.

    `codon_aln` maps member -> aligned codon sequence. `codons` is the alignment length
    yn00 reports (ls), after it removes every codon column that carries a gap or an
    ambiguous base in any sequence. `pairs` holds one dict per pair from yn00's
    "(B) Yang & Nielsen (2000) method" table: seq1, seq2 and YN00_COLS, as the strings
    yn00 printed (omega 99.0000 where yn00 finds dS = 0, nan where it cannot estimate).

    icode 0 is the standard code, whose codon assignments equal those of table 11. An
    in-frame TGA inside a coding sequence occurs only in a gene called under table 4
    (TGA = Trp), for which icode 3 (Mycoplasma/Spiroplasma code) is used; yn00 stops with
    an error on a stop codon under the wrong code.
    """
    workdir = pathlib.Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    names = list(codon_aln)
    seqs = [codon_aln[n].upper() for n in names]
    tga = any(s[i:i + 3] == "TGA" for s in seqs for i in range(0, len(s), 3))
    icode = 3 if tga else 0
    with open(workdir / "aln.phy", "w") as fh:
        fh.write(f"{len(seqs)} {len(seqs[0])}\n")
        for i, s in enumerate(seqs):
            fh.write(f"s{i}  {s}\n")
    (workdir / "yn00.ctl").write_text(
        f"seqfile = aln.phy\noutfile = yn00.out\nverbose = 0\nicode = {icode}\n"
        "weighting = 0\ncommonf3x4 = 0\n")
    proc = subprocess.run(["yn00", "yn00.ctl"], cwd=workdir, capture_output=True, text=True)
    out = workdir / "yn00.out"
    if proc.returncode != 0 or not out.exists():
        return None
    return parse_yn00(out.read_text(), names)


def parse_yn00(text, names):
    """(codons, pairs) from a yn00 output file; None if it holds no Yang-Nielsen table.

    `names` are the members in alignment order; yn00 numbers them from 1.
    """
    codons = None
    pairs, in_table = [], False
    for line in text.splitlines():
        if line.startswith("ns =") and "ls =" in line:
            codons = int(line.split("ls =")[1])
        elif line.startswith("(B) Yang & Nielsen"):
            in_table = True
        elif line.startswith("(C)"):
            in_table = False
        elif in_table:
            f = line.split()
            # i j S N t kappa omega dN +- SE dS +- SE
            if len(f) == 13 and f[0].isdigit() and f[1].isdigit():
                values = [f[k] for k in (2, 3, 4, 5, 6, 7, 9, 10, 12)]
                pairs.append({"seq1": names[int(f[0]) - 1], "seq2": names[int(f[1]) - 1],
                              **dict(zip(YN00_COLS, values))})
    if codons is None or not pairs:
        return None
    return codons, pairs


def back_translate(aligned_protein, cds):
    """Turn a protein alignment column-for-column into a codon alignment.

    mafft aligns proteins, because protein alignment is far more reliable at the identities
    these families show. dN/dS needs codons. Each protein gap becomes exactly three
    nucleotide gaps - never one or two, which would shift the reading frame and turn the
    whole estimate into noise.
    """
    residues = len(aligned_protein.replace("-", ""))
    if len(cds) < residues * 3:
        raise ValueError(
            f"CDS has {len(cds)} nt, too short for {residues} aligned residues "
            f"({residues * 3} nt needed)")
    out, pos = [], 0
    for aa in aligned_protein:
        if aa == "-":
            out.append("---")
        else:
            out.append(cds[pos:pos + 3])
            pos += 3
    return "".join(out)


# ---------------------------------------------------------------------------------
# The family-consensus re-check
# ---------------------------------------------------------------------------------

def consensus(alignment, max_gap_fraction=0.5):
    """The commonest residue per alignment column, as one ungapped sequence.

    A protein can miss every per-sequence threshold while its family is collectively
    recognisable; Pavlopoulos et al. removed 6.5% of their clusters by searching the family
    consensus back against the reference databases (S7c).

    `alignment` maps name -> aligned sequence, all the same length.

    A column that is gap in more than `max_gap_fraction` of members is DROPPED. Such a
    column is an insertion in a minority of members rather than part of what the family
    shares, and keeping it would build a consensus that no member actually carries - so the
    re-search would be of a sequence that does not exist.

    Ties are broken by residue, alphabetically, rather than by whichever the counter
    happened to see first. Without that the consensus, and therefore the re-check verdict,
    would change between runs on identical input.
    """
    if not alignment:
        return ""
    sequences = list(alignment.values())
    length = len(sequences[0])
    if any(len(s) != length for s in sequences):
        raise ValueError("consensus needs an alignment; sequences differ in length")

    out = []
    for i in range(length):
        column = [s[i].upper() for s in sequences]
        residues = [c for c in column if c not in "-."]
        if len(residues) <= len(column) * max_gap_fraction:
            continue
        counts = {}
        for c in residues:
            counts[c] = counts.get(c, 0) + 1
        out.append(min(sorted(counts), key=lambda c: (-counts[c], c)))
    return "".join(out)
