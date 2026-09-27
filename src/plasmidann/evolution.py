"""S7: evolutionary evidence that a dark ORF is a real protein.

WHY dN/dS IS THE STRONGEST FILTER AVAILABLE

Everything else in the pipeline says what is ABSENT: no Pfam hit, no Swiss-Prot hit, no nr
hit. Absence of evidence is exactly what a spurious ORF also produces, and it is why FESNov
discarded 94.3% of its 7,052,473 novel clusters.

dN/dS says something POSITIVE. If replacement changes are suppressed relative to silent
ones, selection is acting on the protein sequence - which means there IS a protein
sequence. It converts "nobody has named it" into "evolution is paying to keep it".

This is also the filter that catches the artefact class nothing else catches. A shadow ORF
on the reverse-complement strand of a real gene is conserved, multi-species, and passes
every absence-based test - but the selection acting on it is acting on the gene on the
OTHER strand, and its own dN/dS reflects that only weakly.

WHY NEI-GOJOBORI COUNTING, AND WHY IT IS THE ONLY ESTIMATOR HERE

Counting needs no tree, no optimiser and no external process. It runs over hundreds of
thousands of small families in-process, is fully unit-testable, and reports a status rather
than a silence when it cannot measure. It is the pipeline's only dN/dS estimate; a codon
model per family was too costly at this scale (docs/PIPELINE_CODE.md), so `dnds_status`
and the pair count carry the qualification of each estimate.
"""
import itertools
import math

# Standard genetic code, table 11 (bacterial). '*' is a stop codon.
GENETIC_CODE = {}
_BASES = "TCAG"
_AAS = ("FFLLSSSSYY**CC*W" "LLLLPPPPHHQQRRRR"
        "IIIMTTTTNNKKSSRR" "VVVVAAAADDEEGGGG")
for _i, _b1 in enumerate(_BASES):
    for _j, _b2 in enumerate(_BASES):
        for _k, _b3 in enumerate(_BASES):
            GENETIC_CODE[_b1 + _b2 + _b3] = _AAS[_i * 16 + _j * 4 + _k]


def synonymous_sites(codon):
    """Split a codon's three positions into synonymous and nonsynonymous fractions.

    For each of the three positions, count how many of the three possible substitutions
    leave the amino acid unchanged; that fraction is the position's synonymous site count.
    A fourfold-degenerate third position contributes 1.0, a first position of most codons
    contributes 0.

    Returns (synonymous, nonsynonymous), summing to 3.0.
    """
    aa = GENETIC_CODE.get(codon)
    if aa is None:
        return 0.0, 0.0
    syn = 0.0
    for pos in range(3):
        for base in "TCAG":
            if base == codon[pos]:
                continue
            mutated = codon[:pos] + base + codon[pos + 1:]
            if GENETIC_CODE.get(mutated) == aa:
                syn += 1 / 3
    return syn, 3.0 - syn


def codon_differences(a, b):
    """Synonymous and nonsynonymous differences between two codons.

    Where two codons differ at more than one position the true mutational path is unknown,
    so all paths are averaged - the standard Nei-Gojobori treatment. Averaging rather than
    picking the parsimonious path avoids a systematic bias toward synonymous change.
    """
    if a == b:
        return 0, 0
    aa_a, aa_b = GENETIC_CODE.get(a), GENETIC_CODE.get(b)
    if aa_a is None or aa_b is None:
        return 0, 0

    positions = [i for i in range(3) if a[i] != b[i]]
    if len(positions) == 1:
        return (1, 0) if aa_a == aa_b else (0, 1)

    # Average over every order in which the differing positions could have changed.
    syn_total = non_total = 0
    paths = 0
    for order in itertools.permutations(positions):
        current = a
        syn = non = 0
        broken = False
        for pos in order:
            nxt = current[:pos] + b[pos] + current[pos + 1:]
            if GENETIC_CODE.get(nxt) == "*":
                # A path through a stop codon is not a path a real lineage took.
                broken = True
                break
            if GENETIC_CODE[nxt] == GENETIC_CODE[current]:
                syn += 1
            else:
                non += 1
            current = nxt
        if not broken:
            syn_total += syn
            non_total += non
            paths += 1
    if not paths:
        return 0, 0
    return syn_total / paths, non_total / paths


def _jukes_cantor(p):
    """Correct an observed proportion of differences for multiple hits at the same site.

    Returns None where the proportion is at or beyond the correction's domain (p >= 0.75),
    which means the sequences are too diverged for the estimate to mean anything.
    """
    if p <= 0:
        return 0.0
    if p >= 0.75:
        return None
    return -0.75 * math.log(1 - (4 / 3) * p)


ABSOLUTE_MIN_CODONS = 3


def dnds_detail(seq_a, seq_b, min_codons=ABSOLUTE_MIN_CODONS):
    """Nei-Gojobori dN/dS with a status code explaining any absent estimate.

    `min_codons` is the number of USABLE codons - gaps, ambiguity and stop codons excluded -
    below which no estimate is returned; S7b passes evolution.min_codons from config. The
    default is the arithmetic floor, three.

    Returns (value, status). Status is one of:

      MEASURED       a usable estimate; value is a float or inf
      NO_DIVERGENCE  the sequences are identical, so there is nothing to measure
      TOO_SHORT      too few usable codons once gaps, ambiguity and stops are excluded
      SATURATED      beyond the Jukes-Cantor domain (pS >= 0.75); no correction exists

    Without a value the status says why. NO_DIVERGENCE (identical sequences, common in
    conserved or clonally redundant families) is not evidence of neutral evolution, and
    evidence.reality_lines does not fire purifying_selection on it; nor is it dN/dS = 0,
    which would read as maximal purifying selection.
    """
    if len(seq_a) != len(seq_b):
        raise ValueError(f"aligned sequences differ in length: {len(seq_a)} vs {len(seq_b)}")
    if len(seq_a) % 3:
        raise ValueError(f"length {len(seq_a)} is not a whole number of codons")

    syn_sites = non_sites = 0.0
    syn_diff = non_diff = 0.0
    usable = 0

    for i in range(0, len(seq_a), 3):
        ca, cb = seq_a[i:i + 3].upper(), seq_b[i:i + 3].upper()
        # Gaps and ambiguity are MISSING DATA, not evidence of conservation. Counting a
        # gapped column as identical would inflate apparent purifying selection exactly
        # where the alignment is least trustworthy.
        if ca not in GENETIC_CODE or cb not in GENETIC_CODE:
            continue
        if GENETIC_CODE[ca] == "*" or GENETIC_CODE[cb] == "*":
            continue
        usable += 1
        sa, na = synonymous_sites(ca)
        sb, nb = synonymous_sites(cb)
        syn_sites += (sa + sb) / 2
        non_sites += (na + nb) / 2
        sd, nd = codon_differences(ca, cb)
        syn_diff += sd
        non_diff += nd

    # Three codons is the arithmetic floor; min_codons is the configured floor for a
    # trustworthy estimate. Below either there is no value, and the status is TOO_SHORT.
    if usable < max(min_codons, ABSOLUTE_MIN_CODONS):
        return None, "TOO_SHORT"
    if syn_diff == 0 and non_diff == 0:
        return None, "NO_DIVERGENCE"
    if not syn_sites or not non_sites:
        return None, "TOO_SHORT"

    ps, pn = syn_diff / syn_sites, non_diff / non_sites
    ds, dn = _jukes_cantor(ps), _jukes_cantor(pn)
    if ds is None or dn is None:
        return None, "SATURATED"
    if ds == 0:
        # Every observed change is a replacement: evidence AGAINST a conserved protein.
        return (float("inf"), "MEASURED") if dn > 0 else (None, "NO_DIVERGENCE")
    return round(dn / ds, 4), "MEASURED"


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
