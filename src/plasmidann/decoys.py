"""Negative controls: decoy sequences that must NOT be annotated (spec section 58.2).

WHAT THEY TEST

Positive controls ask "can the pipeline still recover known biology?". Negative controls
ask the opposite question, and it is the one this project cannot answer without them: does
the pipeline avoid treating a non-biological sequence as a convincing biological protein?

That matters here more than in most pipelines, because the deliverable is the DARK set. A
dark protein is one nothing could name. If the cascade can be induced to name a sequence
that is not a protein at all, then the complement - the set nothing named - is not what it
claims to be, and every downstream statement about dark proteins inherits the error.

Neither control establishes the expected SIZE of the dark population (spec section 59).
They test the two directions of the classifier, not its calibration.

THE TWO CONSTRUCTIONS

Both are derived from real plasmid CDS, so they carry the composition of the collection
rather than of a uniform random model. A decoy drawn from a uniform model is too easy: it
fails to resemble anything for reasons that have nothing to do with the cascade.

  shuffled            the residues of a real protein, permuted. Length and amino-acid
                      composition are preserved EXACTLY; all order information is
                      destroyed. Any hit to a shuffled sequence is a hit to composition
                      alone, which is what a decoy is for.

  reverse_complement  the real CDS read on the opposite strand and translated. This is the
                      shadow-ORF artifact the QC stage exists to flag (spec section 9.3),
                      manufactured deliberately: it has realistic codon statistics and
                      realistic length, and it is not a protein.

EXPECTED BEHAVIOUR (spec section 58.2)

    record_class = negative_control
    no functional annotation
    DARK, or an explicit decoy/QC status

A decoy that comes out FUNCTIONAL is a false positive of the annotation cascade, and the
rate of those is a measurement this pipeline should report rather than assume.
"""

DECOY_PREFIX = "DECOY_"

# Stop codons in the standard genetic code, as translated by the caller.
_STOP = "*"

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

# The standard genetic code, for translating a reverse-complemented CDS. Written out rather
# than imported from a toolkit so the decoy construction is fully specified here: spec
# section 58.2 requires the exact construction to be fixed and documented.
_CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def shuffled_decoy(sequence, rng):
    """A protein's residues in permuted order.

    Length and amino-acid composition are preserved exactly, so the decoy differs from the
    real protein in ORDER alone. That is the point: a hit to this sequence is a hit to
    composition, and composition is not evidence of function.
    """
    residues = list(sequence)
    rng.shuffle(residues)
    return "".join(residues)


def reverse_complement(dna):
    """The reverse complement of a nucleotide sequence."""
    return dna.translate(_COMPLEMENT)[::-1]


def translate(dna):
    """Translate a nucleotide sequence in frame 1, stopping at the first stop codon.

    Stops are truncated rather than kept: a sequence with an internal '*' is not a protein
    the searches can handle, and every tool would either reject it or treat the stop as an
    unknown residue. Truncating gives the longest real open frame on that strand, which is
    what a shadow ORF actually looks like.
    """
    out = []
    for i in range(0, len(dna) - 2, 3):
        residue = _CODONS.get(dna[i:i + 3].upper(), "X")
        if residue == _STOP:
            break
        out.append(residue)
    return "".join(out)


def reverse_complement_decoy(cds):
    """The protein you get by reading a real CDS on the opposite strand.

    This is the shadow-ORF artifact manufactured on purpose: realistic codon statistics,
    realistic length, and not a protein. Returns '' when the opposite strand has no open
    frame worth searching, which the caller treats as "this CDS cannot make a decoy".
    """
    return translate(reverse_complement(cds))


def build_decoys(cds_records, n_shuffled, n_reverse_complement, rng, min_length=50):
    """Build both decoy classes from real CDS.

    `cds_records` is a sequence of (name, cds_nucleotide_sequence). Returns a list of
    (decoy_id, protein_sequence, decoy_class).

    Each source CDS contributes to at most one class, and the two classes are drawn from
    DISJOINT sets of source sequences. Drawing both from the same CDS would make the two
    controls correlated: a composition quirk that fooled the cascade would appear twice and
    read as two independent failures.

    Decoys shorter than `min_length` are skipped. A 12-residue decoy is not a test of the
    cascade - nothing would name it, and it would pass the negative control for a reason
    that has nothing to do with the pipeline behaving correctly.

    A CDS carrying an ambiguity code is skipped: N translates to X, and a run of X is
    neither a real protein nor a decoy with realistic composition.
    """
    usable = [(name, cds) for name, cds in cds_records
              if len(cds) >= min_length * 3 and set(cds.upper()) <= set("ACGT")]
    rng.shuffle(usable)

    out = []
    source = iter(usable)

    for _ in range(n_shuffled):
        for name, cds in source:
            protein = translate(cds)
            if len(protein) >= min_length:
                out.append((f"{DECOY_PREFIX}shuf_{len(out):05d}",
                            shuffled_decoy(protein, rng), "shuffled"))
                break

    for _ in range(n_reverse_complement):
        for name, cds in source:
            protein = reverse_complement_decoy(cds)
            if len(protein) >= min_length:
                out.append((f"{DECOY_PREFIX}rc_{len(out):05d}", protein,
                            "reverse_complement"))
                break

    return out
