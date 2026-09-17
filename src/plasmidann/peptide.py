"""Peptide physicochemistry for stratum assignment.

THE DECLARED METHOD, AND ITS LIMITS

Membrane and secretion signals are predicted from Kyte-Doolittle hydrophobicity windows and
net charge. Licence-restricted topology predictors are deliberately NOT part of this
pipeline: a stage that cannot be installed from the environment file is a stage that cannot
be reproduced, and every result here has to be reproducible from `envs/` alone.

These calculations are coarse, and the coarseness is real: a hydrophobic soluble core will
sometimes read as a membrane helix. They are therefore used to ASSIGN STRATA and to RANK,
never to exclude, and every output row records `topology_method` so a reader knows exactly
what produced the call. A mis-assigned stratum costs a wasted barcode; it never removes a
candidate from consideration.
"""

# Kyte-Doolittle hydropathy. Positive is hydrophobic.
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

# Charge at physiological pH. Histidine is fractional: its pKa is near 6.0, so at pH 7.4
# it is only about a tenth protonated. Counting it as a full positive would misclassify
# His-rich proteins as antimicrobial peptides.
CHARGE = {"K": 1.0, "R": 1.0, "H": 0.1, "D": -1.0, "E": -1.0}

# An antimicrobial peptide is short, cationic and amphipathic. Thresholds follow the
# properties of characterised bacterial AMPs; the 36-residue peptide FESNov validated sits
# comfortably inside them.
AMP_MAX_LENGTH = 100
AMP_MIN_CHARGE = 2.0
AMP_MIN_HYDROPHOBIC_FRACTION = 0.3

# A transmembrane helix is roughly 19 residues of sustained hydrophobicity.
TM_WINDOW = 19
TM_THRESHOLD = 1.6


def net_charge(seq):
    """Net charge at physiological pH, from the standard ionisable residues."""
    return round(sum(CHARGE.get(aa, 0.0) for aa in seq.upper()), 4)


def gravy(seq):
    """Grand average of hydropathy: the mean Kyte-Doolittle value over the sequence."""
    if not seq:
        return 0.0
    vals = [KD.get(aa, 0.0) for aa in seq.upper()]
    return round(sum(vals) / len(vals), 4)


def max_hydrophobic_window(seq, window=TM_WINDOW):
    """Highest mean hydropathy over any window of `window` residues.

    The signature of a transmembrane helix. Sequences shorter than the window are scored
    over their whole length rather than returning nothing, so short hydrophobic peptides
    are not silently missed.
    """
    s = seq.upper()
    if not s:
        return 0.0
    vals = [KD.get(aa, 0.0) for aa in s]
    if len(vals) <= window:
        return round(sum(vals) / len(vals), 4)
    best = None
    running = sum(vals[:window])
    best = running
    for i in range(window, len(vals)):
        running += vals[i] - vals[i - window]
        best = max(best, running)
    return round(best / window, 4)


def has_tm_helix(seq, window=TM_WINDOW, threshold=TM_THRESHOLD):
    """Whether the sequence contains a plausible transmembrane segment.

    Reported alongside `topology_method` so a downstream reader knows this came from a
    hydrophobicity window rather than a trained topology model.
    """
    return max_hydrophobic_window(seq, window) >= threshold


def is_cationic_amphipathic(seq):
    """Whether a peptide looks like an antimicrobial peptide.

    Short, net positive, and carrying enough hydrophobic residues to insert into a
    membrane. This is the NOVOQR9B class - a 36-residue peptide found beside antibiotic
    resistance genes that inhibited every Gram-positive species tested - and it is the
    stratum with the cheapest synthesis and the clearest optical readout, since membrane
    disruption and lysis are visible as morphology alone.
    """
    s = seq.upper()
    if not s or len(s) > AMP_MAX_LENGTH:
        return False
    if net_charge(s) < AMP_MIN_CHARGE:
        return False
    hydrophobic = sum(1 for aa in s if KD.get(aa, 0.0) > 0)
    return hydrophobic / len(s) >= AMP_MIN_HYDROPHOBIC_FRACTION
