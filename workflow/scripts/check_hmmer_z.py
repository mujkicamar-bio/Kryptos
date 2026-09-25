"""S2z: confirm that the declared -Z still describes the protein set, as soon as it exists.

hmmer_z is DERIVED from the analysis set: every unique protein plus the spiked controls and
decoys. Change the ORF set - S1 origin repair and S0 terminal-repeat trimming both do - and
the count moves. A stale -Z silently rescales every E-value in the run, which is precisely
the failure -Z was introduced to prevent.

This check used to run in sweep_cohort, after protein clustering and the cascade selection:
hours into a production run, and after the artefact screen had already searched with the
stale value. It now runs directly after dereplication, and the artefact screen and every
tier depend on it, so nothing searches with a -Z that was not confirmed.

-Z counts EVERY unique protein, including those PlasmidScope annotated and those the
selection does not search: an E-value then means what it would if the whole collection had
been searched, and does not move with PlasmidScope's coverage or with the selection. The
larger -Z is the conservative choice - a smaller one makes weak hits significant, and each
removes a protein from the dark set. The control count is the configured one (positive n,
shuffled and reverse-complement decoys).
"""
import _ctx  # noqa: F401

from plasmidann.cascade import check_hmmer_z

n_unique = sum(1 for l in open(snakemake.input.unique) if l[0] == ">")
n_controls = snakemake.params.n_controls
check_hmmer_z(declared=snakemake.params.hmmer_z, actual=n_unique + n_controls)

with open(snakemake.output[0], "w") as out:
    out.write(f"hmmer_z\t{snakemake.params.hmmer_z}\n"
              f"unique_proteins\t{n_unique}\ncontrols\t{n_controls}\n")
print(f"hmmer_z {snakemake.params.hmmer_z} confirmed: {n_unique} unique proteins + "
      f"{n_controls} controls")
