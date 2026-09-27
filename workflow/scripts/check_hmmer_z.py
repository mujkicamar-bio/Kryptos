"""S2z: confirm that the declared -Z still describes the protein set, as soon as it exists.

hmmer_z is the number of unique proteins. Change the ORF set - S1 origin repair and S0
terminal-repeat trimming both do - and the count moves. A stale -Z silently rescales every
E-value in the run, which is precisely the failure -Z was introduced to prevent.

The check runs directly after dereplication, and the artefact screen and every tier depend
on it, so nothing searches with a -Z that was not confirmed.

-Z counts EVERY unique protein, including those PlasmidScope annotated and those the
selection does not search: an E-value then means what it would if the whole collection had
been searched, and does not move with PlasmidScope's coverage or with the selection. The
larger -Z is the conservative choice - a smaller one makes weak hits significant, and each
removes a protein from the dark set.
"""
import _ctx  # noqa: F401

from plasmidann.cascade import check_hmmer_z

n_unique = sum(1 for l in open(snakemake.input.unique) if l[0] == ">")
check_hmmer_z(declared=snakemake.params.hmmer_z, actual=n_unique)

with open(snakemake.output[0], "w") as out:
    out.write(f"hmmer_z\t{snakemake.params.hmmer_z}\nunique_proteins\t{n_unique}\n")
print(f"hmmer_z {snakemake.params.hmmer_z} confirmed: {n_unique} unique proteins")
