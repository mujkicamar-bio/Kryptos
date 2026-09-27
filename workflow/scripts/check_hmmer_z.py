"""Rule check_hmmer_z: check hmmer_z against the number of unique proteins in
unique_proteins.faa; every search depends on this rule.

-Z counts every unique protein, searched or not, so that E-values do not move with the
selection or with PlasmidScope's coverage. The check allows a 2% difference
(cascade.HMMER_Z_TOLERANCE).
"""
import _ctx  # noqa: F401

from plasmidann.cascade import check_hmmer_z

n_unique = sum(1 for l in open(snakemake.input.unique) if l[0] == ">")
check_hmmer_z(declared=snakemake.params.hmmer_z, actual=n_unique)

with open(snakemake.output[0], "w") as out:
    out.write(f"hmmer_z\t{snakemake.params.hmmer_z}\nunique_proteins\t{n_unique}\n")
print(f"hmmer_z {snakemake.params.hmmer_z} confirmed: {n_unique} unique proteins")
