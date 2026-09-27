"""S7b: evolutionary evidence per dark family.

Everything upstream reports what is ABSENT - no Pfam hit, no Swiss-Prot hit, no nr hit.
Absence is exactly what a spurious ORF also produces, which is why FESNov discarded 94.3%
of its 7,052,473 novel clusters.

This stage reports something POSITIVE:

  dN/dS < 0.5        selection is suppressing replacement changes, so there IS a protein
  RNAcode P < 0.05   the region carries coding signal independent of the gene caller

WHY RNAcode SCORES BOTH STRANDS HERE

RNAcode evaluates all six reading frames and reports the sense and antisense signal
separately. For this project that is the point, not a detail. A shadow ORF is the reverse
complement of a real gene, so the coding signal on its ANTISENSE strand should be the
stronger of the two - and a shadow ORF is exactly the artefact class that survives every
absence-based test in this pipeline. Recording only the sense P would discard the one
number that separates them.

The two are reported, never combined into a verdict: `rnacode_p`, `rnacode_p_antisense`,
and `coding_signal` as the declared threshold applied to the sense P.

dN/dS is also the only filter that catches the artefact class nothing else does. A shadow
ORF on the reverse-complement strand of a real gene is conserved, multi-species, and passes
every absence-based test - but the selection acting on that locus is acting on the gene on
the OTHER strand.

The consensus re-check matters because a protein can miss every per-sequence threshold
while its family is collectively recognisable; this removed 6.5% of clusters in
Pavlopoulos et al.

TWO MEMBER SETS PER FAMILY

Every measurement is made over all of the family's dark members, and again, in columns
prefixed small_, over the dark members that occur on small plasmids (dark_families.tsv
small_members). The consensus is built from all members. Families are processed in
parallel, one alignment per process (mafft --thread 1).
"""
import csv
import multiprocessing
import pathlib

import _ctx  # noqa: F401

from plasmidann import scratch
from plasmidann.evolution_worker import COLS, configure, measure
from plasmidann.fasta import iter_fasta

cfg = snakemake.params.evolution
proteins = dict(iter_fasta([snakemake.input.faa]))
cds = dict(iter_fasta([snakemake.input.cds]))

families = []
with open(snakemake.input.families, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        families.append(r)

tmpdir = str(scratch.scratch_dir(pathlib.Path(snakemake.output[0]).parent))
# The consensus carries the family's shared signal and is re-searched at S7c.
consensus_out = open(snakemake.output.consensus, "w")
jobs = []
# A family whose dark members all sit on small plasmids has one member set, not two: the
# small_ measurement is then the same alignment, and is copied rather than recomputed.
same = {}
for fam in families:
    members = fam["members"].split(",")
    small = [m for m in fam.get("small_members", "").split(",") if m]
    same[fam["family_id"]] = set(small) == set(members)
    jobs.append((fam["family_id"], "all", members))
    if not same[fam["family_id"]]:
        jobs.append((fam["family_id"], "small", small))

n_tested = n_purifying = n_coding = 0
out_cols = COLS + [f"small_{c}" for c in COLS if c != "family_id"]
with open(snakemake.output.tsv, "w", newline="") as out, \
        multiprocessing.Pool(snakemake.threads, initializer=configure,
                             initargs=(proteins, cds, cfg, tmpdir)) as pool:
    w = csv.DictWriter(out, fieldnames=out_cols, delimiter="\t")
    w.writeheader()
    # imap keeps the job order, so each family's "all" result is followed by its "small"
    # when it has a separate one.
    results = pool.imap(measure, jobs, chunksize=4)
    for fam in families:
        fid = fam["family_id"]
        row, family_consensus = next(results)
        small_row = dict(row) if same[fid] else next(results)[0]
        if family_consensus:
            consensus_out.write(f">{fid}\n{family_consensus}\n")
        n_tested += row["dnds_status"] == "MEASURED"
        n_purifying += row["under_purifying_selection"] == 1
        n_coding += row["coding_signal"] == 1
        row.update({f"small_{k}": v for k, v in small_row.items() if k != "family_id"})
        w.writerow(row)

consensus_out.close()

print(f"families={len(families)} with_dnds={n_tested} purifying={n_purifying} "
      f"coding_signal={n_coding} (all dark members; small-plasmid members in small_*)")

# The scratch directory is removed only here, on the ordinary path. A script that raised
# never reaches this line, and its intermediates are what the failure is diagnosed from.
scratch.release(tmpdir)

