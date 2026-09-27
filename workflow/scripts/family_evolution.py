"""Rule family_evolution: evolutionary evidence per dark family (plasmidann.evolution,
plasmidann.evolution_worker).

For each family: the median of the pairwise omega = dN/dS values that PAML yn00 reports
(Yang and Nielsen 2000; a median below evolution.dnds_purifying_max indicates purifying
selection on the protein), with every yn00 pair row written as given to `pairs`; and the
RNAcode coding signal, sense and antisense kept apart (`rnacode_p`, `rnacode_p_antisense`;
`coding_signal` applies evolution.rnacode_max_p to the sense P), because a shadow ORF, the
reverse complement of a real gene, is expected to show the stronger signal on its
antisense strand. The stage also writes one consensus per family from the protein
alignment, which rule consensus_recheck searches back against Pfam.

Every measurement is made over all of the family's dark members and again, in columns
prefixed small_, over the dark members on small plasmids (dark_families.tsv
small_members). Alignments use at most evolution.max_members_aligned members with a
recovered CDS, and the consensus is built from those aligned members. Families are
processed in parallel, one alignment per process (mafft --thread 1).
"""
import csv
import gzip
import multiprocessing
import pathlib

import _ctx  # noqa: F401

from plasmidann import scratch
from plasmidann.evolution import YN00_COLS
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
# The consensus carries the family's shared signal and is re-searched by consensus_recheck.
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
        gzip.open(snakemake.output.pairs, "wt", newline="") as pairs_out, \
        multiprocessing.Pool(snakemake.threads, initializer=configure,
                             initargs=(proteins, cds, cfg, tmpdir)) as pool:
    w = csv.DictWriter(out, fieldnames=out_cols, delimiter="\t")
    w.writeheader()
    pw = csv.DictWriter(pairs_out, delimiter="\t",
                        fieldnames=["family_id", "member_set", "seq1", "seq2"] + YN00_COLS)
    pw.writeheader()
    # imap keeps the job order, so each family's "all" result is followed by its "small"
    # when it has a separate one.
    results = pool.imap(measure, jobs, chunksize=4)
    for fam in families:
        fid = fam["family_id"]
        row, family_consensus, pairs = next(results)
        pw.writerows(pairs)
        if same[fid]:
            small_row = dict(row)
        else:
            small_row, _, small_pairs = next(results)
            pw.writerows(small_pairs)
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

