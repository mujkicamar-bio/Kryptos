"""S7a: recover the nucleotide coding sequence for every dark ORF.

dN/dS needs codons, and the pipeline stores protein only - but orf_index.tsv carries
plasmid_id, start, end, strand and spans_origin, which is enough to recover the CDS from
the analysis-set FASTA written at S0, which is what every stage that needs sequence reads.

Two details that are easy to get wrong and silently corrupt every downstream estimate:

  spans_origin  A gene reconstructed across the cut point of a circular plasmid runs
                start..length THEN 1..end - the GenBank join() convention - so start > end
                and a naive slice returns nothing. These are exactly the genes S1 worked to
                recover, so dropping them here would undo that work.

  strand        A gene on the minus strand must be reverse-complemented before translation.
                Forgetting this yields a sequence that translates to nonsense, which
                dN/dS would happily report a number for.
"""
import csv

import _ctx  # noqa: F401

from plasmidann.fasta import iter_fasta

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


wanted = {l.strip() for l in open(snakemake.input.ids) if l.strip()}

# One representative ORF per unique protein is enough: identical proteins may differ
# synonymously, but the family alignment below draws on many members anyway, and taking
# every ORF would multiply the work by 2.7x for no additional signal.
seq_of_orf = {}
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, ids = line.rstrip("\n").split("\t")
        if sid in wanted:
            seq_of_orf[ids.split(",")[0]] = sid

needed = {}
with open(snakemake.input.index, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        sid = seq_of_orf.get(r["orf_id"])
        if sid is None:
            continue
        needed.setdefault(r["plasmid_id"], []).append((sid, r))

n_written = n_origin = 0
with open(snakemake.output[0], "w") as out:
    for pid, seq in iter_fasta([snakemake.input.fasta]):
        for sid, r in needed.get(pid, ()):
            start, end = int(r["start"]), int(r["end"])
            if r.get("spans_origin") == "1":
                # start..L then 1..end
                nt = seq[start - 1:] + seq[:end]
                n_origin += 1
            else:
                nt = seq[start - 1:end]
            if r["strand"] in ("-1", "-"):
                nt = revcomp(nt)
            if len(nt) < 6:
                continue
            out.write(f">{sid}\n{nt}\n")
            n_written += 1

print(f"CDS written={n_written} of {len(wanted)} dark proteins "
      f"(origin-spanning={n_origin})")
