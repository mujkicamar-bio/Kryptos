"""S0: the analysis set - which plasmids are in scope - as an id list and as sequence.

The scope is decided ONCE, here, from the master table's locked exclusion. Every stage
that needs sequence reads the FASTA written by this script rather than the configured
input, so the configured input may be the whole working set: a record that is not in the
analysis set never reaches gene calling, lineage clustering or the feature files.

The study is about SMALL plasmids (size_bp below input.max_plasmid_size_bp), but every
plasmid is in the analysis set: large-plasmid proteins are clustered with the small ones,
and those in families holding a small-plasmid protein are annotated too (S2s,
cascade_selection). The small ones are listed separately; that list is what "on a small
plasmid" means everywhere downstream.

A record whose topology is circular (darkorf.circular.is_circular) and whose two ends
carry the same >= orf.min_terminal_repeat_bp bases is written with the last copy removed,
so the FASTA holds the molecule once. Left in, the gene caller would join the ends through
both copies and translate a frameshifted protein across the junction, and every other
sequence-reading stage would see the repeated region twice. Every trimmed record is listed
in terminal_repeats.tsv with the bases removed.
"""
import csv

import _ctx  # noqa: F401

from darkorf.circular import is_circular, terminal_repeat_length
from plasmidann.fasta import iter_fasta

exclude = set(snakemake.params.exclude)
max_size = snakemake.params.max_size_bp
min_repeat = snakemake.params.min_terminal_repeat_bp
keep, small, circular = set(), set(), set()
with open(snakemake.input.master, newline="") as fh, open(snakemake.output.ids, "w") as out, \
        open(snakemake.output.small_ids, "w") as small_out:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["hab_top"] not in exclude:
            keep.add(row["plasmid_id"])
            if is_circular(row.get("topology")):
                circular.add(row["plasmid_id"])
            out.write(row["plasmid_id"] + "\n")
            if int(row["size_bp"]) < max_size:
                small.add(row["plasmid_id"])
                small_out.write(row["plasmid_id"] + "\n")

n_seen = n_written = n_trimmed = 0
with open(snakemake.output.fasta, "w") as out, \
        open(snakemake.output.repeats, "w") as rep, \
        open(snakemake.output.lengths, "w") as lengths:
    rep.write("plasmid_id\trecord_bp\trepeat_bp\tmolecule_bp\n")
    lengths.write("plasmid_id\tlength_bp\n")
    for pid, seq in iter_fasta([snakemake.input.fasta]):
        n_seen += 1
        if pid in keep:
            n_written += 1
            k = terminal_repeat_length(seq.upper(), min_repeat) if pid in circular else 0
            if k:
                n_trimmed += 1
                rep.write(f"{pid}\t{len(seq)}\t{k}\t{len(seq) - k}\n")
                seq = seq[:-k]
            out.write(f">{pid}\n{seq}\n")
            lengths.write(f"{pid}\t{len(seq)}\n")

# An in-scope set with no sequence is never a legitimate result: every later stage would
# succeed while writing well-formed empty tables.
if n_written == 0:
    raise SystemExit(
        f"analysis_set: none of the {n_seen} records in {snakemake.input.fasta} is in the "
        f"analysis set of {len(keep)} plasmids. Check that the FASTA identifiers are the "
        "master table's plasmid_id.")
print(f"analysis set: {len(keep)} plasmids in scope, {len(small)} of them small "
      f"(size_bp < {max_size}); {n_written} of {n_seen} FASTA records retained; "
      f"{n_trimmed} circular records carried a terminal repeat of >= {min_repeat} bp, "
      "one copy removed")
