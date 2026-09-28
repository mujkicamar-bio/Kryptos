"""S0: the analysis set - which plasmids are in scope - as an id list and as sequence.

The scope is decided ONCE, here: a plasmid is in the analysis set when the master table's
locked exclusion keeps it and the configured FASTA holds its sequence. Every stage that
needs sequence reads the FASTA written by this script rather than the configured input,
so the configured input may be the whole working set, and a test run may give a FASTA
holding only a sample of the master table.

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

A plasmid whose host is a eukaryote (the yeast 2-micron plasmid, plant and algal plasmids)
is excluded, since the study is of bacterial plasmids. The host is resolved as rule
clonal_registry resolves it (plasmidann.hosts; the GenBank organism is not used for a
metagenomic record), and is eukaryotic when the NCBI taxonomy (input.ncbi_taxdump) places
every taxid of its name in Eukaryota (plasmidann.hosts.eukaryotic). The number removed,
per host, is printed.
"""
import collections
import csv

import _ctx  # noqa: F401

from darkorf.circular import is_circular, terminal_repeat_length
from plasmidann import hosts
from plasmidann.fasta import iter_fasta

exclude = set(snakemake.params.exclude)
max_size = snakemake.params.max_size_bp
min_repeat = snakemake.params.min_terminal_repeat_bp
in_scope = {}                                   # plasmid_id -> (circular, small)
plsdb_species = {}
with open(snakemake.input.master, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["hab_top"] not in exclude:
            in_scope[row["plasmid_id"]] = (is_circular(row.get("topology")),
                                           int(row["size_bp"]) < max_size)
            plsdb_species[row["plasmid_id"]] = row.get("plsdb_species", "")

with open(snakemake.input.ps_hosts, newline="") as fh:
    ps_host = {r["plasmid_id"]: r.get("host", "") for r in csv.DictReader(fh, delimiter="\t")}
with open(snakemake.input.working_set, newline="") as fh:
    organism = {r["plasmid_id"]: r.get("organism", "")
                for r in csv.DictReader(fh, delimiter="\t") if r["lifestyle"] != "metagenomic"}
host = {}
for pid in in_scope:
    species, genus, _ = hosts.resolve([
        ("plsdb", plsdb_species[pid]), ("plasmidscope", ps_host.get(pid, "")),
        ("organism", organism.get(pid, ""))])
    host[pid] = species or genus
eukaryotic = hosts.eukaryotic(set(host.values()) - {""}, snakemake.params.taxdump)
removed = collections.Counter(host[pid] for pid in in_scope if host[pid] in eukaryotic)
for pid in [pid for pid in in_scope if host[pid] in eukaryotic]:
    del in_scope[pid]
print(f"analysis set: {sum(removed.values())} plasmids with a eukaryotic host excluded: "
      + (", ".join(f"{name} {n}" for name, n in removed.most_common()) or "none"))

n_seen = n_small = n_trimmed = 0
written = set()
with open(snakemake.output.fasta, "w") as out, \
        open(snakemake.output.ids, "w") as ids, \
        open(snakemake.output.small_ids, "w") as small_ids, \
        open(snakemake.output.repeats, "w") as rep, \
        open(snakemake.output.lengths, "w") as lengths:
    rep.write("plasmid_id\trecord_bp\trepeat_bp\tmolecule_bp\n")
    lengths.write("plasmid_id\tlength_bp\n")
    for pid, seq in iter_fasta([snakemake.input.fasta]):
        n_seen += 1
        if pid not in in_scope:
            continue
        if pid in written:
            raise SystemExit(f"analysis_set: {pid} occurs twice in {snakemake.input.fasta}")
        written.add(pid)
        circular, small = in_scope[pid]
        k = terminal_repeat_length(seq.upper(), min_repeat) if circular else 0
        if k:
            n_trimmed += 1
            rep.write(f"{pid}\t{len(seq)}\t{k}\t{len(seq) - k}\n")
            seq = seq[:-k]
        out.write(f">{pid}\n{seq}\n")
        lengths.write(f"{pid}\t{len(seq)}\n")
        ids.write(pid + "\n")
        if small:
            n_small += 1
            small_ids.write(pid + "\n")

# An in-scope set with no sequence is never a legitimate result: every later stage would
# succeed while writing well-formed empty tables.
if not written:
    raise SystemExit(
        f"analysis_set: none of the {n_seen} records in {snakemake.input.fasta} is among "
        f"the {len(in_scope)} in-scope plasmids. Check that the FASTA identifiers are the "
        "master table's plasmid_id.")
print(f"analysis set: {len(written)} plasmids, {n_small} of them small "
      f"(size_bp < {max_size}); {len(in_scope) - len(written)} in-scope plasmids of the "
      f"master table have no record in the FASTA; {n_seen} FASTA records read; "
      f"{n_trimmed} circular records carried a terminal repeat of >= {min_repeat} bp, "
      "one copy removed")
