"""Discovering the input shards the pipeline will process.

WHY THE PIPELINE NO LONGER CREATES ITS SHARDS

It used to: a rule read the whole working-set FASTA and dealt records round-robin into a
fixed number of files, and the count lived in configuration. Two consequences made that the
wrong model for this project.

  * The count was load-bearing. Changing it invalidated every shard and therefore every
    downstream stage, so the partition could not be adjusted without discarding work.
  * Nothing could be added. A new batch of plasmids meant re-dealing the whole collection,
    because a round-robin partition has no stable membership.

Shards are now prepared OUTSIDE the pipeline and discovered from a directory. Adding a
batch means dropping one more file in and re-running: the per-shard stages run for the new
shard alone, and only the stages that aggregate across shards re-run. That is the property
the round-robin partition could not give.

WHAT COUNTS AS A SHARD

Any FASTA in the configured directory, plain or gzipped. The shard NAME is the file name
with its extension removed, and it becomes a wildcard value, so it must contain no path
separator and no whitespace.

WHY AN EMPTY DIRECTORY IS AN ERROR

With no shards, every aggregating rule has an empty input list, every stage succeeds
trivially, and the run produces well-formed empty tables. This pipeline has already been
burned by exactly that class of failure - four stages once reported success while writing
nothing - so an empty shard set is refused at workflow load time, before any job is
scheduled.
"""
import pathlib
import re

# Extensions accepted for a shard, longest first so '.fna.gz' is stripped whole rather than
# leaving a stranded '.fna' behind the shorter match.
SHARD_SUFFIXES = (".fna.gz", ".fasta.gz", ".fa.gz", ".fna", ".fasta", ".fa")

# A shard name becomes a Snakemake wildcard and a path component. Anything outside this set
# either breaks the DAG or silently changes where an output lands.
_VALID_NAME = re.compile(r"^[A-Za-z0-9._-]+$")


def shard_name(path):
    """The shard name for a file: its basename with the FASTA extension removed."""
    name = pathlib.Path(path).name
    for suffix in SHARD_SUFFIXES:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


def discover_shards(directory):
    """Every shard in `directory`, as {name: path}, sorted by name.

    Sorted so that the DAG is identical between runs on the same inputs: an unsorted
    directory listing is filesystem order, which would change the job order and every log
    that records it for no biological reason.

    Raises when the directory is missing or holds no FASTA. An empty shard set would give
    every aggregating rule an empty input, and the run would report success having written
    well-formed empty tables - the failure this pipeline has already had four times.
    """
    directory = pathlib.Path(directory)
    if not directory.is_dir():
        raise SystemExit(
            f"shard directory {directory} does not exist. Shards are prepared outside the "
            "pipeline and discovered from this directory; create it and add at least one "
            "FASTA before running.")

    found = {}
    for path in sorted(directory.iterdir()):
        if not path.is_file() or not path.name.endswith(SHARD_SUFFIXES):
            continue
        name = shard_name(path)
        if not _VALID_NAME.match(name):
            raise SystemExit(
                f"shard name {name!r} (from {path.name}) contains characters that cannot "
                "be a Snakemake wildcard. Use letters, digits, dot, dash and underscore.")
        if name in found:
            raise SystemExit(
                f"two files in {directory} both give the shard name {name!r}: "
                f"{found[name].name} and {path.name}. A shard name must identify one file, "
                "or the two would overwrite each other's outputs.")
        found[name] = path

    if not found:
        raise SystemExit(
            f"no FASTA shards found in {directory} (looked for {', '.join(SHARD_SUFFIXES)}). "
            "An empty shard set would let every stage succeed while writing nothing.")
    return {name: str(found[name]) for name in sorted(found)}
