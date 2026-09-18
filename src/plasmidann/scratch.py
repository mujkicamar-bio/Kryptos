"""Working directories for stages that shell out to a tool.

Six stages created a working directory beside their output - mmseqs, foldseek, hmmsearch,
DIAMOND and the alignment steps all need somewhere to put intermediate databases - and
none of them removed it. On the 100-plasmid test set that left seven directories behind.
At production scale the same directories hold the intermediate databases of a 3.5M-protein
clustering, so the leak is measured in hundreds of gigabytes rather than in tidiness.

REMOVED ON SUCCESS, KEPT ON FAILURE

The obvious repair is try/finally, and it is the wrong one. A failed mmseqs or foldseek run
is diagnosed from exactly these files: the partial database, the tool's own log, the query
set as the tool saw it. A stage that deleted them on the way out would leave a cluster
failure with nothing to look at, which is a worse problem than disk usage.

So `release` is called at the END of a script, on the ordinary path. A Snakemake script
that raises never reaches it, and the directory that survives a run is a directory with
something wrong in it.

TWO FUNCTIONS RATHER THAN A CONTEXT MANAGER

These scripts use their scratch path at module level, from creation near the top to the
last subprocess call at the bottom. A `with` block would mean indenting the body of six
scripts to gain the one behaviour - deletion on failure - that is deliberately not wanted.
"""
import pathlib
import shutil
import tempfile


def scratch_dir(parent, name=None):
    """Create and return a working directory under `parent`.

    `name` gives a stable, readable path - mmseqs_tmp_broad rather than tmpaixkezqn - for a
    stage that produces one output. Omit it for a stage that runs concurrently across
    shards: two shards sharing one scratch path would overwrite each other's intermediates,
    and the corruption would present as a search error rather than as a collision.

    A stale directory from a previous failed attempt is cleared before use. Keeping failed
    directories guarantees those exist, and mmseqs refuses to start against a tmp directory
    whose contents do not match the current run - so the rerun would fail for a reason that
    has nothing to do with the original one.
    """
    parent = pathlib.Path(parent)
    parent.mkdir(parents=True, exist_ok=True)
    if name is None:
        return pathlib.Path(tempfile.mkdtemp(dir=parent))
    path = parent / name
    if path.exists():
        shutil.rmtree(path)
    path.mkdir()
    return path


def release(path):
    """Remove a scratch directory. Called at the end of a script, on success only."""
    shutil.rmtree(path, ignore_errors=True)
