"""Working directories for stages that shell out to a tool.

A script creates its directory with scratch_dir and calls release as its last statement,
so the directory is removed on success and kept on failure: a failed mmseqs or foldseek
run is diagnosed from exactly these files (the partial database, the tool's log, the query
as the tool saw it). Two functions rather than a context manager: what a context manager
would add is deletion on failure, which is the one behaviour not wanted.
"""
import pathlib
import shutil
import tempfile


def scratch_dir(parent, name=None):
    """Create and return a working directory under `parent`.

    `name` gives a stable, readable path - mmseqs_tmp_broad rather than tmpaixkezqn. Without
    it the directory gets a unique name, for a stage whose jobs may share an output
    directory and must not overwrite each other's intermediates.

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
