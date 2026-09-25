"""Scratch directories are removed when a stage succeeds and kept when it fails.

Six scripts created a working directory beside their output and none removed it. On the
100-plasmid test set that left mmseqs_tmp_close, mmseqs_tmp_intermediate,
mmseqs_tmp_broad, foldseek_tmp and two anonymous mkdtemp directories behind; at production
scale the same directories hold the intermediate databases of a 3.5M-protein clustering,
so the leak is measured in hundreds of gigabytes rather than in tidiness.
"""
from plasmidann.scratch import release, scratch_dir


def test_a_named_directory_is_stable_and_readable(tmp_path):
    """mmseqs_tmp_broad rather than tmpaixkezqn: a directory left behind by a failure has
    to say which stage and which resolution produced it."""
    path = scratch_dir(tmp_path, "mmseqs_tmp_broad")

    assert path == tmp_path / "mmseqs_tmp_broad"
    assert path.is_dir()


def test_release_removes_the_directory_and_its_contents(tmp_path):
    path = scratch_dir(tmp_path, "work")
    (path / "intermediate.db").write_text("x")

    release(path)

    assert not path.exists(), (
        "the scratch directory survived; at production scale these hold the intermediate "
        "databases of a 3.5M-protein clustering")


def test_a_stale_directory_from_a_failed_attempt_is_cleared_first(tmp_path):
    """Failed runs keep their scratch directory on purpose - it is what a cluster failure
    is diagnosed from - so stale directories are guaranteed to exist. mmseqs refuses to
    start against a tmp directory whose contents do not match the current run, so a rerun
    would otherwise fail for a reason unrelated to the original one."""
    stale = tmp_path / "work"
    stale.mkdir()
    (stale / "leftover").write_text("from the run that crashed")

    path = scratch_dir(tmp_path, "work")

    assert not (path / "leftover").exists()


def test_an_unnamed_directory_is_unique_per_caller(tmp_path):
    """Sharded rules run concurrently in the same output directory. Two shards sharing one
    scratch path would overwrite each other's intermediates."""
    a = scratch_dir(tmp_path)
    b = scratch_dir(tmp_path)

    assert a != b
    assert a.is_dir() and b.is_dir()


def test_releasing_a_directory_that_is_already_gone_is_not_an_error(tmp_path):
    """release is the last statement of a script. It must not turn a successful run into
    a failed one because a filesystem removed the directory first."""
    release(tmp_path / "never_existed")
