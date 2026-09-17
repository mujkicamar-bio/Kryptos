"""Shared preamble for every workflow script: import path, and log capture.

Two lines of setup that every script needs and none should repeat.

IMPORT PATH
    Snakemake executes a script with its own directory on sys.path, not the project's, so
    `from plasmidann...` fails without this.

LOG CAPTURE
    Every rule declares `log:`, and for `script:` rules Snakemake creates that path but
    does NOT redirect the script's stdout into it - that redirection only happens for
    `shell:` and `run:` directives. All 26 log files therefore stayed empty while the
    diagnostics the scripts print (hit counts, rejected hits, background rates, per-tier
    control resolution) went to one SLURM output file, interleaved across up to 1,481
    concurrent jobs and impossible to attribute to a rule.

    Output is TEED rather than moved: the log gets the rule's own record, and the SLURM
    file keeps a live progress trace. A script that declares no log is unaffected.
"""
import atexit
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "src"))


class _Tee:
    """Write to the terminal and to the rule's log file at once."""

    def __init__(self, *streams):
        self._streams = streams

    def write(self, text):
        for s in self._streams:
            s.write(text)
        return len(text)

    def flush(self):
        for s in self._streams:
            s.flush()

    def isatty(self):
        return False


def _find_snakemake():
    """The `snakemake` object lives in the globals of whichever script imported us.

    Snakemake injects it there rather than into a module, so it cannot be imported. Walking
    back up the frame stack from this import finds the importing script's globals; the
    import machinery adds a handful of frames in between, so the walk is bounded rather
    than open-ended.
    """
    frame = sys._getframe()
    for _ in range(30):
        frame = frame.f_back
        if frame is None:
            return None
        found = frame.f_globals.get("snakemake")
        if found is not None:
            return found
    return None


def _install_log_capture():
    logs = getattr(_find_snakemake(), "log", None)
    if not logs:
        return
    path = pathlib.Path(str(logs[0]))
    path.parent.mkdir(parents=True, exist_ok=True)
    # Append, so a Snakemake retry of the same job keeps the failed attempt's output next
    # to the successful one rather than erasing the evidence.
    handle = open(path, "a", buffering=1)
    sys.stdout = _Tee(sys.__stdout__, handle)
    sys.stderr = _Tee(sys.__stderr__, handle)

    @atexit.register
    def _restore():
        sys.stdout, sys.stderr = sys.__stdout__, sys.__stderr__
        handle.close()


_install_log_capture()
