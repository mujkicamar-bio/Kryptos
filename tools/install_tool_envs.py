"""Build the CONJScan, AMRFinderPlus and geNomad environments and install their reference
data.

WHAT IS INSTALLED

  envs/conjscan         from workflow/envs/conjscan.yaml (MacSyFinder 2.1.6, HMMER 3.4)
  data/refs/conjscan    CONJScan 2.1.0 models, `msf_data install` from that environment
  envs/amrfinder        from workflow/envs/amrfinder.yaml (NCBI AMRFinderPlus)
  data/refs/amrfinder   the AMRFinderPlus database, `amrfinder_update` from that environment
  envs/genomad          from workflow/envs/genomad.yaml (geNomad)
  data/refs/genomad     the geNomad database, `genomad download-database` from that
                        environment

Each tool runs from its own environment, named by path, as pharokka does, because its pins
conflict with the main environment (see the comments in the YAML files).

IDEMPOTENT

An environment is rebuilt only when a pinned package is missing or at another version, the
models only when metadata.yml does not report 2.1.0, the AMRFinderPlus database only when
no AMRProt file is present, and the geNomad database only when genomad_db/version.txt is
absent. A second run therefore installs nothing.

The package and repodata caches are a temporary directory inside envs/, removed
afterwards, so nothing is written to the home directory (which has a file quota) and
packages are hard-linked rather than copied.

Run from the repository root:  python tools/install_tool_envs.py
"""
import argparse
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

# Must equal conjugation.version in config/targets.yaml (tests/test_config.py checks it);
# this script uses the standard library only, so it cannot read the YAML.
CONJSCAN_VERSION = "2.1.0"


def pins(spec):
    """The `name=version` pins of an environment YAML, as {name: version}."""
    return dict(re.findall(r"^\s*-\s*([A-Za-z0-9_.-]+)=([^\s#]+)", spec.read_text(), re.M))


def installed(prefix):
    """{name: version} of the packages in a conda prefix, empty when there is none."""
    meta = prefix / "conda-meta"
    if not meta.is_dir():
        return {}
    out = {}
    for p in meta.glob("*.json"):
        d = json.loads(p.read_text())
        out[d["name"]] = d["version"]
    return out


def matches(have, want):
    """A conda pin `x=3.12` matches 3.12 and 3.12.4, not 3.1 or 3.120."""
    return have is not None and (have == want or have.startswith(want + "."))


def solver():
    for name in ("micromamba", "mamba", "conda"):
        path = shutil.which(name)
        if path:
            return name, path
    sys.exit("none of micromamba, mamba or conda is on PATH")


def build_env(root, name):
    spec = root / "workflow/envs" / f"{name}.yaml"
    prefix = root / "envs" / name
    want, have = pins(spec), installed(prefix)
    stale = {k: v for k, v in want.items() if not matches(have.get(k), v)}
    if not stale:
        print(f"envs/{name}: present, all {len(want)} pins match - skipped")
        return
    print(f"envs/{name}: building ({', '.join(f'{k}={v}' for k, v in stale.items())} "
          f"{'missing or different' if have else 'absent'})")
    if prefix.exists():
        shutil.rmtree(prefix)
    kind, exe = solver()
    with tempfile.TemporaryDirectory(dir=root / "envs", prefix=".pkgs-") as cache:
        env = dict(os.environ, CONDA_PKGS_DIRS=cache, MAMBA_ROOT_PREFIX=cache,
                   XDG_CACHE_HOME=cache)  # micromamba keeps its repodata shards there
        if kind == "conda":
            cmd = [exe, "env", "create", "-y", "-p", str(prefix), "-f", str(spec)]
        else:
            # --no-rc: the site configuration redirects channels to a local mirror; the
            # environment file names conda-forge and bioconda, and those are what is solved.
            cmd = [exe, "create", "-y", "--no-rc", "--override-channels",
                   "-c", "conda-forge", "-c", "bioconda", "--strict-channel-priority",
                   "-p", str(prefix), "-f", str(spec)]
        subprocess.run(cmd, env=env, check=True)
    have = installed(prefix)
    bad = {k: have.get(k) for k, v in want.items() if not matches(have.get(k), v)}
    if bad:
        sys.exit(f"envs/{name}: pins not satisfied after install: {bad}")


def install_conjscan(root):
    target = root / "data/refs/conjscan"
    meta = target / "CONJScan/metadata.yml"
    if meta.is_file() and re.search(rf"^vers:\s*{re.escape(CONJSCAN_VERSION)}\s*$",
                                    meta.read_text(), re.M):
        print(f"data/refs/conjscan: CONJScan {CONJSCAN_VERSION} present - skipped")
        return
    target.mkdir(parents=True, exist_ok=True)
    # msf_data is MacSyFinder 2.1.6's name for macsydata; it needs hmmsearch on PATH.
    bindir = root / "envs/conjscan/bin"
    env = dict(os.environ, PATH=f"{bindir}{os.pathsep}{os.environ['PATH']}")
    subprocess.run([str(bindir / "msf_data"), "install", "--force", "--target", str(target),
                    f"CONJScan=={CONJSCAN_VERSION}"], env=env, check=True)
    if not meta.is_file() or CONJSCAN_VERSION not in meta.read_text():
        sys.exit(f"CONJScan {CONJSCAN_VERSION} not found in {meta} after install")


def install_amrfinder_db(root):
    target = root / "data/refs/amrfinder"
    if any(target.glob("**/AMRProt*")):
        print("data/refs/amrfinder: database present - skipped")
    else:
        target.mkdir(parents=True, exist_ok=True)
        subprocess.run([str(root / "envs/amrfinder/bin/amrfinder_update"),
                        "--database", str(target)], check=True)
    # The config names data/refs/amrfinder/latest; amrfinder_update creates that link, and
    # it is recreated here if missing, pointing at the newest installed version directory.
    latest = target / "latest"
    if not latest.exists():
        dirs = sorted(p.parent.name for p in target.glob("*/version.txt")
                      if not p.parent.is_symlink())
        if not dirs:
            sys.exit(f"no installed database version under {target}")
        latest.unlink(missing_ok=True)
        latest.symlink_to(dirs[-1])
    version = latest / "version.txt"
    if not version.is_file():
        sys.exit(f"no latest/version.txt under {target}")
    print(f"data/refs/amrfinder: database version {version.read_text().strip()}")


def install_genomad_db(root):
    target = root / "data/refs/genomad"
    version = target / "genomad_db/version.txt"
    if not version.is_file():
        target.mkdir(parents=True, exist_ok=True)
        subprocess.run([str(root / "envs/genomad/bin/genomad"), "download-database",
                        str(target)], check=True)
    if not version.is_file():
        sys.exit(f"no {version} after download")
    print(f"data/refs/genomad: database version {version.read_text().strip()}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".", help="repository root")
    root = pathlib.Path(ap.parse_args().root).resolve()
    build_env(root, "conjscan")
    install_conjscan(root)
    build_env(root, "amrfinder")
    install_amrfinder_db(root)
    build_env(root, "genomad")
    install_genomad_db(root)


if __name__ == "__main__":
    main()
