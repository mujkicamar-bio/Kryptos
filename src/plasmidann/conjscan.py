"""Conjugation and mobilisation systems from CONJScan, and the plasmid mobility class.

CONJScan (Cury et al. 2020, Methods Mol. Biol. 2075:265) is a set of MacSyFinder models.
Its `Plasmids` set is the one built for plasmids (Coluzzi et al. 2022, Mol. Biol. Evol.
39:msac115): T4SS_type<X> is a complete mating-pair formation (MPF) machinery, dCONJ_type<X>
a decayed one, and MOB a relaxase without an MPF. Thresholds are the shipped model
definitions, not overridden, for the same reason DefenseFinder's are not: the definitions
ARE the published models.

THE CLASS

Coluzzi et al. 2022 ("Assessment of Plasmid Mobility") call a plasmid conjugative (pCONJ)
when it encodes a relaxase, a VirB4, a T4CP and a minimum number of further MPF proteins;
decayed conjugative (pdCONJ) when it encodes a relaxase and an incomplete MPF; mobilisable
(pMOB) when it encodes a relaxase with no or very few MPF genes; and pMOBless when it
encodes no relaxase. In the CONJScan 2.1.0 Plasmids definitions the relaxase is mandatory
in every model, and the T4SS, dCONJ and MOB models implement the three relaxase-carrying
classes with their own quorums (dCONJ_typeF, for instance, needs five genes including the
relaxase); those quorums are the models', not re-derived here. A plasmid takes the most
complete system it carries, and a plasmid with no system carries no detected relaxase:
pMOBless.
"""
import pathlib

MODEL_SET = "CONJScan/Plasmids"

COLUMNS = ["orf_id", "plasmid_id", "system", "system_id", "component", "hit_status",
           "sys_wholeness", "conjscan_version", "status"]
CLASS_COLUMNS = ["plasmid_id", "class"]

# Model-name prefix -> class, most complete first.
_CLASS_OF_PREFIX = (("T4SS_type", "pCONJ"), ("dCONJ_type", "pdCONJ"), ("MOB", "pMOB"))


def installed_version(models_dir):
    """The CONJScan version under `models_dir`, from the package's metadata.yml `vers:`.

    None when the package is not installed there. The metadata is read directly rather
    than through macsydata, so the check needs no MacSyFinder environment.
    """
    meta = pathlib.Path(models_dir) / "CONJScan" / "metadata.yml"
    if not meta.is_file():
        return None
    for line in meta.read_text().splitlines():
        if line.startswith("vers:"):
            return line.split(":", 1)[1].strip().strip("'\"")
    return None


def read_best_solution(path):
    """Rows of one MacSyFinder best_solution.tsv, as dicts keyed by its header.

    The file opens with '#' comment lines and a blank line before the header, and a run
    that found no system is a file of comments alone.
    """
    with open(path, newline="") as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    if not lines:
        return []
    header = lines[0].rstrip("\n").split("\t")
    return [dict(zip(header, l.rstrip("\n").split("\t"))) for l in lines[1:]]


def plasmid_class(system_types):
    """pCONJ / pdCONJ / pMOB / pMOBless from the CONJScan system types on one plasmid.

    `system_types` are model names (`T4SS_typeF`, `dCONJ_typeG`, `MOB`). An unknown type
    raises: a model renamed in a later release must stop the stage, not fall through to
    pMOBless.
    """
    found = set()
    for t in system_types:
        cls = next((c for p, c in _CLASS_OF_PREFIX if t.startswith(p)), None)
        if cls is None:
            raise ValueError(f"unknown CONJScan system type {t!r}")
        found.add(cls)
    for _, cls in _CLASS_OF_PREFIX:
        if cls in found:
            return cls
    return "pMOBless"
