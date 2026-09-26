"""tools/calibrate_context.py: how often a context term predicts a family's own function.

Benchmark per term: known families whose own members carry the term (majority); negatives:
labelled known families no member of which carries it. Precision at each conservation
level, and the lowest level reaching 50% and 90% precision.
"""
import importlib.util
import pathlib

import pytest
from conftest import read_tsv, write_tsv

TOOL = pathlib.Path(__file__).resolve().parents[1] / "tools" / "calibrate_context.py"


@pytest.fixture(scope="module")
def cal():
    spec = importlib.util.spec_from_file_location("calibrate_context", TOOL)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --- precision curve and thresholds -------------------------------------------------------

def test_precision_is_counted_at_and_above_each_conservation_level(cal):
    scores = {"p1": 0.9, "p2": 0.9, "n1": 0.9, "p3": 0.3, "n2": 0.3, "n3": 0.0}
    curve = cal.precision_curve(scores, positives={"p1", "p2", "p3"})
    assert [(c["conservation"], c["n_predicted"], c["n_true"]) for c in curve] == [
        (0.3, 5, 3), (0.9, 3, 2)]
    assert curve[1]["precision"] == pytest.approx(2 / 3)


def test_the_threshold_is_the_lowest_level_whose_precision_holds_above_it(cal):
    curve = [{"conservation": 0.2, "n_predicted": 40, "n_true": 16, "precision": 0.4},
             {"conservation": 0.5, "n_predicted": 30, "n_true": 18, "precision": 0.6},
             {"conservation": 0.7, "n_predicted": 20, "n_true": 19, "precision": 0.95},
             {"conservation": 0.8, "n_predicted": 12, "n_true": 11, "precision": 0.917}]
    assert cal.threshold(curve, 0.5, min_families=10)["conservation"] == 0.5
    assert cal.threshold(curve, 0.9, min_families=10)["conservation"] == 0.7
    assert cal.threshold(curve, 0.99, min_families=10) is None


def test_levels_with_too_few_families_neither_set_nor_block_a_threshold(cal):
    """One family at 1.0 that is a false positive is noise, not evidence against 0.7; and
    a perfect level holding three families is not evidence for itself."""
    curve = [{"conservation": 0.7, "n_predicted": 20, "n_true": 19, "precision": 0.95},
             {"conservation": 0.9, "n_predicted": 3, "n_true": 3, "precision": 1.0},
             {"conservation": 1.0, "n_predicted": 1, "n_true": 0, "precision": 0.0}]
    assert cal.threshold(curve, 0.9, min_families=10)["conservation"] == 0.7
    assert cal.threshold(curve[1:], 0.9, min_families=10) is None


# --- end to end on synthetic families with known answers ------------------------------------

TERM_COLS = ["family_id", "family_set", "term_type", "term", "n_lineages",
             "n_lineages_with_term", "conservation", "status",
             "window_covers_plasmid_fraction"]


def _write_inputs(d):
    """12 amr families (P) and 12 ta families (N). For amr:x, 10 P and 1 N sit at 0.9, the
    other 2 P and 11 N at 0.3: precision 10/11 at 0.9 and 12/24 at 0.3. Noise that must not
    count: a TOO_FEW_LINEAGES positive, a minority-labelled family, a dark family, and a
    term (mge:z) with 3 benchmark families."""
    terms, labels, families, pmap = [], [], [], []
    for i in range(12):
        p, n = f"intermediate:P{i}", f"intermediate:N{i}"
        families += [[p, "intermediate", f"P{i}", f"P{i}a,P{i}b"],
                     [n, "intermediate", f"N{i}", f"N{i}a"]]
        labels += [[f"P{i}a", "card", "card_amr_family", "x", ""],
                   [f"P{i}b", "card", "card_amr_family", "x", ""],
                   [f"N{i}a", "tadb", "tadb_ta", "y", ""]]
        terms += [[p, "known", "amr", "amr:x", 3, 3 if i < 10 else 1,
                   0.9 if i < 10 else 0.3, "SUCCESS", 0.0],
                  [n, "known", "amr", "amr:x", 10, 9 if i < 1 else 3,
                   0.9 if i < 1 else 0.3, "SUCCESS", 0.0]]
    families += [["intermediate:Q", "intermediate", "Qa", "Qa"],
                 ["intermediate:M", "intermediate", "Ma", "Ma,Mb,Mc"]]
    labels += [["Qa", "card", "card_amr_family", "x", ""],
               ["Ma", "card", "card_amr_family", "x", ""]]
    terms += [["intermediate:Q", "known", "amr", "amr:x", 1, 1, 1.0, "TOO_FEW_LINEAGES", 0],
              ["intermediate:M", "known", "amr", "amr:x", 2, 2, 1.0, "SUCCESS", 0],
              ["intermediate:D", "dark", "amr", "amr:x", 5, 5, 1.0, "SUCCESS", 0]]
    for i in range(3):
        fid = f"intermediate:Z{i}"
        families.append([fid, "intermediate", f"Z{i}", f"Z{i}a"])
        labels.append([f"Z{i}a", "mobileog", "mobileog_category", "z", ""])
        terms.append([fid, "known", "mge", "mge:z", 4, 4, 1.0, "SUCCESS", 0])
    # A family whose member is a DefenseFinder component carries defence:Clover itself.
    families.append(["intermediate:S", "intermediate", "Sa", "Sa"])
    pmap.append("Sa\tpl9|4")
    terms.append(["intermediate:S", "known", "defence", "defence:Clover", 2, 2, 1.0,
                  "SUCCESS", 0])

    write_tsv(d / "terms.tsv", TERM_COLS, terms)
    write_tsv(d / "labels.tsv", ["protein_id", "source", "kind", "label", "sub_label"],
              labels)
    write_tsv(d / "families.tsv", ["family_id", "family_resolution", "representative",
                                   "members"], families)
    (d / "map.tsv").write_text("\n".join(pmap) + "\n")
    write_tsv(d / "defence.tsv", ["orf_id", "plasmid_id", "system"],
              [["pl9|4", "pl9", "defense-finder-models/DefenseFinder/Clover/Clover"]])
    write_tsv(d / "conj.tsv", ["orf_id", "plasmid_id", "system"], [])


def _run(cal, d):
    _write_inputs(d)
    cal.main(["--terms", str(d / "terms.tsv"), "--labels", str(d / "labels.tsv"),
              "--families", str(d / "families.tsv"), "--map", str(d / "map.tsv"),
              "--defence", str(d / "defence.tsv"), "--conjugation", str(d / "conj.tsv"),
              "--out", str(d / "calibration.tsv"), "--curve", str(d / "curve.tsv"),
              "--summary", str(d / "summary.txt")])
    return {r["term"]: r for r in read_tsv(d / "calibration.tsv")}


def test_thresholds_are_found_at_the_right_levels(cal, tmp_path):
    amr = _run(cal, tmp_path)["amr:x"]
    assert amr["status"] == "CALIBRATED"
    # Negatives: the 12 ta families, the 3 mge families and the defence family - every
    # evaluable labelled family without amr:x, whether or not it has an amr:x context row.
    # Q (too few lineages), M (minority) and the dark family D are neither.
    assert (amr["n_benchmark"], amr["n_negative"]) == ("12", "16")
    assert amr["threshold_50"] == "0.3" and amr["n_at_50"] == "24"
    assert amr["threshold_90"] == "0.9" and amr["n_at_90"] == "11"
    assert float(amr["precision_90"]) == pytest.approx(10 / 11, abs=1e-4)


def test_too_small_a_benchmark_is_uncalibrated(cal, tmp_path):
    rows = _run(cal, tmp_path)
    assert rows["mge:z"]["status"] == "UNCALIBRATED"
    assert rows["mge:z"]["n_benchmark"] == "3"
    assert rows["mge:z"]["threshold_50"] == "" and rows["mge:z"]["threshold_90"] == ""
    # A system term is carried by a family whose member is itself a component.
    assert rows["defence:Clover"]["n_benchmark"] == "1"
    assert rows["defence:Clover"]["status"] == "UNCALIBRATED"


def test_the_curve_and_summary_are_written(cal, tmp_path):
    _run(cal, tmp_path)
    curve = [r for r in read_tsv(tmp_path / "curve.tsv") if r["term"] == "amr:x"]
    assert [(r["conservation"], r["n_predicted"], r["n_true"]) for r in curve] == [
        ("0.3", "24", "12"), ("0.9", "11", "10")]
    summary = (tmp_path / "summary.txt").read_text()
    assert "amr:x" in summary and "UNCALIBRATED" in summary
