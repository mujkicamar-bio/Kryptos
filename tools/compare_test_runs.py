"""Compare the test run before the PlasmidScope change (results_test) with the one after
(results_test_ps). Run from the repository root; prints to stdout.
"""
import collections
import csv


OLD, NEW = "results_test", "results_test_ps"


def tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ids(path):
    return [l.strip() for l in open(path) if l.strip()]


def fasta_n(path):
    return sum(1 for l in open(path) if l.startswith(">"))


for run in (OLD, NEW):
    print(f"== {run}")
    print("  cascade_input sequences", fasta_n(f"{run}/03_dereplication/cascade_input.faa"))
    print("  quality gate:", open(f"{run}/09_quality_gate/quality_gate.txt").readline().strip())
    ann = tsv(f"{run}/06_annotation_tables/plasmid_annotation.tsv")
    print("  plasmid_annotation rows", len(ann),
          "| empty functional_class", sum(not r["functional_class"] for r in ann))
    prot = tsv(f"{run}/05_annotation_cascade/protein_annotation.tsv")
    real = [r for r in prot if not r["seq_id"].startswith(("CTRL_", "DECOY_"))]
    print("  protein classes", dict(collections.Counter(r["functional_class"] for r in real)))
    print("  tiers", dict(collections.Counter(r["annot_tier"] for r in real)))

import glob


def dark(run):
    p = glob.glob(f"{run}/**/dark_ids.txt", recursive=True)
    return set(ids(p[0])), p[0]

(od, op), (nd, np_) = dark(OLD), dark(NEW)
ps = {r["seq_id"]: r["ps_class"] for r in tsv(f"{NEW}/03_dereplication/plasmidscope_proteins.tsv")}
print(f"\n== dark set  old {len(od)} ({op})  new {len(nd)}")
lost, gained = od - nd, nd - od
print("  lost", len(lost), dict(collections.Counter(ps.get(s, "NOT_IN_PS") for s in lost)))
print("  gained", len(gained), dict(collections.Counter(ps.get(s, "NOT_IN_PS") for s in gained)))
print("  kept", len(od & nd))

orth = tsv(f"{NEW}/07_orthology/orthology.tsv")
print("\n== orthology source", dict(collections.Counter(r.get("orthology_source", "") for r in orth)))
