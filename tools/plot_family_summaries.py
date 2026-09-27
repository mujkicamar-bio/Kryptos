"""Summary figures for the protein families of one pipeline run.

Every figure is drawn from the run's own tables, so it regenerates exactly:
  family_sizes         family size (members) per clustering resolution, log-log counts,
                       over the families holding a small-plasmid protein
  dark_fraction        percentage of dark members per family, families with >= 2 members
  dark_family_classes  dark-only vs mixed dark families, and ORPHAN vs FAMILY
  rarefaction          dark families against lineages sampled (mean and range)
  network              degree distribution and community sizes of the family network
  dark_orf_lengths     protein length (aa) of dark ORFs, FUNCTIONAL ORFs for comparison

Usage:
  python tools/plot_family_summaries.py --run results_test_ps [--out DIR]
Figures go to <run>/15_report/figures/ as PDF and SVG.
"""
import argparse
import collections
import csv
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

# Okabe-Ito, colour-blind safe.
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def save(fig, out, name):
    fig.tight_layout()
    for ext in ("pdf", "svg"):
        fig.savefig(out / f"{name}.{ext}")
    plt.close(fig)
    print(out / f"{name}.pdf")


def family_sizes(run, out):
    by_res = collections.defaultdict(collections.Counter)
    for r in read_tsv(run / "10_clustering/protein_families.tsv"):
        by_res[r["family_resolution"]][int(r["family_size"])] += 1
    fig, ax = plt.subplots(figsize=(5, 4))
    for colour, res in zip(PALETTE, sorted(by_res)):
        sizes = sorted(by_res[res])
        ax.plot(sizes, [by_res[res][s] for s in sizes], "o-", ms=3, color=colour,
                label=f"{res} ({sum(by_res[res].values()):,} families)")
    ax.set_title("families with a small-plasmid protein", fontsize=9)
    ax.set(xscale="log", yscale="log", xlabel="family size (unique proteins)",
           ylabel="number of families")
    ax.legend(frameon=False, fontsize=8)
    save(fig, out, "family_sizes")


def dark_fraction(run, out):
    rows = [r for r in read_tsv(run / "10_clustering/dark_families.tsv")
            if int(r["n_members"]) >= 2]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.hist([float(r["percentage_dark_in_family"]) for r in rows], bins=20, range=(0, 100),
            color=PALETTE[0])
    ax.set(xlabel="dark members in family (%)", ylabel="dark families (>= 2 members)",
           title=f"{len(rows):,} families with >= 2 members")
    save(fig, out, "dark_fraction")


def dark_family_classes(run, out):
    rows = read_tsv(run / "10_clustering/dark_families.tsv")
    counts = collections.Counter(
        (r["family_class"], "dark-only" if r["dark_only"] == "1" else "mixed") for r in rows)
    labels = sorted(counts)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.barh([f"{a}, {b}" for a, b in labels], [counts[k] for k in labels],
            color=PALETTE[:len(labels)])
    for i, k in enumerate(labels):
        ax.text(counts[k], i, f" {counts[k]:,}", va="center", fontsize=8)
    ax.set(xlabel="dark families", title=f"{len(rows):,} dark families")
    save(fig, out, "dark_family_classes")


def rarefaction(run, out):
    rows = read_tsv(run / "15_report/dark_family_rarefaction.tsv")
    x = [int(r["n_lineages"]) for r in rows]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.fill_between(x, [float(r["min_families"]) for r in rows],
                    [float(r["max_families"]) for r in rows], color=PALETTE[0], alpha=0.25,
                    lw=0, label="range over replicates")
    ax.plot(x, [float(r["mean_families"]) for r in rows], color=PALETTE[0], label="mean")
    ax.set(xlabel="lineages sampled", ylabel="dark families observed")
    ax.legend(frameon=False, fontsize=8)
    save(fig, out, "rarefaction")


def network(run, out):
    nodes = read_tsv(run / "10_clustering/network_nodes.tsv")
    degree = collections.Counter(int(r["degree"]) for r in nodes)
    comm = collections.Counter(collections.Counter(r["community_id"] for r in nodes).values())
    fig, (a, b) = plt.subplots(1, 2, figsize=(9, 3.8))
    d = sorted(degree)
    a.bar(d, [degree[k] for k in d], color=PALETTE[0])
    a.set(yscale="log", xlabel="degree", ylabel="nodes (50%-identity clusters)")
    a.xaxis.set_major_locator(MaxNLocator(integer=True))
    c = sorted(comm)
    b.plot(c, [comm[k] for k in c], "o", ms=3, color=PALETTE[3])
    b.set(xscale="log", yscale="log", xlabel="community size (nodes)",
          ylabel="communities")
    save(fig, out, "network")


def dark_orf_lengths(run, out):
    """Every ORF of a dark-set protein (dark_ids.txt: no artefacts, no partials), against
    the ORFs classed FUNCTIONAL. Counted per ORF, so a protein on many plasmids counts
    once per copy."""
    dark = {line.strip() for line in open(run / "10_clustering/dark_ids.txt") if line.strip()}
    dark_orfs = set()
    with open(run / "03_dereplication/protein_map.tsv") as fh:
        for line in fh:
            sid, orf_ids = line.rstrip("\n").split("\t")
            if sid in dark:
                dark_orfs.update(orf_ids.split(","))
    # Both tables have a row per ORF, with sequences in orf_index, so they are streamed and
    # only the ORF ids and lengths needed are kept.
    functional = set()
    with open(run / "06_annotation_tables/plasmid_annotation.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["functional_class"] == "FUNCTIONAL":
                functional.add(r["orf_id"])
    dark_lengths, functional_lengths = [], []
    with open(run / "02_orf_calling/orf_index.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["orf_id"] in dark_orfs:
                dark_lengths.append(len(r["seq"]))
            elif r["orf_id"] in functional:
                functional_lengths.append(len(r["seq"]))
    groups = [("dark", sorted(dark_lengths), PALETTE[3]),
              ("FUNCTIONAL", sorted(functional_lengths), PALETTE[0])]
    top = max(max(v) for _, v, _ in groups)
    edges = [20 * (top / 20) ** (i / 40) for i in range(41)]   # 40 log-spaced bins
    fig, ax = plt.subplots(figsize=(5.5, 4))
    for name, values, colour in groups:
        median = values[len(values) // 2]
        ax.hist(values, bins=edges, histtype="step" if name == "FUNCTIONAL" else "stepfilled",
                color=colour, alpha=1 if name == "FUNCTIONAL" else 0.6, lw=1.5,
                label=f"{name} ORFs (n = {len(values):,}, median {median} aa)")
        ax.axvline(median, color=colour, ls="--", lw=1, ymax=0.75)
    ax.set(xscale="log", xlabel="protein length (aa)", ylabel="ORFs")
    ax.set_ylim(top=ax.get_ylim()[1] * 1.3)    # headroom so the legend clears the bars
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    save(fig, out, "dark_orf_lengths")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True, type=pathlib.Path)
    ap.add_argument("--out", type=pathlib.Path)
    a = ap.parse_args()
    out = a.out or a.run / "15_report/figures"
    out.mkdir(parents=True, exist_ok=True)
    for figure in (family_sizes, dark_fraction, dark_family_classes, rarefaction, network,
                   dark_orf_lengths):
        figure(a.run, out)


if __name__ == "__main__":
    main()
