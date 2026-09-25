"""Gene-neighbourhood diagrams for dark families.

For each requested family, every occurrence of a member (dark or annotated) is drawn as
one row: the family's gene in the centre, turned so that it points right, with the genes around it drawn as arrows.
Rows are placed one above the other so that a conserved neighbourhood shows as a column
pattern.

Colours use only evidence that needs no text matching:
    the family's own gene      the focal dark gene
    dark                       NONE or UNCHARACTERIZED_HOMOLOG
    annotated                  any other class; the arrow carries its annot_label
    defence / integron / IS    the gene lies in a DefenseFinder system, an IntegronFinder
                               array or an ISEScan element (drawn as a band behind it)
Functional categories (replication, relaxase ...) are not coloured: the pipeline does not
group labels into categories.

Which families to draw is the reader's choice - this is a tool, not a pipeline stage.
Occurrences shown: one per plasmid, plasmids in id order, at most --max-rows; the figure
title states how many of how many are shown.

Usage:
  python tools/draw_neighbourhoods.py --run results_test_ps --family intermediate:<seq_id> \
      [--family ...] [--families-file ids.txt] [--window 3] [--max-rows 30] [--out DIR] \
      [--master data/plasmidscope_primary/analysis_set.tsv]
"""
import argparse
import collections
import csv
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrow, Rectangle  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from plasmidann.context import overlapping_islands  # noqa: E402
from plasmidann.labels import parse_ncbi_title  # noqa: E402

COLOURS = {"focal": "#D55E00", "dark": "#444444", "annotated": "#88CCEE",
           "not_searched": "#DDDDDD",
           "defence": "#CC79A7", "integron": "#009E73", "is_element": "#E69F00"}
DARK_CLASSES = {"NONE", "UNCHARACTERIZED_HOMOLOG"}


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load(run, families, master):
    """Everything needed for the requested families, reading the large tables once."""
    dark_families = {r["family_id"] for r in read_tsv(run / "10_clustering/dark_families.tsv")}
    missing = set(families) - dark_families
    if missing:
        sys.exit(f"not dark families in {run}: {sorted(missing)}")
    # dark_families.tsv lists only the dark members; every member, annotated homologues
    # included, is drawn, so the rows match the family's n_plasmids.
    fam_members = {r["family_id"]: r["members"].split(",")
                   for r in read_tsv(run / "10_clustering/protein_families.tsv")
                   if r["family_id"] in families}
    wanted = {m for ms in fam_members.values() for m in ms}

    orfs_of = {}
    with open(run / "03_dereplication/protein_map.tsv") as fh:
        for line in fh:
            sid, orf_ids = line.rstrip("\n").split("\t")
            if sid in wanted:
                orfs_of[sid] = orf_ids.split(",")
    plasmids = {o.rsplit("|", 1)[0] for ms in orfs_of.values() for o in ms}

    genes = collections.defaultdict(list)
    with open(run / "06_annotation_tables/plasmid_annotation.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["plasmid_id"] in plasmids:
                genes[r["plasmid_id"]].append(r)
    for g in genes.values():
        g.sort(key=lambda r: int(r["start"]))

    islands = collections.defaultdict(list)
    ctx = run / "12_context_and_structure"
    for r in read_tsv(ctx / "integrons.tsv"):
        if r["plasmid_id"] in plasmids:
            islands[r["plasmid_id"]].append(
                {"name": "integron", "start": int(r["start"]), "end": int(r["end"])})
    for r in read_tsv(ctx / "is_elements.tsv"):
        if r["plasmid_id"] in plasmids:
            islands[r["plasmid_id"]].append(
                {"name": "is_element", "start": int(r["start"]), "end": int(r["end"])})
    defence_orfs = {r["orf_id"] for r in read_tsv(ctx / "defence_systems.tsv")
                    if r.get("orf_id")}
    size = {r["plasmid_id"]: int(r["size_bp"]) for r in read_tsv(master)
            if r["plasmid_id"] in plasmids}
    # S0 removed a terminal repeat from some circular records; their genes wrap at the
    # molecule's length, not at size_bp.
    repeats = run / "01_analysis_set/terminal_repeats.tsv"
    if repeats.exists():
        size.update({r["plasmid_id"]: int(r["molecule_bp"]) for r in read_tsv(repeats)
                     if r["plasmid_id"] in plasmids})
    return fam_members, orfs_of, genes, islands, defence_orfs, size


def occurrences(members, orfs_of, max_rows):
    """One focal ORF per plasmid, plasmids in id order, capped."""
    first = {}
    for sid in members:
        for orf in orfs_of.get(sid, []):
            first.setdefault(orf.rsplit("|", 1)[0], orf)
    picked = [first[p] for p in sorted(first)]
    return picked[:max_rows], len(picked)


def draw_family(fid, members, orfs_of, genes, islands, defence_orfs, size, window,
                max_rows, out):
    focal_orfs, n_total = occurrences(members, orfs_of, max_rows)
    fig_h = 0.55 * len(focal_orfs) + 1.2
    fig, ax = plt.subplots(figsize=(13, fig_h))
    member_set = set(members)
    orf_to_seq = {o: s for s in members for o in orfs_of.get(s, [])}

    for row, focal in enumerate(focal_orfs):
        pid = focal.rsplit("|", 1)[0]
        g = genes[pid]
        i = next(k for k, r in enumerate(g) if r["orf_id"] == focal)
        f = g[i]
        # Turn the row so the focal gene points right; x is bp from the focal start. A
        # focal gene across the origin ends past the molecule's length, as below.
        flip = f["strand"] in ("-1", "-")
        fs, fe = int(f["start"]), int(f["end"])
        if fe < fs:
            fe += size[pid]
        anchor = fe if flip else fs
        y = -row
        for r in g[max(0, i - window): i + window + 1]:
            s, e = int(r["start"]), int(r["end"])
            if e < s:                      # origin-spanning: continue past the end
                e += size[pid]
            x0, x1 = ((anchor - e, anchor - s) if flip else (s - anchor, e - anchor))
            forward = (r["strand"] in ("1", "+")) != flip
            gene = {"start": s, "end": e}
            for isl in overlapping_islands(gene, islands.get(pid, [])):
                ax.add_patch(Rectangle((x0, y - 0.32), x1 - x0, 0.64,
                                       color=COLOURS[isl["name"]], alpha=0.35, lw=0))
            if r["orf_id"] in defence_orfs:
                ax.add_patch(Rectangle((x0, y - 0.32), x1 - x0, 0.64,
                                       color=COLOURS["defence"], alpha=0.35, lw=0))
            if r["orf_id"] == focal or orf_to_seq.get(r["orf_id"]) in member_set:
                colour = COLOURS["focal"]
            elif r["functional_class"] in DARK_CLASSES:
                colour = COLOURS["dark"]
            elif r["functional_class"] == "NOT_SEARCHED":
                # Never searched (outside every selected family): neither dark nor named.
                colour = COLOURS["not_searched"]
            else:
                colour = COLOURS["annotated"]
            length = x1 - x0
            head = min(abs(length) * 0.3, 150)
            start, dx = (x0, length) if forward else (x1, -length)
            ax.add_patch(FancyArrow(start, y, dx, 0, width=0.36, head_width=0.5,
                                    head_length=head, length_includes_head=True,
                                    color=colour, lw=0))
            if colour == COLOURS["annotated"] and r.get("annot_label"):
                # Swiss-Prot and nr labels are whole NCBI titles ('E0SIS4.1 RecName:
                # Full=... [organism]'); the product name is what fits on an arrow.
                label = parse_ncbi_title(r["annot_label"])["product"]
                ax.text((x0 + x1) / 2, y + 0.3, label[:22], ha="center",
                        va="bottom", fontsize=6)
        ax.annotate(pid, xy=(0, y), xycoords=("axes fraction", "data"), xytext=(-5, 0),
                    textcoords="offset points", ha="right", va="center", fontsize=7)

    ax.autoscale_view()
    ax.set_ylim(-len(focal_orfs) + 0.3, 1.0)
    ax.set_yticks([])
    ax.set_xlabel("bp from the start of the family's gene")
    for spine in ("left", "right", "top"):
        ax.spines[spine].set_visible(False)
    handles = [Rectangle((0, 0), 1, 1, color=c, alpha=a) for c, a in
               [(COLOURS["focal"], 1), (COLOURS["dark"], 1), (COLOURS["annotated"], 1),
                (COLOURS["not_searched"], 1), (COLOURS["defence"], .35),
                (COLOURS["integron"], .35), (COLOURS["is_element"], .35)]]
    ax.legend(handles, ["this family", "dark", "annotated", "not searched",
                        "defence system", "integron", "IS element"], loc="upper right",
              fontsize=7, ncol=7, frameon=False, bbox_to_anchor=(1, 1.12))
    ax.set_title(f"{fid}: {len(focal_orfs)} of {n_total} plasmids shown, "
                 f"+/-{window} genes", fontsize=9, loc="left")
    fig.tight_layout()
    path = out / f"{fid.replace(':', '_')}.pdf"
    fig.savefig(path)
    fig.savefig(path.with_suffix(".svg"))
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True, type=pathlib.Path)
    ap.add_argument("--family", action="append", default=[])
    ap.add_argument("--families-file", type=pathlib.Path)
    ap.add_argument("--window", type=int, default=3)
    ap.add_argument("--max-rows", type=int, default=30)
    ap.add_argument("--out", type=pathlib.Path)
    ap.add_argument("--master", type=pathlib.Path,
                    default=pathlib.Path("data/plasmidscope_primary/analysis_set.tsv"),
                    help="plasmid sizes, for genes that span the origin")
    a = ap.parse_args()
    families = list(a.family)
    if a.families_file:
        families += [l.strip() for l in open(a.families_file) if l.strip()]
    if not families:
        ap.error("give at least one --family or --families-file")
    out = a.out or a.run / "15_report/figures/neighbourhoods"
    out.mkdir(parents=True, exist_ok=True)
    fam_members, orfs_of, genes, islands, defence_orfs, size = load(
        a.run, set(families), a.master)
    for fid in families:
        print(draw_family(fid, fam_members[fid], orfs_of, genes, islands, defence_orfs,
                          size, a.window, a.max_rows, out))


if __name__ == "__main__":
    main()
