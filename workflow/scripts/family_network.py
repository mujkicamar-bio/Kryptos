"""Stage 5b: the family network - 50%-identity clusters linked by sequence similarity.

Method and its source are in src/plasmidann/network.py (Durairaj et al., Nature 2023).
The outputs are two flat tables that Cosmograph and Cytoscape both read directly, plus a
summary written in the terms of the paper, so the numbers can be set side by side:

  network_nodes.tsv    one row per 50%-identity cluster, with its annotation state
  network_edges.tsv    source, target, evalue, weight
  network_summary.tsv  metric, value

Nodes are the clusters that hold at least one small-plasmid protein - the clusters
protein_families.tsv reports at this resolution. Clusters of large-plasmid proteins alone
are left out. `scope` is the family's small-only/mixed and known/unknown call from
protein_families.tsv.

This node set is a SUBSET, and share_dark_connected_to_bright is therefore a lower bound: a
dark small-plasmid cluster whose bright relative sits only on large plasmids has no edge to
it here. The subset is deliberate rather than an oversight. Those large-only clusters were
mostly never searched (NOT_SEARCHED: the cascade annotates only families holding a small-
plasmid protein), so their brightness is unknown and adding them would add nodes whose
state cannot be read. Durairaj et al. also built their network over a subset (UniRef50
clusters with pLDDT > 90), and the summary states the node scope beside the numbers.

MMseqs2 runs at its own default search sensitivity: Durairaj et al. do not state theirs,
and a value we chose would be a parameter with no source.
"""
import collections
import csv
import math
import pathlib
import subprocess

import _ctx  # noqa: F401

from plasmidann import scratch
from plasmidann.fasta import iter_fasta
from plasmidann.cascade import is_informative
from plasmidann.network import communities, edges_from_hits, node_brightness

cfg = snakemake.params.network
out_nodes = pathlib.Path(snakemake.output.nodes)
work = scratch.scratch_dir(out_nodes.parent, "network_tmp")

# ---- nodes: the 50%-identity clusters with a small-plasmid member ---------------------
family_of, node_family = {}, {}
with open(snakemake.input.families, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["family_resolution"] == snakemake.params.primary:
            for m in r["members"].split(","):
                family_of[m] = (r["family_id"], r["family_class"])
        if r["family_resolution"] == snakemake.params.node_resolution:
            node_family[r["representative"]] = r

members = collections.defaultdict(list)
with open(snakemake.input.clusters) as fh:
    for line in fh:
        rep, mem = line.rstrip("\n").split("\t")
        if rep in node_family:
            members[rep].append(mem)

reps_path = work / "node_reps.fasta"
with open(reps_path, "w") as out:
    for sid, seq in iter_fasta([snakemake.input.reps]):
        if sid in members:
            out.write(f">{sid}\n{seq}\n")

annotation = {}
with open(snakemake.input.prot, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        annotation[r["seq_id"]] = r

dark = {line.strip() for line in open(snakemake.input.dark_ids) if line.strip()}

n_orfs, plasmids = collections.Counter(), collections.defaultdict(set)
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, orf_ids = line.rstrip("\n").split("\t")
        for orf_id in orf_ids.split(","):
            n_orfs[sid] += 1
            plasmids[sid].add(orf_id.rsplit("|", 1)[0])

# ---- edges: all-against-all over the representatives ---------------------------------
hits_path = work / "hits.tsv"
subprocess.run(
    f"mmseqs easy-search {reps_path} {reps_path} {hits_path} "
    f"{work / 'tmp'} -e {cfg['max_evalue']} -c 0 --threads {snakemake.threads} "
    f"--format-output query,target,evalue,qcov,tcov -v 1",
    shell=True, check=True)
with open(hits_path) as fh:
    hits = [dict(zip(("query", "target", "evalue", "qcov", "tcov"), line.split()))
            for line in fh]
edges = edges_from_hits(hits, max_out=cfg["max_out_edges"], min_cov=cfg["min_cov"],
                        max_evalue=cfg["max_evalue"])
community = communities(edges, list(members), seed=snakemake.params.seed)

# ---- write nodes and edges -----------------------------------------------------------
neighbours = collections.defaultdict(set)
for a, b in edges:
    neighbours[a].add(b)
    neighbours[b].add(a)

node_rows = {}
for rep, mem in members.items():
    rows = [annotation.get(m, {}) for m in mem]
    labels = collections.Counter(r.get("annot_label") for r in rows
                                 if is_informative(r.get("annot_label")))
    bright = node_brightness(rows)
    node_rows[rep] = {
        "node_id": rep,
        "n_members": len(mem),
        "n_orfs": sum(n_orfs[m] for m in mem),
        "n_plasmids": len(set().union(*(plasmids[m] for m in mem))),
        "brightness": round(bright, 4),
        "dark": int(bright <= cfg["dark_brightness"]),
        "dark_fraction": round(sum(m in dark for m in mem) / len(mem), 4),
        "family_id": family_of.get(rep, ("", ""))[0],
        "scope": node_family[rep]["scope"],
        "label": labels.most_common(1)[0][0] if labels else "dark",
        "community_id": community[rep],
        "degree": len(neighbours[rep]),
    }

with open(out_nodes, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(next(iter(node_rows.values()))), delimiter="\t")
    w.writeheader()
    w.writerows(node_rows.values())

with open(snakemake.output.edges, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["source", "target", "evalue", "weight"])
    for (a, b), ev in sorted(edges.items()):
        w.writerow([a, b, ev, round(-math.log10(ev), 3) if ev > 0 else 300.0])

# ---- summary, in the paper's terms ---------------------------------------------------
dark_nodes = [n for n, r in node_rows.items() if r["dark"]]
dark_connected = [n for n in dark_nodes if neighbours[n]]
dark_to_bright = [n for n in dark_connected
                  if any(not node_rows[m]["dark"] for m in neighbours[n])]
# Dark singleton families (primary resolution): is the protein alone because nothing
# resembles it, or because the family thresholds were not met?
dark_singletons = [m for m, (_, cls) in family_of.items()
                   if cls == "ORPHAN" and m in dark]
singletons_with_edge = [m for m in dark_singletons if neighbours.get(m)]
n_comm = collections.Counter(r["community_id"] for r in node_rows.values())


def share(part, whole):
    return round(len(part) / len(whole), 4) if whole else ""


summary = [
    ("node_scope", f"{snakemake.params.node_resolution} clusters holding a small-plasmid "
                   "protein; share_dark_connected_to_bright is a lower bound"),
    ("nodes", len(node_rows)),
    ("edges", len(edges)),
    ("communities_2plus", sum(1 for c in n_comm.values() if c >= 2)),
    ("dark_nodes", len(dark_nodes)),
    ("dark_nodes_connected", len(dark_connected)),
    ("share_dark_nodes_connected", share(dark_connected, dark_nodes)),
    ("dark_connected_to_bright", len(dark_to_bright)),
    ("share_dark_connected_to_bright", share(dark_to_bright, dark_connected)),
    ("dark_family_singletons", len(dark_singletons)),
    ("dark_family_singletons_with_edge", len(singletons_with_edge)),
    ("mmseqs_sensitivity", "MMseqs2 default"),
    ("min_cov_either_protein", cfg["min_cov"]),
    ("max_evalue", cfg["max_evalue"]),
    ("max_out_edges", cfg["max_out_edges"]),
    ("dark_brightness", cfg["dark_brightness"]),
]
with open(snakemake.output.summary, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["metric", "value"])
    w.writerows(summary)

print("family_network: " + " ".join(f"{k}={v}" for k, v in summary[1:11]))
scratch.release(work)
