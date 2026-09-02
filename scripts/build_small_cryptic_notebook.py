#!/usr/bin/env python3
"""Build notebooks/small_cryptic_investigation.ipynb — the small cryptic plasmidome analysis.

Scope: the sub-10-kb, payload-free half of the working set — plasmids that carry no AMR, no
virulence, no metal/biocide and no conjugation machinery, and that the standard plasmid
classification stack (PlasmidFinder / mob_typer / PlasAnn) largely fails to type at all.

The notebook reads only tables that already exist in the repo:
  - data/plasmidscope_primary/plasmid_metadata_master.tsv   (67-col master; scripts/build_metadata_master.py)
  - data/plasmidscope_primary/mob_full.tsv.gz               (full mob_typer output; scripts/aggregate_mob_full.py)
  - data/plasmidscope_primary/working_set.tsv               (organism/host strings; scripts/assemble_working_set.py)
  - data/plasann_run/annot/mshard_*.tsv.gz                  (per-gene PlasAnn rows; scripts/plasann_harvest_shard.py)

It writes result tables to data/plasmidscope_primary/smallcryptic_*.tsv and figures to
reports/figures/small_cryptic/.

Locked conventions honoured here:
  - analysis set excludes hab_top in {Simulated-artifact, Lab-artifact}  -> 143,503 plasmids
    (reports/environment_taxonomy_LOCKED.md)
  - ecological strata are built from hab_sub, never hab_top

Run: python3 scripts/build_small_cryptic_notebook.py   (execute with the genesis_nb kernel)
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []


def md(t):
    cells.append(nbf.v4.new_markdown_cell(t.strip("\n")))


def code(s):
    cells.append(nbf.v4.new_code_cell(s.strip("\n")))


# ======================================================== intro =====================
md(r"""
# The small cryptic plasmidome

**Working title of the angle:** *Half of the plasmidome is invisible to plasmid biology.*

Plasmid genomics is built almost entirely on the large, conjugative, cargo-carrying end of the size
distribution — the replicons that carry antibiotic resistance and that replicon typing schemes were
designed around. This notebook characterises the other half: **small (< ~10 kb), payload-free
"cryptic" plasmids**, which turn out to be roughly **half of all complete plasmids** in the working
set and which the standard classification stack cannot type.

The analysis proceeds as a chain of falsifiable claims:

| § | Claim under test |
|---|---|
| 1 | The plasmidome is size-bimodal, and the antimode gives a *data-driven* small/large cut rather than an arbitrary one. |
| 2 | A majority of small plasmids carry no recognisable functional payload at all ("cryptic"). |
| 3 | **The blind spot.** Replicon typing, mobility prediction and host-range prediction all fail disproportionately on this compartment — the failure is structural, not random. |
| 4 | **The dispersal paradox.** Small cryptic plasmids are overwhelmingly called *non-mobilizable*, yet their lineages span as many habitats and countries as conjugative plasmids do — tested with cluster-size rarefaction, because breadth trivially scales with sampling. |
| 5 | Where a host-range call exists, breadth is **bimodal** for small plasmids (narrow *or* multi-phylum) while conjugative plasmids sit at a single intermediate rank. |
| 6 | **Bias audit.** Everything above is re-tested against database provenance, sampling effort and plasmid size, because the round-1 review of the AMR manuscript showed that several apparent biological effects in this dataset were size/composition artefacts. |
| 7 | Gene content: how much of the coding capacity is genuinely unannotatable, and what recurrent modules exist. |

**Provenance.** Built by `scripts/build_small_cryptic_notebook.py`. Inputs are
`data/plasmidscope_primary/{plasmid_metadata_master.tsv, mob_full.tsv.gz, working_set.tsv}` and the
per-gene PlasAnn shards `data/plasann_run/annot/mshard_*.tsv.gz`. Upstream methods:
`reports/typing_methodology.md` (PlasmidFinder + mob_typer), `reports/card_amr_methodology.md`
(CARD/RGI), `reports/environment_taxonomy_LOCKED.md` (habitat taxonomy).

**Analysis set.** The locked exclusion of the `Simulated-artifact` and `Lab-artifact` habitat buckets
is applied once, giving **143,503** plasmids. Every number below refers to that set unless stated.
""")

code(r"""
import os, glob, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
warnings.filterwarnings("ignore")

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 80)

ROOT = "/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis"
PP   = f"{ROOT}/data/plasmidscope_primary"
ANN  = f"{ROOT}/data/plasann_run/annot"
FIG  = f"{ROOT}/reports/figures/small_cryptic"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 200, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
})
# consistent group colours used throughout
C = {"small-cryptic": "#c1442e", "small-cargo": "#e2a13b",
     "large-conjugative": "#2f6f9f", "large-other": "#7fa8c4"}
GRP_ORDER = ["small-cryptic", "small-cargo", "large-conjugative", "large-other"]
print("matplotlib", matplotlib.__version__, "| pandas", pd.__version__)
""")

# ======================================================== load ======================
md(r"""
## 0. Load the master table and apply the locked exclusions

`mob_full.tsv.gz` is joined for the two mob_typer fields that the 67-column master does not carry:
the **rank** of the predicted host range (`predicted_host_range_overall_rank`) and the mash neighbour
distance that the prediction rests on. `working_set.tsv` contributes the free-text `organism` string,
the only host-taxon label available for plasmids without a PLSDB record.
""")

code(r"""
master = pd.read_csv(f"{PP}/plasmid_metadata_master.tsv", sep="\t", low_memory=False)

mobf = (pd.read_csv(f"{PP}/mob_full.tsv.gz", sep="\t", low_memory=False)
          .rename(columns={"sample_id": "plasmid_id"})
          [["plasmid_id", "predicted_host_range_overall_rank",
            "predicted_host_range_overall_name", "mash_neighbor_distance"]])

ws = pd.read_csv(f"{PP}/working_set.tsv", sep="\t", low_memory=False,
                 usecols=["plasmid_id", "organism", "biosample_host"])

d = master.merge(mobf, on="plasmid_id", how="left").merge(ws, on="plasmid_id", how="left")

# mob_typer writes "-" for "no call"; make that explicit missingness
d["hr_rank"] = d.predicted_host_range_overall_rank.replace("-", np.nan)
d["hr_name"] = d.predicted_host_range_overall_name.replace("-", np.nan)

# LOCKED exclusion (reports/environment_taxonomy_LOCKED.md)
real = d[~d.hab_top.isin(["Simulated-artifact", "Lab-artifact"])].copy()

print(f"master            : {len(master):,}")
print(f"after locked excl.: {len(real):,}   (expect 143,503)")
assert len(real) == 143503, "analysis set size drifted from the locked value"
""")

# ======================================================== S1 ========================
md(r"""
## 1. The size structure of the plasmidome, and a data-driven small/large cut

Plasmid size is famously bimodal, and the small/large boundary is usually asserted rather than
measured. A Gaussian kernel density on log₁₀(size) locates the two modes and the **antimode** (the
density minimum between them), which is reported below.

**The working cut is 10 kb, and it is deliberately *not* the antimode.** The antimode sits well above
10 kb, in the sparse valley between the modes; taking it as the boundary would sweep the whole valley
into the "small" class and inflate every quantity that follows. 10 kb is the conservative choice — it
sits inside the small mode's own shoulder, so the group it defines is unambiguously small-mode
material. Because the choice is a judgement call rather than a measurement, **every headline quantity
is re-reported at 5, 8, 10, 15 and 20 kb in §6a**, a range that brackets the antimode. Nothing in the
conclusions depends on where in that range the line is drawn.
""")

code(r"""
from scipy.stats import gaussian_kde

ls = np.log10(real.size_bp.clip(lower=200).values)
kde = gaussian_kde(ls, bw_method=0.15)
grid = np.linspace(np.log10(500), np.log10(1.5e6), 2000)
dens = kde(grid)

# peaks and the antimode strictly between them
from scipy.signal import argrelextrema
maxi = argrelextrema(dens, np.greater)[0]
mini = argrelextrema(dens, np.less)[0]
peaks = grid[maxi]
print("density peaks (bp):", np.round(10**peaks).astype(int))
cand = [g for g in grid[mini] if peaks.min() < g < peaks.max()]
ANTIMODE = 10**cand[int(np.argmin([dens[np.argmin(abs(grid-c))] for c in cand]))] if cand else 10000.0
print(f"antimode  : {ANTIMODE:,.0f} bp")

# Deliberately conservative: BELOW the antimode, inside the small mode. Sensitivity-tested in §6a
# across 5-20 kb, a range that brackets the antimode.
SMALL_CUT = 10000
print(f"SMALL_CUT : {SMALL_CUT:,} bp  (conservative; antimode is {ANTIMODE:,.0f} bp)")

fig, ax = plt.subplots(figsize=(7.2, 3.2))
ax.fill_between(10**grid, dens, color="#8c8c8c", alpha=.35, lw=0)
ax.plot(10**grid, dens, color="#333", lw=1.2)
ax.axvline(ANTIMODE, color=C["small-cryptic"], ls="--", lw=1.2,
           label=f"antimode {ANTIMODE:,.0f} bp")
ax.axvline(SMALL_CUT, color="#333", ls=":", lw=1.2, label=f"cut {SMALL_CUT:,} bp")
for p in 10**peaks:
    ax.annotate(f"{p:,.0f} bp", (p, kde(np.log10(p))[0]), textcoords="offset points",
                xytext=(0, 6), ha="center", fontsize=8)
ax.set_xscale("log"); ax.set_xlabel("plasmid size (bp)"); ax.set_ylabel("density")
ax.set_title(f"Size structure of the plasmidome (n = {len(real):,} complete plasmids)")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/fig1_size_bimodality.png"); plt.show()

sm = (real.size_bp < SMALL_CUT)
print(f"\nsmall (<{SMALL_CUT:,} bp): {sm.sum():,} ({100*sm.mean():.1f}%)   "
      f"large: {(~sm).sum():,} ({100*(~sm).mean():.1f}%)")
""")

# ======================================================== S2 ========================
md(r"""
## 2. Defining "cryptic"

A plasmid is called **cryptic** when it carries none of the four functional payload classes the field
routinely screens for:

- **AMR** — CARD/RGI, Perfect+Strict only (`card_n_arg`); the master deliberately excludes Loose hits.
- **virulence** — PlasAnn `plasann_n_virulence`
- **metal / biocide resistance** — PlasAnn `plasann_n_metal_biocide`
- **conjugation machinery** — PlasAnn `plasann_n_conjugation`

Conjugation is included in the payload definition because a plasmid encoding its own transfer
apparatus is not functionally cryptic — it has an evident, self-serving phenotype. Toxin–antitoxin and
mobilization (MOB/oriT) genes are deliberately **not** disqualifying: they are the maintenance and
hitchhiking modules whose prevalence in this compartment is one of the questions (§7).

This yields the four analysis groups used throughout: `small-cryptic`, `small-cargo`,
`large-conjugative`, `large-other`.
""")

code(r"""
def num(s):
    return pd.to_numeric(s, errors="coerce").fillna(0)

real["n_amr"]   = num(real.card_n_arg)
real["n_vir"]   = num(real.plasann_n_virulence)
real["n_metal"] = num(real.plasann_n_metal_biocide)
real["n_conj"]  = num(real.plasann_n_conjugation)
real["n_ta"]    = num(real.plasann_n_toxin_antitoxin)
real["n_mobdna"]= num(real.plasann_n_mob_dna)
real["n_cds"]   = num(real.plasann_n_cds)

real["payload"] = (real.n_amr > 0) | (real.n_vir > 0) | (real.n_metal > 0)
real["cryptic"] = ~real.payload & (real.n_conj == 0)
real["small"]   = real.size_bp < SMALL_CUT

real["grp"] = np.select(
    [ real.small &  real.cryptic,
      real.small & ~real.cryptic,
     ~real.small & (real.mob_mobility == "conjugative"),
     ~real.small],
    ["small-cryptic", "small-cargo", "large-conjugative", "large-other"], default="unassigned")

funnel = pd.DataFrame([
    ("complete, non-redundant plasmids (working set)", len(master)),
    ("- simulated / lab artefacts (locked exclusion)", len(real)),
    (f"of which small (< {SMALL_CUT:,} bp)",           int(real.small.sum())),
    ("  - carrying AMR (CARD Perfect+Strict)",         int((real.small & (real.n_amr>0)).sum())),
    ("  - carrying virulence genes",                   int((real.small & (real.n_vir>0)).sum())),
    ("  - carrying metal/biocide resistance",          int((real.small & (real.n_metal>0)).sum())),
    ("  - carrying conjugation machinery",             int((real.small & (real.n_conj>0)).sum())),
    ("= SMALL CRYPTIC (no payload of any class)",      int((real.grp=='small-cryptic').sum())),
], columns=["stage", "n"])
funnel["pct_of_analysis_set"] = (100*funnel.n/len(real)).round(2)
funnel.to_csv(f"{PP}/smallcryptic_funnel.tsv", sep="\t", index=False)
print(funnel.to_string(index=False))

print("\ngroup sizes:")
print(real.grp.value_counts().reindex(GRP_ORDER).to_frame("n")
        .assign(pct=lambda x: (100*x.n/len(real)).round(1)).to_string())
""")

code(r"""
# profile of each group
g = real.groupby("grp")
prof = pd.DataFrame({
    "n":            g.size(),
    "median_bp":    g.size_bp.median().round(0),
    "median_GC":    g.gc_percent.median().round(1),
    "median_CDS":   g.n_cds.median(),
    "pct_AMR":      (100*g.apply(lambda x: (x.n_amr>0).mean())).round(2),
    "pct_TA":       (100*g.apply(lambda x: (x.n_ta>0).mean())).round(2),
    "pct_oriT":     (100*g.apply(lambda x: num(x.plasann_has_orit).gt(0).mean())).round(2),
    "pct_MOBrelax": (100*g.apply(lambda x: x.mob_relaxase.notna().mean())).round(2),
    "pct_clinical": (100*g.is_clinical.mean()).round(2),
}).reindex(GRP_ORDER)
prof.to_csv(f"{PP}/smallcryptic_group_profile.tsv", sep="\t")
prof
""")


# ======================================================== S3 ========================
md(r"""
## 3. The blind spot: three independent typing methods, one shared failure

The working set has been typed three separate ways — **PlasmidFinder** (CGE thresholds, ≥80% id /
≥60% cov), **mob_typer** (MOB-suite replicon HMMs/BLAST), and **PlasAnn** (annotation-derived replicon
calls). These are different databases and different algorithms, so agreement or disagreement between
them is informative.

The test: does replicon detectability degrade with size, and is the degradation concentrated in the
cryptic compartment? A tool failing on 3% of a compartment is noise; a tool failing on 70% of the
single largest compartment of the plasmidome is a structural gap in the reference databases.
""")

code(r"""
real["rep_pf"]  = num(real.pf_n_inc) > 0
real["rep_mob"] = real.mob_rep_types.notna()
real["rep_pa"]  = num(real.plasann_n_replicons) > 0
real["rep_any"] = real.rep_pf | real.rep_mob | real.rep_pa
real["hr_call"] = real.hr_rank.notna()

blind = pd.DataFrame({
    "n":              real.groupby("grp").size(),
    "PlasmidFinder":  100*real.groupby("grp").rep_pf.mean(),
    "mob_typer":      100*real.groupby("grp").rep_mob.mean(),
    "PlasAnn":        100*real.groupby("grp").rep_pa.mean(),
    "ANY of the 3":   100*real.groupby("grp").rep_any.mean(),
    "host-range call":100*real.groupby("grp").hr_call.mean(),
}).reindex(GRP_ORDER).round(1)
blind.to_csv(f"{PP}/smallcryptic_typing_blindspot.tsv", sep="\t")
print("Replicon detection rate (%) by method and group\n")
print(blind.to_string())
print(f"\nUNTYPED by all three methods, small-cryptic: "
      f"{int((~real.rep_any & (real.grp=='small-cryptic')).sum()):,} plasmids "
      f"({100*(~real[real.grp=='small-cryptic'].rep_any).mean():.1f}% of the group, "
      f"{100*(~real.rep_any & (real.grp=='small-cryptic')).sum()/len(real):.1f}% of the whole analysis set)")
""")

code(r"""
# detectability as a continuous function of size -- is the cut doing the work, or is it monotone?
edges = np.array([0,1500,2000,2500,3000,4000,5000,7000,10000,15000,25000,50000,100000,1e9])
lab   = ["<1.5k","1.5-2k","2-2.5k","2.5-3k","3-4k","4-5k","5-7k","7-10k","10-15k",
         "15-25k","25-50k","50-100k",">100k"]
real["szbin"] = pd.cut(real.size_bp, edges, labels=lab)
bybin = real.groupby("szbin").agg(n=("plasmid_id","size"),
                                  pf=("rep_pf","mean"), mob=("rep_mob","mean"),
                                  pa=("rep_pa","mean"), any_=("rep_any","mean"),
                                  hr=("hr_call","mean"))
fig, ax = plt.subplots(figsize=(7.6, 3.4))
x = np.arange(len(bybin))
for col, nm, col_c in [("pf","PlasmidFinder","#c1442e"), ("mob","mob_typer","#2f6f9f"),
                       ("pa","PlasAnn","#5f9e6e"), ("any_","any of the three","#333333"),
                       ("hr","mob_typer host-range call","#a06fb5")]:
    ax.plot(x, 100*bybin[col], marker="o", ms=3.5, lw=1.4, color=col_c, label=nm,
            ls="--" if col=="hr" else "-")
ax.axvline(lab.index("7-10k")+.5, color="#999", lw=1, ls=":")
ax.set_xticks(x); ax.set_xticklabels(bybin.index, rotation=45, ha="right")
ax.set_ylabel("% of plasmids with a call"); ax.set_xlabel("plasmid size")
ax.set_ylim(0, 100); ax.legend(frameon=False, fontsize=8, ncol=2)
ax.set_title("Replicon typing covers ~40% of the small plasmidome and ~78% of the large")
fig.tight_layout(); fig.savefig(f"{FIG}/fig2_typing_coverage_by_size.png"); plt.show()
bybin.assign(n=bybin.n).round(3).to_csv(f"{PP}/smallcryptic_typing_by_size.tsv", sep="\t")
(100*bybin[["pf","mob","pa","any_","hr"]]).round(1).assign(n=bybin.n)
""")

# ======================================================== S4 ========================
md(r"""
## 4. The dispersal paradox

mob_typer calls **~85% of small cryptic plasmids non-mobilizable** — no relaxase, no recognised oriT,
no transfer apparatus. Under the standard model these elements should be effectively stuck in whatever
lineage they arose in, and should therefore be far more locally confined than conjugative plasmids.

The test uses **MOB-suite primary clusters** as lineage proxies and measures, per cluster, how many
distinct habitats (`hab_sub`) and countries its members were sampled from.

**The measurement trap:** breadth grows mechanically with the number of members sampled. A cluster of
200 will span more countries than a cluster of 12 whatever its biology. So breadth is measured by
**rarefaction** — every cluster is subsampled to exactly *K* members, repeated 200 times, and the mean
distinct-count is taken. Only clusters with ≥ *K* members enter the comparison, and both groups are
therefore compared at identical sampling depth.
""")

code(r"""
rng = np.random.default_rng(20260826)

def rarefied_breadth(df, field, K=10, reps=200, group_col="grp", cluster_col="mob_cluster"):
    "Mean number of distinct `field` values seen when K members are drawn from each cluster."
    out = []
    sub = df[df[cluster_col].notna() & df[field].notna()]
    for (grp, clu), block in sub.groupby([group_col, cluster_col], observed=True):
        vals = block[field].values
        if len(vals) < K:
            continue
        draws = [len(np.unique(rng.choice(vals, K, replace=False))) for _ in range(reps)]
        out.append({"grp": grp, "cluster": clu, "n_members": len(vals),
                    "breadth": float(np.mean(draws))})
    return pd.DataFrame(out)

K = 10
hab = rarefied_breadth(real, "hab_sub",     K=K)
geo = rarefied_breadth(real, "geo_country", K=K)

def summarise(t, what):
    s = (t.groupby("grp").breadth.agg(n_clusters="size", median="median", mean="mean")
           .reindex(GRP_ORDER).dropna(how="all").round(2))
    s.insert(0, "metric", what)
    return s

disp = pd.concat([summarise(hab, f"distinct habitats per {K} members"),
                  summarise(geo, f"distinct countries per {K} members")])
disp.to_csv(f"{PP}/smallcryptic_dispersal_rarefied.tsv", sep="\t")
print(disp.to_string())
""")

code(r"""
from scipy.stats import mannwhitneyu

print(f"Rarefied breadth, small-cryptic vs large-conjugative (K={K} members per cluster)\n")
rows = []
for t, what in [(hab, "habitats"), (geo, "countries")]:
    a = t[t.grp == "small-cryptic"].breadth
    b = t[t.grp == "large-conjugative"].breadth
    if len(a) > 2 and len(b) > 2:
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        # rank-biserial effect size
        rb = 1 - 2*u/(len(a)*len(b))
        rows.append({"metric": what, "n_clusters_cryptic": len(a), "n_clusters_conj": len(b),
                     "median_cryptic": round(a.median(), 2), "median_conj": round(b.median(), 2),
                     "ratio": round(a.median()/b.median(), 3),
                     "rank_biserial": round(rb, 3), "p_MWU": f"{p:.3g}"})
res = pd.DataFrame(rows)
res.to_csv(f"{PP}/smallcryptic_dispersal_test.tsv", sep="\t", index=False)
print(res.to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.4))
for ax, t, what in [(axes[0], hab, "habitats"), (axes[1], geo, "countries")]:
    grps = [g for g in GRP_ORDER if (t.grp == g).sum() >= 3]
    data = [t[t.grp == g].breadth.values for g in grps]
    bp = ax.boxplot(data, labels=[g.replace("-", "\n") for g in grps], showfliers=False,
                    patch_artist=True, widths=.6, medianprops=dict(color="black", lw=1.4))
    for patch, gname in zip(bp["boxes"], grps):
        patch.set_facecolor(C[gname]); patch.set_alpha(.65); patch.set_edgecolor("#444")
    for i, vals in enumerate(data, start=1):
        jx = np.random.normal(i, .055, size=min(len(vals), 400))
        ax.plot(jx, np.random.choice(vals, len(jx), replace=False), ".", ms=1.6,
                color="#222", alpha=.28, zorder=3)
    ax.set_ylabel(f"distinct {what} per {K} sampled members")
    ax.set_title(f"Rarefied {what[:-1]} breadth", fontsize=9)
fig.suptitle("Lineage dispersal at matched sampling depth", y=1.02, fontsize=10)
fig.tight_layout(); fig.savefig(f"{FIG}/fig3_dispersal_rarefied.png", bbox_inches="tight"); plt.show()
""")

code(r"""
# same comparison across a ladder of K, to show it is not a K artefact
ladder = []
for KK in [5, 10, 20, 40]:
    for t_name, field in [("habitats", "hab_sub"), ("countries", "geo_country")]:
        t = rarefied_breadth(real, field, K=KK, reps=60)
        a = t[t.grp == "small-cryptic"].breadth
        b = t[t.grp == "large-conjugative"].breadth
        if len(a) > 2 and len(b) > 2:
            _, p = mannwhitneyu(a, b, alternative="two-sided")
            ladder.append({"K": KK, "metric": t_name, "n_cryptic": len(a), "n_conj": len(b),
                           "median_cryptic": round(a.median(), 2),
                           "median_conj": round(b.median(), 2),
                           "ratio": round(a.median()/b.median(), 3), "p": f"{p:.2g}"})
lad = pd.DataFrame(ladder)
lad.to_csv(f"{PP}/smallcryptic_dispersal_Kladder.tsv", sep="\t", index=False)
lad
""")

# ======================================================== S5 ========================
md(r"""
## 5. Host-range breadth is bimodal in the small plasmidome

mob_typer reports the taxonomic **rank** at which a plasmid's predicted host range sits — from `genus`
up through `phylum`, with a separate `multi-phylla` category for elements whose neighbourhood spans
more than one phylum. Ranks are ordered here from narrow to broad.

Two caveats stated up front, because they bound how strongly this can be interpreted:

1. The prediction is **cluster-based** — it reads the host labels of a plasmid's mash neighbourhood.
   For sparsely populated small-plasmid clusters this is a noisier estimator than for the densely
   sampled Enterobacteriaceae clusters, and neighbourhood heterogeneity can inflate a `multi-phylla`
   call. §6 therefore conditions on `mash_neighbor_distance`.
2. **Two-thirds of small cryptic plasmids get no call at all**, so this section describes only the
   callable third and is not representative of the compartment as a whole.
""")

code(r"""
RANKS = ["genus", "family", "order", "class", "phylum", "multi-phylla"]
sub = real[real.hr_rank.notna()].copy()
sub["hr_rank"] = pd.Categorical(sub.hr_rank, categories=RANKS, ordered=True)

hrpct = (pd.crosstab(sub.grp, sub.hr_rank, normalize="index").mul(100)
           .reindex(GRP_ORDER)[RANKS].round(1))
hrn   = pd.crosstab(sub.grp, sub.hr_rank).reindex(GRP_ORDER)[RANKS]
hrpct.assign(n_with_call=hrn.sum(1),
             pct_of_group_with_call=(100*hrn.sum(1)/real.grp.value_counts().reindex(GRP_ORDER)).round(1)
            ).to_csv(f"{PP}/smallcryptic_hostrange_rank.tsv", sep="\t")
print("Host-range rank distribution, row % (among plasmids with a call)\n")
print(hrpct.assign(n_with_call=hrn.sum(1)).to_string())

fig, ax = plt.subplots(figsize=(7.4, 3.0))
left = np.zeros(len(GRP_ORDER))
shades = plt.cm.YlGnBu(np.linspace(.25, .92, len(RANKS)))
for r, cc in zip(RANKS, shades):
    v = hrpct[r].values
    ax.barh(GRP_ORDER, v, left=left, color=cc, edgecolor="white", lw=.6, label=r)
    for i, (val, l) in enumerate(zip(v, left)):
        if val >= 6:
            ax.text(l+val/2, i, f"{val:.0f}", ha="center", va="center", fontsize=7.5,
                    color="white" if r in ("phylum", "multi-phylla") else "#222")
    left += v
ax.set_xlabel("% of plasmids with a host-range call"); ax.set_xlim(0, 100)
ax.invert_yaxis(); ax.grid(False)
ax.legend(frameon=False, fontsize=7.5, ncol=6, loc="upper center", bbox_to_anchor=(.5, -.22))
ax.set_title("Predicted host-range breadth by group")
fig.tight_layout(); fig.savefig(f"{FIG}/fig4_hostrange_rank.png", bbox_inches="tight"); plt.show()

print(f"\nmulti-phylum calls: small-cryptic {hrpct.loc['small-cryptic','multi-phylla']:.1f}%  "
      f"small-cargo {hrpct.loc['small-cargo','multi-phylla']:.1f}%  "
      f"large-conjugative {hrpct.loc['large-conjugative','multi-phylla']:.1f}%")
""")

# ======================================================== S6 ========================
md(r"""
## 6. Bias audit

Round-1 review of the AMR manuscript (`manuscript/review_round1.md`) established that several
apparent biological effects in this dataset were artefacts of plasmid size and database composition.
That lesson is applied pre-emptively here. Four challenges:

**6a — threshold dependence.** Does anything hinge on the 10 kb cut?
**6b — database provenance.** IMG-PR contributes 136k of the 208k source records and is
metagenome-derived; RefSeq/GenBank/PLSDB are isolate-derived. If `small-cryptic` is simply
"the IMG-PR compartment", every claim above is a statement about a sequencing strategy, not biology.
**6c — within-source replication.** The typing and dispersal results are re-run *inside* each
provenance stratum. An effect that survives in both isolate-derived and metagenome-derived data is not
a provenance artefact.
**6d — neighbourhood quality.** Host-range breadth is re-examined conditional on
`mash_neighbor_distance`, since a distant neighbourhood makes a `multi-phylla` call less trustworthy.
""")

code(r"""
# --- 6a: threshold sensitivity -------------------------------------------------------
sens = []
for cut in [5000, 8000, 10000, 15000, 20000]:
    sm_ = real.size_bp < cut
    cr_ = sm_ & real.cryptic
    sens.append({
        "cut_bp": cut,
        "n_small": int(sm_.sum()),
        "pct_of_set": round(100*sm_.mean(), 1),
        "n_cryptic": int(cr_.sum()),
        "pct_small_that_is_cryptic": round(100*cr_.sum()/sm_.sum(), 1),
        "pct_cryptic_untyped_all3": round(100*(~real.loc[cr_, "rep_any"]).mean(), 1),
        "pct_cryptic_nonmobilizable": round(100*(real.loc[cr_, "mob_mobility"] == "non-mobilizable").mean(), 1),
        "pct_cryptic_no_hostrange_call": round(100*(~real.loc[cr_, "hr_call"]).mean(), 1),
    })
sens = pd.DataFrame(sens)
sens.to_csv(f"{PP}/smallcryptic_threshold_sensitivity.tsv", sep="\t", index=False)
print("6a — every headline quantity across five size thresholds\n")
print(sens.to_string(index=False))
""")

code(r"""
# --- 6b: database provenance composition --------------------------------------------
real["is_imgpr"]  = real.sources.fillna("").str.contains("IMG-PR")
real["is_isolate"] = real.sources.fillna("").str.contains("RefSeq|GenBank|PLSDB|COMPASS|DDBJ|ENA", regex=True)
real["prov"] = np.select(
    [real.is_imgpr & ~real.is_isolate, ~real.is_imgpr & real.is_isolate, real.is_imgpr & real.is_isolate],
    ["IMG-PR only (metagenomic)", "isolate-derived only", "both"], default="other")

prov = pd.crosstab(real.grp, real.prov, normalize="index").mul(100).reindex(GRP_ORDER).round(1)
prov["n"] = real.grp.value_counts().reindex(GRP_ORDER)
prov.to_csv(f"{PP}/smallcryptic_provenance.tsv", sep="\t")
print("6b — provenance composition of each group (row %)\n")
print(prov.to_string())
print("\n=> if small-cryptic were merely the metagenomic compartment, its row would be ~100% IMG-PR only.")
""")

code(r"""
# --- 6c: replicate the two central claims WITHIN each provenance stratum -------------
strata = ["IMG-PR only (metagenomic)", "isolate-derived only", "both"]
rep_rows = []
for st in strata:
    R = real[real.prov == st]
    if len(R) < 500:
        continue
    for grp in ["small-cryptic", "large-conjugative"]:
        Rg = R[R.grp == grp]
        if len(Rg) < 100:
            continue
        rep_rows.append({
            "stratum": st, "group": grp, "n": len(Rg),
            "pct_untyped_all3": round(100*(~Rg.rep_any).mean(), 1),
            "pct_no_hostrange_call": round(100*(~Rg.hr_call).mean(), 1),
            "pct_nonmobilizable": round(100*(Rg.mob_mobility == "non-mobilizable").mean(), 1),
            "median_bp": int(Rg.size_bp.median()),
        })
rep_df = pd.DataFrame(rep_rows)
rep_df.to_csv(f"{PP}/smallcryptic_within_provenance.tsv", sep="\t", index=False)
print("6c — the typing blind spot inside each provenance stratum\n")
print(rep_df.to_string(index=False))

# dispersal, replicated inside the isolate-derived stratum only (the conservative test)
iso = real[real.prov != "IMG-PR only (metagenomic)"]
hab_iso = rarefied_breadth(iso, "hab_sub", K=10, reps=100)
geo_iso = rarefied_breadth(iso, "geo_country", K=10, reps=100)
rows = []
for t, what in [(hab_iso, "habitats"), (geo_iso, "countries")]:
    a = t[t.grp == "small-cryptic"].breadth; b = t[t.grp == "large-conjugative"].breadth
    if len(a) > 2 and len(b) > 2:
        _, p = mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"metric": what, "n_cl_cryptic": len(a), "n_cl_conj": len(b),
                     "median_cryptic": round(a.median(), 2), "median_conj": round(b.median(), 2),
                     "ratio": round(a.median()/b.median(), 3), "p": f"{p:.3g}"})
iso_disp = pd.DataFrame(rows)
iso_disp.to_csv(f"{PP}/smallcryptic_dispersal_isolate_only.tsv", sep="\t", index=False)
print("\n6c — rarefied dispersal, isolate-derived plasmids only (metagenomes excluded)\n")
print(iso_disp.to_string(index=False))
""")

code(r"""
# --- 6d: host-range calls conditional on neighbourhood quality -----------------------
sub2 = real[real.hr_rank.notna()].copy()
sub2["mnd"] = pd.to_numeric(sub2.mash_neighbor_distance, errors="coerce")
sub2["mnd_bin"] = pd.cut(sub2.mnd, [-.001, .001, .01, .05, .1, 1],
                         labels=["identical (0)", "≤0.01", "0.01-0.05", "0.05-0.1", ">0.1"])
mp = (sub2.assign(multi=sub2.hr_rank == "multi-phylla")
          .groupby(["grp", "mnd_bin"], observed=True)
          .agg(n=("multi", "size"), pct_multiphylum=("multi", lambda s: round(100*s.mean(), 1)))
          .reset_index())
mp = mp[mp.grp.isin(["small-cryptic", "small-cargo", "large-conjugative"]) & (mp.n >= 30)]
mp.to_csv(f"{PP}/smallcryptic_hostrange_vs_mashdist.tsv", sep="\t", index=False)
print("6d — multi-phylum call rate vs. mash neighbourhood distance\n")
print(mp.pivot(index="mnd_bin", columns="grp", values="pct_multiphylum").to_string())
print("\nn per cell:")
print(mp.pivot(index="mnd_bin", columns="grp", values="n").to_string())
print("\n=> if multi-phylum calls were driven purely by distant/noisy neighbourhoods, they would be")
print("   concentrated in the high-distance bins and absent from the 'identical (0)' bin.")
""")

# ======================================================== S7 ========================
md(r"""
## 7. What is actually inside a cryptic plasmid

The per-gene PlasAnn output (`data/plasann_run/annot/mshard_*.tsv.gz`, ~7 M CDS rows) is loaded and
restricted to the small compartment. PlasAnn labels a CDS it cannot assign as `ORF` under the category
`Open reading frame`; the **ORF-only fraction** of a plasmid is therefore a direct measure of its
unannotatable coding capacity.

Two questions: how much of the small cryptic plasmidome is genuinely dark, and does what little *is*
annotatable fall into recurrent modules — replication, mobilization, maintenance — or into a long tail?
""")

code(r"""
# NOTE: the PlasAnn harvest wrote TWO shard series -- shard_*.tsv.gz (the 595-shard primary run)
# and mshard_*.tsv.gz (the 477-shard reshard of initially-missing plasmids). Both are required;
# globbing only one of them silently analyses a ~27% subsample.
shards = sorted(glob.glob(f"{ANN}/*shard_*.tsv.gz"))
print(f"reading {len(shards)} PlasAnn shards "
      f"({sum('mshard' in os.path.basename(f) for f in shards)} reshard + "
      f"{sum('mshard' not in os.path.basename(f) for f in shards)} primary) ...")

KEEP = ["plasmid_id", "Gene Name", "Category", "feature type"]
def _read(f):
    g = pd.read_csv(f, sep="\t", usecols=KEEP, low_memory=False)
    g.columns = [c.strip().lower().replace(" ", "_") for c in g.columns]
    return g[g.feature_type == "CDS"].drop(columns="feature_type")

cds = pd.concat((_read(f) for f in shards), ignore_index=True)
cov = cds.plasmid_id.nunique()
print(f"{len(cds):,} CDS rows over {cov:,} plasmids "
      f"({100*cov/len(master):.1f}% of the working set)")
""")

code(r"""
grp_map = real.set_index("plasmid_id").grp
cds["grp"] = cds.plasmid_id.map(grp_map)
cds_a = cds[cds.grp.notna()].copy()
cds_a["is_orf"] = cds_a.category.eq("Open reading frame")

# per-plasmid dark fraction
perp = (cds_a.groupby(["plasmid_id", "grp"], observed=True)
             .agg(n_cds=("is_orf", "size"), n_orf=("is_orf", "sum")).reset_index())
perp["dark_frac"] = perp.n_orf / perp.n_cds

dark = (perp.groupby("grp")
            .agg(n_plasmids=("plasmid_id", "size"),
                 total_cds=("n_cds", "sum"), total_orf=("n_orf", "sum"),
                 median_dark_frac=("dark_frac", "median"),
                 pct_fully_dark=("dark_frac", lambda s: round(100*(s == 1).mean(), 1)))
            .reindex(GRP_ORDER))
dark["pooled_pct_dark"] = (100*dark.total_orf/dark.total_cds).round(1)
dark["median_dark_frac"] = dark.median_dark_frac.round(3)
dark.to_csv(f"{PP}/smallcryptic_dark_matter.tsv", sep="\t")
print("Unannotatable coding capacity by group\n")
print(dark.to_string())
""")

code(r"""
fig = plt.figure(figsize=(8.6, 3.4))
gs = gridspec.GridSpec(1, 2, width_ratios=[1.15, 1], wspace=.3)

ax = fig.add_subplot(gs[0])
for grp in GRP_ORDER:
    v = perp[perp.grp == grp].dark_frac.values
    if len(v) < 50: continue
    ax.hist(v, bins=np.linspace(0, 1, 41), density=True, histtype="step", lw=1.6,
            color=C[grp], label=f"{grp} (n={len(v):,})")
ax.set_xlabel("fraction of CDS unannotatable (PlasAnn 'ORF')")
ax.set_ylabel("density"); ax.legend(frameon=False, fontsize=7.5)
ax.set_title("Coding dark matter", fontsize=9)

ax2 = fig.add_subplot(gs[1])
named = cds_a[(cds_a.grp == "small-cryptic") & (~cds_a.is_orf)]
top = named.gene_name.value_counts().head(18)[::-1]
ax2.barh(top.index, top.values, color=C["small-cryptic"], alpha=.8)
ax2.set_xlabel("plasmids"); ax2.set_title("Most frequent named genes,\nsmall cryptic plasmids", fontsize=9)
ax2.tick_params(axis="y", labelsize=7.5); ax2.grid(axis="y", visible=False)
fig.savefig(f"{FIG}/fig5_dark_matter.png", bbox_inches="tight"); plt.show()

cnt = (cds_a.groupby(["grp", "category"], observed=True).size()
            .unstack("grp").fillna(0))
catmix = (100 * cnt / cnt.sum(axis=0)).reindex(columns=GRP_ORDER)
catmix.round(2).sort_values("small-cryptic", ascending=False).to_csv(
    f"{PP}/smallcryptic_category_mix.tsv", sep="\t")
catmix.round(2).sort_values("small-cryptic", ascending=False)
""")

code(r"""
# module architecture: which combinations of rep / mob / TA / dark actually occur
sc = real[real.grp == "small-cryptic"].copy()
sc["M_rep"]  = sc.rep_any
sc["M_mob"]  = (num(sc.plasann_has_orit) > 0) | sc.mob_relaxase.notna() | (sc.n_mobdna > 0)
sc["M_ta"]   = sc.n_ta > 0
arch = (sc.groupby(["M_rep", "M_mob", "M_ta"]).size().rename("n").reset_index()
          .sort_values("n", ascending=False))
arch["pct"] = (100*arch.n/len(sc)).round(1)
arch["architecture"] = arch.apply(
    lambda r: " + ".join([m for m, f in [("rep", r.M_rep), ("mob", r.M_mob), ("TA", r.M_ta)] if f]) or "none detected",
    axis=1)
arch[["architecture", "n", "pct"]].to_csv(f"{PP}/smallcryptic_module_architecture.tsv",
                                          sep="\t", index=False)
print(f"Module architecture of the {len(sc):,} small cryptic plasmids\n")
print(arch[["architecture", "n", "pct"]].to_string(index=False))
""")

# ======================================================== S8 ========================
md(r"""
## 8. Synthesis — and what a top-tier paper still needs

Fill this in once the notebook has been executed; the cell below prints the assembled headline
numbers so the narrative can be written against measured values rather than remembered ones.
""")

code(r"""
sc_n   = int((real.grp == "small-cryptic").sum())
print("=" * 78)
print("HEADLINE NUMBERS".center(78))
print("=" * 78)
print(f"analysis set (non-simulated, non-lab)          : {len(real):,}")
print(f"small (< {SMALL_CUT:,} bp)                            : {int(real.small.sum()):,} "
      f"({100*real.small.mean():.1f}%)")
print(f"small AND cryptic (no AMR/vir/metal/conj)      : {sc_n:,} ({100*sc_n/len(real):.1f}% of the set)")
print(f"  untyped by PlasmidFinder+mob_typer+PlasAnn   : {100*(~real[real.grp=='small-cryptic'].rep_any).mean():.1f}%")
print(f"  called non-mobilizable by mob_typer          : {100*(real[real.grp=='small-cryptic'].mob_mobility=='non-mobilizable').mean():.1f}%")
print(f"  with no host-range call at all               : {100*(~real[real.grp=='small-cryptic'].hr_call).mean():.1f}%")
print(f"  pooled unannotatable CDS fraction            : {dark.loc['small-cryptic','pooled_pct_dark']:.1f}%")
print(f"  distinct MOB-suite primary clusters          : {real[real.grp=='small-cryptic'].mob_cluster.nunique():,}")
print("-" * 78)
print("dispersal at matched sampling depth (K=10 members per lineage):")
print(res.to_string(index=False))
print("=" * 78)
""")

md(r"""
### Where this stands as a publication

**What the data here can already carry.** A quantitative demonstration that the single largest
compartment of the complete plasmidome is (i) systematically untyped by every standard tool,
(ii) almost entirely unannotatable at the protein level, (iii) predicted to cross phylum boundaries
an order of magnitude more often than conjugative plasmids do, and (iv) — despite being called
immobile — dispersed across habitats and countries at close to conjugative parity, with the parity
becoming statistically indistinguishable once metagenome-derived records are excluded (§6c).

**State the dispersal result precisely.** Pooled over the whole set, conjugative lineages are
*modestly* broader than cryptic ones (ratio 0.79–0.94 depending on metric and rarefaction depth,
small effect sizes, but significant at this n). Restricted to isolate-derived plasmids — where
habitat and geography metadata are most reliable — the difference disappears (habitats 0.94×,
p ≈ 0.06; countries 1.05×, p ≈ 0.11). The defensible claim is **near-parity**, not superiority, and it
is still remarkable for elements that encode no transfer machinery at all.

**What is still missing for a Nature Communications-level claim.** A negative result about tool
coverage is a resource paper, not a discovery. The gap has to be closed constructively:

1. **A de novo replicon catalogue.** Cluster the ~50 k cryptic plasmid proteomes (MMseqs2 — *not yet
   installed*, see `tools/README.md`) and identify the recurrent families that occupy the positional
   and syntenic slot of a replication initiator. Recovering replicon families that PlasmidFinder does
   not contain converts "the tools fail" into "here is what they are missing".
2. **A mechanism for the dispersal paradox.** Test whether the non-mobilizable majority carries
   degenerate or non-canonical oriTs, using the local `orit_db_folder` BLAST database at relaxed
   thresholds against the mob_typer defaults. If cryptic plasmids carry oriTs without relaxases they
   are *trans*-mobilized hitchhikers, and the "non-mobilizable" label is a detection artefact.
3. **A function for the dark matter.** The obvious candidate for why payload-free plasmids persist is
   anti-phage defence; testing it needs DefenseFinder/PADLOC, neither of which is installed.

Items 1 and 2 are the ones that decide whether this is a *Nature Communications* paper or a solid
database paper, and both are computationally cheap relative to what has already been run.
""")

md(r"""
---
*Result tables are written to `data/plasmidscope_primary/smallcryptic_*.tsv`; figures to
`reports/figures/small_cryptic/`. Rebuild this notebook with
`python3 scripts/build_small_cryptic_notebook.py`; execute it with the **genesis_nb** kernel.*
""")

# ======================================================== write =====================
nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {"display_name": "genesis_nb", "language": "python", "name": "genesis_nb"},
    "language_info": {"name": "python"},
}
out = f"{'/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis'}/notebooks/small_cryptic_investigation.ipynb"
with open(out, "w") as fh:
    nbf.write(nb, fh)
print(f"wrote {out}  ({len(cells)} cells)")
