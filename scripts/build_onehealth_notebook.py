#!/usr/bin/env python3
"""Build notebooks/onehealth_amr_investigation.ipynb — the Phase-3 One Health AMR analysis behind
the manuscript `manuscript/amr_onehealth_NAR.md`, as a single runnable notebook.

This consolidates the two scripts that were run for the paper into one notebook, preserving their
exact computations and figures:
  - scripts/analyze_amr_onehealth_mobility.py  (main analysis: 10 tables + 6 figures)
  - scripts/reassess_review_findings.py        (round-1 review reassessment: 6 tables + 1 figure)

The notebook reads data/plasmidscope_primary/plasmid_metadata_master.tsv (+ all_metadata.tsv and
complete_provenance.tsv for the host-taxon join), recomputes every result table under
data/plasmidscope_primary/{onehealth_*,reassess_*}.tsv, and renders the seven figures under
reports/figures/onehealth/ inline. Compartments are built from hab_sub + is_clinical (NOT hab_top);
see reports/environment_taxonomy_LOCKED.md. Method summary: reports/card_amr_methodology.md.

Run: python3 scripts/build_onehealth_notebook.py   (then execute the notebook with the genesis_nb kernel)
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []


def md(t):
    cells.append(nbf.v4.new_markdown_cell(t.strip("\n")))


def code(s):
    cells.append(nbf.v4.new_code_cell(s.strip("\n")))


# ============================================================ intro ==================
md(r"""
# Plasmid-borne AMR across One Health compartments and mobility classes

The Phase-3 analysis behind the manuscript **_Mobility rivals provenance: plasmid-borne antimicrobial
resistance across One Health compartments_** (`manuscript/amr_onehealth_NAR.md`, target venue NAR).
It stratifies the CARD/RGI resistome (`card_*` columns, Perfect+Strict) of the 208,248-plasmid working
set by a **One Health compartment** crossed with **MOB-suite mobility**.

**Compartments are built from `hab_sub` + `is_clinical`, never `hab_top`.** `hab_top`'s
"Host-associated" level collapses human, livestock, plant and insect into one bucket, which destroys
the human / animal / environment structure a One Health analysis needs. `hab_top` is used only to honour
the locked exclusion of the Simulated- and Lab-artifact buckets (`reports/environment_taxonomy_LOCKED.md`),
and that exclusion is expressed here over `hab_sub`. Compartments assign **104,169 plasmids (50.0%)**;
the rest carry no curated habitat and are excluded from every compartment view.

**Headline results.** Resistance prevalence spans **24-fold**, from **44.6%** in the human-clinical
compartment to **1.8%** in wildlife (livestock/poultry second at 38.9%). Mobility is an axis of
comparable magnitude — **conjugative 60.1%** vs **non-mobilizable 6.7%** — and the ordering is monotonic
within every compartment. After adjusting for plasmid size and GC, both axes survive, and some
univariate compartment gaps (plant/food, wildlife) turn out to be size/mobility-composition artefacts.

**Methods.** Wilson 95% CIs for prevalences; χ² with Cramér's *V*; pairwise Fisher vs the natural-
environment baseline with Benjamini–Hochberg FDR; logistic regression of AMR status on
compartment + mobility + log₁₀(size) + GC (references: natural environment / non-mobilizable); Jaccard
over ARO sets; and a clonal-redundancy check redrawing one plasmid per MOB cluster under a fixed seed.

**Provenance.** Built by `scripts/build_onehealth_notebook.py` from
`scripts/analyze_amr_onehealth_mobility.py` (§1–§8) and `scripts/reassess_review_findings.py`
(round-1 review reassessment, R1–R8). CARD v4.0.1 / RGI 6.0.8; see `reports/card_amr_methodology.md`.
""")

# ============================================================ setup ==================
md("## 0 · Setup — data, compartment map, and helpers")

code(r"""
%matplotlib inline
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint

# run from the project root regardless of where the kernel started
while not Path("data/plasmidscope_primary/plasmid_metadata_master.tsv").exists() and Path.cwd() != Path.cwd().parent:
    os.chdir("..")

pd.set_option("display.width", 160)
pd.set_option("display.max_columns", 60)
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 150, "savefig.bbox": "tight",
                     "figure.facecolor": "white", "font.size": 9})

PP = "data/plasmidscope_primary"
MASTER = f"{PP}/plasmid_metadata_master.tsv"
FIGDIR = Path("reports/figures/onehealth"); FIGDIR.mkdir(parents=True, exist_ok=True)
SEED = 20260727
""")

code(r"""
# hab_sub -> (One Health domain, compartment). Human sub-habitats split by is_clinical at runtime;
# anything absent from this map carries no compartment and is excluded.
HUMAN_SUBS = {
    "Human: unspecified", "Human: gut/faeces", "Human: blood", "Human: urine",
    "Human: respiratory", "Human: skin/wound", "Human: other clinical", "Human: oral",
}
COMPARTMENT_MAP = {
    "Livestock: pig": ("Animal", "Animal - livestock & poultry"),
    "Livestock: cattle": ("Animal", "Animal - livestock & poultry"),
    "Livestock: other": ("Animal", "Animal - livestock & poultry"),
    "Poultry/bird": ("Animal", "Animal - livestock & poultry"),
    "Companion animal": ("Animal", "Animal - companion & aquaculture"),
    "Fish/aquaculture": ("Animal", "Animal - companion & aquaculture"),
    "Other mammal": ("Animal", "Animal - wildlife & other"),
    "Rodent": ("Animal", "Animal - wildlife & other"),
    "Insect/arthropod": ("Animal", "Animal - wildlife & other"),
    "Other invertebrate": ("Animal", "Animal - wildlife & other"),
    "Gut/faeces: unspecified host": ("Animal", "Animal - wildlife & other"),
    "Wastewater/sewage": ("Environment", "Environment - engineered interface"),
    "Bioreactor": ("Environment", "Environment - engineered interface"),
    "Built environment": ("Environment", "Environment - engineered interface"),
    "Solid waste/compost": ("Environment", "Environment - engineered interface"),
    "Industrial/remediation": ("Environment", "Environment - engineered interface"),
    "Terrestrial/soil": ("Environment", "Environment - natural"),
    "Aquatic: freshwater/other": ("Environment", "Environment - natural"),
    "Aquatic: marine": ("Environment", "Environment - natural"),
    "Air": ("Environment", "Environment - natural"),
    "Extreme/other": ("Environment", "Environment - natural"),
    "Environmental: unspecified": ("Environment", "Environment - natural"),
    "Plant": ("Plant & food", "Plant & food"),
    "Food/fermentation": ("Plant & food", "Plant & food"),
    "Algae": ("Plant & food", "Plant & food"),
    "Fungi": ("Plant & food", "Plant & food"),
}
# Fixed display / model order; "Environment - natural" is the regression reference level.
ORDER = [
    "Human - clinical", "Human - community",
    "Animal - livestock & poultry", "Animal - companion & aquaculture",
    "Animal - wildlife & other",
    "Environment - engineered interface", "Environment - natural",
    "Plant & food",
]
REF = "Environment - natural"
MOB_ORDER = ["conjugative", "mobilizable", "non-mobilizable"]


def assign_compartment(row):
    sub = row["hab_sub"]
    if sub in HUMAN_SUBS:
        return "Human - clinical" if row["is_clinical"] == 1 else "Human - community"
    return COMPARTMENT_MAP.get(sub, (None, None))[1]


def assign_domain(row):
    sub = row["hab_sub"]
    if sub in HUMAN_SUBS:
        return "Human"
    return COMPARTMENT_MAP.get(sub, (None, None))[0]


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan)
    return proportion_confint(k, n, alpha=0.05, method="wilson")


def split_semi(s):
    "';'-joined list column -> set of non-empty items."
    if not isinstance(s, str) or not s:
        return set()
    return {x.strip() for x in s.split(";") if x.strip()}


def prevalence_table(df, group_col, order=None):
    rows = []
    groups = order if order is not None else sorted(df[group_col].dropna().unique())
    for g in groups:
        sub = df[df[group_col] == g]
        n = len(sub)
        if n == 0:
            continue
        k = int(sub.amr_pos.sum())
        km = int(sub.card_multidrug.sum())
        lo, hi = wilson(k, n)
        mlo, mhi = wilson(km, n)
        carriers = sub.loc[sub.amr_pos, "card_n_arg"]
        rows.append({
            group_col: g, "n": n,
            "amr_pos": k, "amr_prevalence": k / n, "amr_ci_lo": lo, "amr_ci_hi": hi,
            "multidrug": km, "multidrug_prevalence": km / n, "md_ci_lo": mlo, "md_ci_hi": mhi,
            "median_args_in_carriers": float(carriers.median()) if k else 0.0,
            "mean_args_in_carriers": float(carriers.mean()) if k else 0.0,
            "max_args": int(sub.card_n_arg.max()),
        })
    return pd.DataFrame(rows)


def cramers_v(ct):
    chi2 = stats.chi2_contingency(ct)[0]
    n = ct.values.sum()
    return np.sqrt(chi2 / (n * (min(ct.shape) - 1)))
""")

code(r"""
df = pd.read_csv(MASTER, sep="\t", low_memory=False)
df["amr_pos"] = df.card_n_arg > 0
df["compartment"] = df.apply(assign_compartment, axis=1)
df["oh_domain"] = df.apply(assign_domain, axis=1)

oh = df[df.compartment.notna()].copy()
oh["compartment"] = pd.Categorical(oh.compartment, categories=ORDER, ordered=True)
print(f"master:              {len(df):,} plasmids x {df.shape[1]} columns")
print(f"One Health set:      {len(oh):,} ({100*len(oh)/len(df):.1f}% of working set)")
print(f"excluded (no comp.): {len(df)-len(oh):,}")
""")

# ============================================================ §1 =====================
md(r"""
## 1 · A 24-fold resistance gradient across One Health compartments

AMR / multidrug prevalence per compartment with Wilson 95% CIs, a global χ² with Cramér's *V*, and
pairwise Fisher tests against the natural-environment baseline (Benjamini–Hochberg FDR).
→ `onehealth_prevalence_by_compartment.tsv`, `onehealth_pairwise_vs_natural.tsv`. **(Paper Table 1.)**
""")

code(r"""
comp_tab = prevalence_table(oh, "compartment", ORDER)
comp_tab.to_csv(f"{PP}/onehealth_prevalence_by_compartment.tsv", sep="\t", index=False)

ct = pd.crosstab(oh.compartment, oh.amr_pos)
chi2, p, dof, _ = stats.chi2_contingency(ct)
print(f"chi2={chi2:,.1f}  dof={dof}  p={p:.3e}  Cramers_V={cramers_v(ct):.3f}")
comp_tab[["compartment", "n", "amr_pos", "amr_prevalence", "amr_ci_lo", "amr_ci_hi",
          "multidrug_prevalence", "median_args_in_carriers"]].round(4)
""")

code(r"""
# pairwise vs the natural-environment baseline, BH-corrected
ref_sub = oh[oh.compartment == REF]
rows = []
for g in ORDER:
    if g == REF:
        continue
    s = oh[oh.compartment == g]
    table = [[int(s.amr_pos.sum()), len(s) - int(s.amr_pos.sum())],
             [int(ref_sub.amr_pos.sum()), len(ref_sub) - int(ref_sub.amr_pos.sum())]]
    orv, pv = stats.fisher_exact(table)
    rows.append({"compartment": g, "vs": REF, "odds_ratio": orv, "p_raw": pv})
pw = pd.DataFrame(rows)
pw["p_bh"] = multipletests(pw.p_raw, method="fdr_bh")[1]
pw.to_csv(f"{PP}/onehealth_pairwise_vs_natural.tsv", sep="\t", index=False)
pw.round(4)
""")

code(r"""
# Figure 1 — prevalence by compartment (Wilson 95% CI)
colors = plt.cm.tab10(np.linspace(0, 1, 10))
fig, ax = plt.subplots(figsize=(7.2, 4))
y = np.arange(len(comp_tab))
ax.barh(y, comp_tab.amr_prevalence * 100, color=colors[:len(comp_tab)],
        xerr=[(comp_tab.amr_prevalence - comp_tab.amr_ci_lo) * 100,
              (comp_tab.amr_ci_hi - comp_tab.amr_prevalence) * 100],
        error_kw={"ecolor": "0.3", "lw": 1})
ax.set_yticks(y)
ax.set_yticklabels([f"{c}\n(n={n:,})" for c, n in zip(comp_tab.compartment, comp_tab.n)])
ax.invert_yaxis()
ax.set_xlabel("plasmids carrying >=1 confident ARG (%)")
ax.set_title("Plasmid-borne AMR prevalence by One Health compartment\n"
             "(CARD v4.0.1, Perfect+Strict; Wilson 95% CI)")
fig.savefig(f"{FIGDIR}/amr_prevalence_by_compartment.png")
plt.show()
""")

# ============================================================ §2 =====================
md(r"""
## 2 · Mobility is an axis of comparable magnitude

AMR / multidrug prevalence by MOB-suite mobility class. → `onehealth_prevalence_by_mobility.tsv`.
**(Paper Table 2.)**
""")

code(r"""
mob = oh[oh.mob_mobility.notna()].copy()
mob_tab = prevalence_table(mob, "mob_mobility", MOB_ORDER)
mob_tab.to_csv(f"{PP}/onehealth_prevalence_by_mobility.tsv", sep="\t", index=False)
mob_tab[["mob_mobility", "n", "amr_pos", "amr_prevalence", "amr_ci_lo", "amr_ci_hi",
         "multidrug_prevalence"]].round(4)
""")

# ============================================================ §3 =====================
md(r"""
## 3 · The mobility ordering holds inside every compartment

Compartment × mobility joint prevalence — a Simpson's-paradox check. The conjugative > mobilizable >
non-mobilizable ordering is monotonic in every compartment. → `onehealth_compartment_x_mobility.tsv`.
**(Paper Figure 2.)**
""")

code(r"""
joint = (mob.groupby(["compartment", "mob_mobility"], observed=True)
            .agg(n=("amr_pos", "size"), amr_pos=("amr_pos", "sum")).reset_index())
joint["amr_prevalence"] = joint.amr_pos / joint.n
joint.to_csv(f"{PP}/onehealth_compartment_x_mobility.tsv", sep="\t", index=False)
pivot = joint.pivot(index="compartment", columns="mob_mobility",
                    values="amr_prevalence").reindex(index=ORDER, columns=MOB_ORDER)
npivot = joint.pivot(index="compartment", columns="mob_mobility",
                     values="n").reindex(index=ORDER, columns=MOB_ORDER)
(pivot * 100).round(1)
""")

code(r"""
# Figure 2 — compartment x mobility heatmap
fig, ax = plt.subplots(figsize=(5.6, 4.2))
im = ax.imshow(pivot.values * 100, cmap="YlOrRd", aspect="auto")
ax.set_xticks(range(len(MOB_ORDER))); ax.set_xticklabels(MOB_ORDER, rotation=20, ha="right")
ax.set_yticks(range(len(ORDER))); ax.set_yticklabels(ORDER)
for i in range(len(ORDER)):
    for j in range(len(MOB_ORDER)):
        val = pivot.values[i, j]
        if not np.isnan(val):
            ax.text(j, i, f"{val*100:.0f}%\nn={int(npivot.values[i,j]):,}", ha="center",
                    va="center", fontsize=6.5, color="white" if val > 0.45 else "black")
fig.colorbar(im, ax=ax, label="AMR prevalence (%)")
ax.set_title("AMR prevalence: compartment x mobility")
fig.savefig(f"{FIGDIR}/amr_compartment_x_mobility.png")
plt.show()
""")

# ============================================================ §4 =====================
md(r"""
## 4 · Adjustment separates genuine compartment effects from composition

Logistic regression of AMR status on compartment + mobility + log₁₀(size) + GC, references natural
environment / non-mobilizable. Odds ratios with 95% CIs. → `onehealth_logit_odds_ratios.tsv`.
**(Paper Figure 3.)**
""")

code(r"""
reg = mob[mob.size_bp > 0].copy()
reg["amr_pos"] = reg.amr_pos.astype(int)   # statsmodels reads bool endog as categorical
reg["log_size"] = np.log10(reg.size_bp)
reg["comp"] = reg.compartment.astype(str)
reg["mobc"] = reg.mob_mobility.astype(str)
model = smf.logit(
    f'amr_pos ~ C(comp, Treatment(reference="{REF}")) '
    f'+ C(mobc, Treatment(reference="non-mobilizable")) + log_size + gc_percent',
    data=reg).fit(disp=0)
res = pd.DataFrame({"term": model.params.index, "coef": model.params.values, "p": model.pvalues.values})
ci = model.conf_int()
res["or"] = np.exp(res.coef); res["or_lo"] = np.exp(ci[0].values); res["or_hi"] = np.exp(ci[1].values)
res.to_csv(f"{PP}/onehealth_logit_odds_ratios.tsv", sep="\t", index=False)
print(f"n={len(reg):,}  pseudo-R2={model.prsquared:.3f}")
res[["term", "or", "or_lo", "or_hi", "p"]].round(4)
""")

code(r"""
# Figure 3 — adjusted odds-ratio forest plot
plot_terms = res[res.term.str.contains(r"C\(comp|C\(mobc")].copy()
plot_terms["label"] = (plot_terms.term.str.replace(r".*\[T\.", "", regex=True)
                                 .str.replace(r"\]$", "", regex=True))
fig, ax = plt.subplots(figsize=(6.4, 4))
yy = np.arange(len(plot_terms))
ax.errorbar(plot_terms["or"], yy,
            xerr=[plot_terms["or"] - plot_terms.or_lo, plot_terms.or_hi - plot_terms["or"]],
            fmt="o", color="#22456b", ecolor="0.4", capsize=3)
ax.axvline(1, color="crimson", ls="--", lw=1)
ax.set_yticks(yy); ax.set_yticklabels(plot_terms.label, fontsize=8)
ax.invert_yaxis(); ax.set_xscale("log")
ax.set_xlabel("adjusted odds ratio for carrying >=1 ARG (log scale)")
ax.set_title(f"AMR odds, adjusted for plasmid size and GC\n(ref: {REF} / non-mobilizable)")
fig.savefig(f"{FIGDIR}/amr_adjusted_odds_ratios.png")
plt.show()
""")

# ============================================================ §5 =====================
md(r"""
## 5 · Compartment-specific drug-class signatures

Drug-class spectrum of the resistome per compartment (fraction of that compartment's AMR+ carriers
carrying each class), top-15 classes. → `onehealth_drug_class_by_compartment.tsv`. **(Paper Figure 4.)**
""")

code(r"""
pos = oh[oh.amr_pos].copy()
dc_rows = []
for g in ORDER:
    s = pos[pos.compartment == g]
    if len(s) == 0:
        continue
    counter = {}
    for v in s.card_drug_classes:
        for c in split_semi(v):
            counter[c] = counter.get(c, 0) + 1
    for c, k in counter.items():
        dc_rows.append({"compartment": g, "drug_class": c, "n_plasmids": k, "frac_of_carriers": k / len(s)})
dc = pd.DataFrame(dc_rows)
dc.to_csv(f"{PP}/onehealth_drug_class_by_compartment.tsv", sep="\t", index=False)
top_classes = (dc.groupby("drug_class").n_plasmids.sum().sort_values(ascending=False).head(15).index.tolist())
dc_pivot = (dc[dc.drug_class.isin(top_classes)]
            .pivot(index="drug_class", columns="compartment", values="frac_of_carriers")
            .reindex(index=top_classes, columns=ORDER))
(dc_pivot * 100).round(1)
""")

code(r"""
# Figure 4 — drug-class spectrum heatmap
fig, ax = plt.subplots(figsize=(7.6, 5))
im = ax.imshow(dc_pivot.values, cmap="viridis", aspect="auto")
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=7.5)
ax.set_yticks(range(len(dc_pivot))); ax.set_yticklabels(dc_pivot.index, fontsize=7.5)
fig.colorbar(im, ax=ax, label="fraction of AMR+ plasmids in compartment")
ax.set_title("Drug-class spectrum of the plasmid resistome, by compartment")
fig.savefig(f"{FIGDIR}/drug_class_by_compartment.png")
plt.show()
""")

# ============================================================ §6 =====================
md(r"""
## 6 · Resistance repertoires form a human block and an agricultural–environmental block

Cross-compartment sharing of ARG repertoires: Jaccard index over the sets of distinct ARO accessions,
plus unique-ARO richness per compartment. → `onehealth_arg_jaccard.tsv`, `onehealth_arg_richness.tsv`.
**(Paper Figure 5.)**
""")

code(r"""
aro_sets = {}
for g in ORDER:
    s = pos[pos.compartment == g]
    acc = set()
    for v in s.card_aro_list:
        acc |= split_semi(v)
    aro_sets[g] = acc
jac = pd.DataFrame(index=ORDER, columns=ORDER, dtype=float)
for a in ORDER:
    for b in ORDER:
        ua, ub = aro_sets[a], aro_sets[b]
        jac.loc[a, b] = len(ua & ub) / len(ua | ub) if (ua | ub) else np.nan
jac.to_csv(f"{PP}/onehealth_arg_jaccard.tsv", sep="\t")
rich = pd.DataFrame({"compartment": ORDER,
                     "n_carriers": [int((pos.compartment == g).sum()) for g in ORDER],
                     "n_unique_aro": [len(aro_sets[g]) for g in ORDER]})
rich.to_csv(f"{PP}/onehealth_arg_richness.tsv", sep="\t", index=False)
display(rich)
jac.round(3)
""")

code(r"""
# Figure 5 — ARG-repertoire Jaccard heatmap
fig, ax = plt.subplots(figsize=(6.2, 5))
im = ax.imshow(jac.values.astype(float), cmap="magma", vmin=0, vmax=1)
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=7.5)
ax.set_yticks(range(len(ORDER))); ax.set_yticklabels(ORDER, fontsize=7.5)
for i in range(len(ORDER)):
    for j in range(len(ORDER)):
        v_ = jac.values[i, j]
        if not np.isnan(v_):
            ax.text(j, i, f"{v_:.2f}", ha="center", va="center", fontsize=6.5,
                    color="white" if v_ < 0.6 else "black")
fig.colorbar(im, ax=ax, label="Jaccard index (shared AROs)")
ax.set_title("Cross-compartment sharing of ARG repertoires")
fig.savefig(f"{FIGDIR}/arg_repertoire_jaccard.png")
plt.show()
""")

# ============================================================ §7 =====================
md(r"""
## 7 · Multireplicon architecture and resistance outside the clinic

Within each compartment, does a multireplicon plasmid (PlasmidFinder ≥2 distinct Inc) carry AMR more
often than a single/none plasmid? Fisher per compartment, BH-corrected. → `onehealth_multireplicon_amr.tsv`.
""")

code(r"""
mr = oh[oh.pf_n_inc.notna()].copy()
mr["multireplicon"] = mr.pf_n_inc >= 2
mr_rows = []
for g in ORDER:
    s = mr[mr.compartment == g]
    a = s[s.multireplicon]; b = s[~s.multireplicon]
    if len(a) == 0 or len(b) == 0:
        continue
    orv, pv = stats.fisher_exact(
        [[int(a.amr_pos.sum()), len(a) - int(a.amr_pos.sum())],
         [int(b.amr_pos.sum()), len(b) - int(b.amr_pos.sum())]])
    mr_rows.append({"compartment": g, "n_multireplicon": len(a),
                    "amr_prev_multireplicon": a.amr_pos.mean(), "n_single_or_none": len(b),
                    "amr_prev_other": b.amr_pos.mean(), "odds_ratio": orv, "p_raw": pv})
mrt = pd.DataFrame(mr_rows)
mrt["p_bh"] = multipletests(mrt.p_raw, method="fdr_bh")[1]
mrt.to_csv(f"{PP}/onehealth_multireplicon_amr.tsv", sep="\t", index=False)
mrt.round(4)
""")

# ============================================================ §8 =====================
md(r"""
## 8 · Clonal redundancy inflates absolute prevalence but preserves the gradient

Robustness check: redraw one plasmid at random per MOB cluster (fixed seed) and recompute prevalence.
Absolute prevalence roughly halves, but the compartment ordering is intact.
→ `onehealth_dereplication_sensitivity.tsv`. **(Paper Figure 6.)**
""")

code(r"""
dedup = (oh[oh.mob_cluster.notna()].groupby("mob_cluster", observed=True).sample(n=1, random_state=SEED))
ded_tab = prevalence_table(dedup, "compartment", ORDER)[["compartment", "n", "amr_prevalence"]].rename(
    columns={"n": "n_dedup", "amr_prevalence": "amr_prevalence_dedup"})
comp_cmp = comp_tab[["compartment", "n", "amr_prevalence"]].merge(ded_tab, on="compartment", how="left")
comp_cmp.to_csv(f"{PP}/onehealth_dereplication_sensitivity.tsv", sep="\t", index=False)
print(f"one plasmid per MOB cluster: n={len(dedup):,}")
comp_cmp.round(4)
""")

code(r"""
# Figure 6 — dereplication sensitivity
fig, ax = plt.subplots(figsize=(7, 3.8))
xx = np.arange(len(comp_cmp))
ax.bar(xx - 0.2, comp_cmp.amr_prevalence * 100, 0.4, label="all plasmids")
ax.bar(xx + 0.2, comp_cmp.amr_prevalence_dedup * 100, 0.4, label="1 per MOB cluster")
ax.set_xticks(xx); ax.set_xticklabels(comp_cmp.compartment, rotation=35, ha="right", fontsize=7.5)
ax.set_ylabel("AMR prevalence (%)"); ax.legend(frameon=False)
ax.set_title("Clonal-redundancy robustness check")
fig.savefig(f"{FIGDIR}/dereplication_sensitivity.png")
plt.show()
""")

# ============================================ round-1 review reassessment ============
md(r"""
---
# Round-1 peer-review reassessment (R1–R8)

The computable items of the round-1 revision roadmap (`manuscript/review_round1.md`), reproduced from
`scripts/reassess_review_findings.py`. Host taxonomy is joined from PlasmidScope `all_metadata.tsv`
via the provenance representative key (the join used by `scripts/build_analysis_table.py`).

- **R1** host-taxon coverage per compartment (can the confounder even be tested?)
- **R2** collinearity — VIF + mobility-OR stability with/without log size
- **R3** cluster-robust SEs (groups = MOB cluster)
- **R4** compartment ORs before/after host-taxon adjustment
- **R5** adjusted headline contrast at matched size
- **R6** assignment bias — AMR in assigned vs unassigned plasmids
- **R7** PlasmidFinder typing coverage per compartment
- **R8** what the `is_clinical` flag rests on

→ `data/plasmidscope_primary/reassess_*.tsv`, `reports/figures/onehealth/reassess_compartment_or_shift.png`.
""")

code(r"""
import csv, sys
csv.field_size_limit(sys.maxsize)

# Enterobacteriaceae (sensu lato) genera used for the detection-bias argument.
ENTERO = {
    "Escherichia", "Klebsiella", "Salmonella", "Enterobacter", "Citrobacter", "Serratia",
    "Proteus", "Shigella", "Cronobacter", "Providencia", "Morganella", "Yersinia",
    "Kluyvera", "Raoultella", "Leclercia", "Edwardsiella", "Hafnia", "Pluralibacter",
    "Lelliottia", "Atlantibacter", "Phytobacter", "Kosakonia",
}
NULL_HOST = {"-", "", "synthetic construct", "unclassified", "uncultured bacterium"}

# join PlasmidScope host taxonomy via the provenance representative key
rep2full = {}
with open(f"{PP}/complete_provenance.tsv", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        rep2full[r["plasmid_id"].split(",")[0]] = r["plasmid_id"]
host = {}
with open(f"{PP}/all_metadata.tsv", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        host[r["Plasmid_ID"]] = r.get("Host", "")
df["host_raw"] = df.plasmid_id.map(lambda p: host.get(rep2full.get(p, ""), "")).fillna("")
df["has_host"] = ~df.host_raw.str.strip().isin(NULL_HOST)
df["host_genus"] = np.where(df.has_host, df.host_raw.str.split().str[0], None)
df["is_entero"] = df.host_genus.isin(ENTERO).astype(int)

ohr = df[df.compartment.notna()].copy()

def or_table(m, label):
    ci = m.conf_int()
    out = pd.DataFrame({"term": m.params.index, "or": np.exp(m.params.values),
                        "or_lo": np.exp(ci[0].values), "or_hi": np.exp(ci[1].values), "p": m.pvalues.values})
    out["model"] = label
    out["k"] = out.term.str.replace(r".*\[T\.", "", regex=True).str.replace(r"\]$", "", regex=True)
    return out

print(f"working set {len(df):,}; One Health set {len(ohr):,}; with host taxonomy {int(df.has_host.sum()):,}")
""")

md("### R1 · Host-taxonomy coverage by compartment → `reassess_host_coverage.tsv`")
code(r"""
cov = (ohr.groupby("compartment").agg(n=("amr_pos", "size"), with_host=("has_host", "sum"),
                                      entero_frac=("is_entero", "mean")).reindex(ORDER))
cov["host_coverage"] = cov.with_host / cov.n
cov["entero_frac_of_hosted"] = ohr[ohr.has_host].groupby("compartment").is_entero.mean().reindex(ORDER)
cov.to_csv(f"{PP}/reassess_host_coverage.tsv", sep="\t")
cov[["n", "with_host", "host_coverage", "entero_frac_of_hosted"]].round(3)
""")

md("### R6 · Assignment bias — AMR prevalence, assigned vs unassigned → `reassess_assignment_bias.tsv`")
code(r"""
df["assigned"] = df.compartment.notna()
df["unassigned_kind"] = np.where(
    df.assigned, "assigned",
    np.where(df.hab_sub == "Simulated communities", "simulated",
             np.where(df.hab_sub.isin(["Unlabelled", "Geography only (no habitat)", "Unresolved"]),
                      "no label", "other excluded")))
ak = df.groupby("unassigned_kind").agg(n=("amr_pos", "size"), amr=("amr_pos", "sum"))
ak["prevalence"] = ak.amr / ak.n
ak.to_csv(f"{PP}/reassess_assignment_bias.tsv", sep="\t")
ak.round(4)
""")

md("### R7 · PlasmidFinder typing coverage by compartment → `reassess_typing_coverage.tsv`")
code(r"""
ohr["typed"] = (ohr.pf_n_inc.fillna(0) >= 1).astype(int)
tc = (ohr.groupby("compartment").agg(n=("typed", "size"), typed=("typed", "sum")).reindex(ORDER))
tc["typing_coverage"] = tc.typed / tc.n
tc.to_csv(f"{PP}/reassess_typing_coverage.tsv", sep="\t")
tc.round(3)
""")

md("### R8 · What the `is_clinical` flag rests on → `reassess_clinical_basis.tsv`")
code(r"""
clin = ohr[ohr.compartment == "Human - clinical"]
body_site = clin.hab_sub.isin(["Human: blood", "Human: urine", "Human: respiratory",
                               "Human: skin/wound", "Human: other clinical"])
has_disease = clin.plsdb_disease_tags.notna() & (clin.plsdb_disease_tags != "")
basis = pd.DataFrame({
    "basis": ["clinical body site", "no body site, has PLSDB disease tag",
              "no body site, no disease tag (context only)"],
    "n": [int(body_site.sum()), int((~body_site & has_disease).sum()),
          int((~body_site & ~has_disease).sum())]})
basis["frac"] = basis.n / len(clin)
basis.to_csv(f"{PP}/reassess_clinical_basis.tsv", sep="\t", index=False)
basis.round(3)
""")

md("### Modelling set + R2 collinearity (VIF, mobility-OR stability) → `reassess_vif.tsv`")
code(r"""
reg2 = ohr[ohr.mob_mobility.notna() & (ohr.size_bp > 0)].copy()
reg2["amr_pos"] = reg2.amr_pos.astype(int)   # statsmodels reads bool endog as categorical
reg2["log_size"] = np.log10(reg2.size_bp)
reg2["comp"] = reg2.compartment.astype(str)
reg2["mobc"] = reg2.mob_mobility.astype(str)
reg2["cluster"] = reg2.mob_cluster.fillna("NA")
BASE = (f'amr_pos ~ C(comp, Treatment(reference="{REF}")) '
        f'+ C(mobc, Treatment(reference="non-mobilizable")) + log_size + gc_percent')

from statsmodels.stats.outliers_influence import variance_inflation_factor
dm = pd.get_dummies(reg2[["comp", "mobc"]], drop_first=True).astype(float)
dm["log_size"] = reg2.log_size.values; dm["gc_percent"] = reg2.gc_percent.values
dm.insert(0, "const", 1.0)
vif = pd.DataFrame({"term": dm.columns,
                    "VIF": [variance_inflation_factor(dm.values, i) for i in range(dm.shape[1])]})
vif = vif[vif.term != "const"]
vif.to_csv(f"{PP}/reassess_vif.tsv", sep="\t", index=False)

m_full = smf.logit(BASE, data=reg2).fit(disp=0)
m_nosize = smf.logit(BASE.replace(" + log_size", ""), data=reg2).fit(disp=0)
print("mobility OR with vs without log_size in the model:")
for k in ["conjugative", "mobilizable"]:
    a = or_table(m_full, "with").set_index("k").loc[k]
    b = or_table(m_nosize, "without").set_index("k").loc[k]
    print(f"  {k:16s} with size OR={a['or']:.2f}   without size OR={b['or']:.2f}")
vif.sort_values("VIF", ascending=False).round(2)
""")

md("### R3 · Cluster-robust standard errors (groups = MOB cluster) → `reassess_cluster_robust.tsv`")
code(r"""
m_clu = smf.logit(BASE, data=reg2).fit(disp=0, cov_type="cluster", cov_kwds={"groups": reg2["cluster"]})
a = or_table(m_full, "nominal").set_index("k"); b = or_table(m_clu, "clustered").set_index("k")
cmp_rows = []
for k in a.index:
    if k == "Intercept":
        continue
    cmp_rows.append({
        "term": k, "or": round(a.loc[k, "or"], 3),
        "nominal_lo": round(a.loc[k, "or_lo"], 3), "nominal_hi": round(a.loc[k, "or_hi"], 3),
        "clustered_lo": round(b.loc[k, "or_lo"], 3), "clustered_hi": round(b.loc[k, "or_hi"], 3),
        "clustered_p": b.loc[k, "p"],
        "ci_width_ratio": round((b.loc[k, "or_hi"] - b.loc[k, "or_lo"]) / (a.loc[k, "or_hi"] - a.loc[k, "or_lo"]), 1),
        "still_sig_0.05": bool(b.loc[k, "p"] < 0.05)})
clu = pd.DataFrame(cmp_rows)
clu.to_csv(f"{PP}/reassess_cluster_robust.tsv", sep="\t", index=False)
print(f"n clusters = {reg2['cluster'].nunique():,}")
clu
""")

md("### R4 · Compartment ORs before vs after host-taxon adjustment → `reassess_host_adjustment.tsv`")
code(r"""
hs = reg2[reg2.has_host].copy()
top_genera = hs.host_genus.value_counts().head(20).index.tolist()
hs["genus_grp"] = np.where(hs.host_genus.isin(top_genera), hs.host_genus, "other")
m_h0 = smf.logit(BASE, data=hs).fit(disp=0)                                            # no host
m_h1 = smf.logit(BASE + " + is_entero", data=hs).fit(disp=0)                           # + Entero indicator
m_h2 = smf.logit(BASE + ' + C(genus_grp, Treatment(reference="other"))', data=hs).fit(disp=0)  # + top-20 genus
t0 = or_table(m_h0, "no host").set_index("k"); t1 = or_table(m_h1, "+entero").set_index("k")
t2 = or_table(m_h2, "+genus").set_index("k")
rows = []
for k in [c for c in ORDER if c != REF] + ["conjugative", "mobilizable", "log_size"]:
    if k not in t0.index:
        continue
    rows.append({"term": k, "OR_no_host": round(t0.loc[k, "or"], 2),
                 "OR_plus_entero": round(t1.loc[k, "or"], 2), "OR_plus_genus": round(t2.loc[k, "or"], 2),
                 "pct_change_genus": round(100 * (t2.loc[k, "or"] - t0.loc[k, "or"]) / t0.loc[k, "or"], 1),
                 "p_genus_model": t2.loc[k, "p"]})
host_cmp = pd.DataFrame(rows)
host_cmp.to_csv(f"{PP}/reassess_host_adjustment.tsv", sep="\t", index=False)
print(f"subset with host taxonomy: n={len(hs):,}   pseudo-R2: no host {m_h0.prsquared:.3f} | "
      f"+entero {m_h1.prsquared:.3f} | +genus {m_h2.prsquared:.3f}")
host_cmp
""")

code(r"""
# Figure (supplementary) — compartment OR shift under host adjustment
plot = host_cmp[host_cmp.term.isin([c for c in ORDER if c != REF])]
fig, ax = plt.subplots(figsize=(7, 4))
y = np.arange(len(plot))
ax.barh(y - 0.2, plot.OR_no_host, 0.4, label="compartment only")
ax.barh(y + 0.2, plot.OR_plus_genus, 0.4, label="+ host genus")
ax.axvline(1, color="crimson", ls="--", lw=1)
ax.set_yticks(y); ax.set_yticklabels(plot.term, fontsize=8)
ax.invert_yaxis(); ax.set_xscale("log")
ax.set_xlabel("adjusted odds ratio for carrying >=1 ARG (log scale)")
ax.set_title("Compartment effects before and after host-taxon adjustment")
ax.legend(frameon=False)
fig.savefig(f"{FIGDIR}/reassess_compartment_or_shift.png")
plt.show()
""")

md("### R5 · Adjusted headline contrast at matched plasmid size → `reassess_headline_contrast.tsv`")
code(r"""
med_gc = float(reg2.gc_percent.median())
scen = []
for size_bp in [5000, 10000, 50000, 100000]:
    for comp_, mob_ in [("Environment - natural", "conjugative"),
                        ("Human - clinical", "non-mobilizable"),
                        ("Human - clinical", "conjugative"),
                        ("Environment - natural", "non-mobilizable")]:
        scen.append({"comp": comp_, "mobc": mob_, "log_size": np.log10(size_bp),
                     "gc_percent": med_gc, "size_bp": size_bp})
sc = pd.DataFrame(scen)
sc["pred_prob"] = m_full.predict(sc)
sc.to_csv(f"{PP}/reassess_headline_contrast.tsv", sep="\t", index=False)
piv = sc.pivot_table(index="size_bp", columns=["comp", "mobc"], values="pred_prob")

cell_a = reg2[(reg2.comp == "Environment - natural") & (reg2.mobc == "conjugative")]
cell_b = reg2[(reg2.comp == "Human - clinical") & (reg2.mobc == "non-mobilizable")]
print(f"env-natural/conjugative : n={len(cell_a):,} median={cell_a.size_bp.median():,.0f} bp  raw prev={cell_a.amr_pos.mean():.3f}")
print(f"clinical/non-mobilizable: n={len(cell_b):,} median={cell_b.size_bp.median():,.0f} bp  raw prev={cell_b.amr_pos.mean():.3f}")
print("\npredicted P(AMR+) at matched size, median GC (%):")
(piv * 100).round(1)
""")

md(r"""
---
*All result tables are written under `data/plasmidscope_primary/` (`onehealth_*.tsv`, `reassess_*.tsv`)
and the seven figures under `reports/figures/onehealth/`. Rebuild this notebook with
`python3 scripts/build_onehealth_notebook.py`; execute it with the **genesis_nb** kernel.*
""")

# ============================================================ write =================
nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {"display_name": "genesis_nb", "language": "python", "name": "genesis_nb"},
    "language_info": {"name": "python"},
}
out = "notebooks/onehealth_amr_investigation.ipynb"
with open(out, "w") as fh:
    nbf.write(nb, fh)
print(f"wrote {out}  ({len(cells)} cells)")
