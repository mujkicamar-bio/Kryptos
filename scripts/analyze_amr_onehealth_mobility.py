#!/usr/bin/env python3
"""Phase-3 analysis: plasmid-borne AMR across One Health compartments and mobility classes.

Stratifies the CARD/RGI resistome (`card_*` columns) of the 208,248-plasmid working set by a
**One Health compartment** built from `hab_sub` + `is_clinical`, crossed with MOB-suite mobility
(`mob_mobility`). `hab_top` is deliberately NOT used as a stratifier: its "Host-associated" level
collapses human, livestock, plant and insect into one bucket, which destroys the human/animal/
environment structure a One Health analysis depends on. `hab_top` values are used only to honour the
locked exclusion of the Simulated-artifact and Lab-artifact buckets (see
`reports/environment_taxonomy_LOCKED.md`), and that exclusion is expressed here over `hab_sub`.

Analyses:
  1. AMR / multidrug prevalence per compartment, with Wilson 95% CIs
  2. AMR prevalence per mobility class
  3. Compartment x mobility joint prevalence (Simpson's-paradox check)
  4. Logistic regression AMR+ ~ compartment + mobility + log10(size) + GC (odds ratios)
  5. Drug-class spectrum per compartment
  6. Cross-compartment ARG-repertoire sharing (Jaccard over ARO sets)
  7. Multireplicon (PlasmidFinder >=2 Inc) vs AMR, per compartment
  8. Robustness: prevalence recomputed on one plasmid per MOB cluster (clonal de-duplication)

Inputs:  data/plasmidscope_primary/plasmid_metadata_master.tsv
Outputs: data/plasmidscope_primary/onehealth_*.tsv  (10 result tables)
         reports/figures/onehealth/*.png            (6 figures)
Run:     python3 scripts/analyze_amr_onehealth_mobility.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint

PP = "data/plasmidscope_primary"
MASTER = f"{PP}/plasmid_metadata_master.tsv"
FIGDIR = "reports/figures/onehealth"
SEED = 20260727

# hab_sub -> (One Health domain, compartment). Human sub-habitats are split by `is_clinical`
# at runtime. Anything absent from this map carries no compartment and is excluded.
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
    """';'-joined list column -> set of non-empty items."""
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
            "multidrug": km, "multidrug_prevalence": km / n,
            "md_ci_lo": mlo, "md_ci_hi": mhi,
            "median_args_in_carriers": float(carriers.median()) if k else 0.0,
            "mean_args_in_carriers": float(carriers.mean()) if k else 0.0,
            "max_args": int(sub.card_n_arg.max()),
        })
    return pd.DataFrame(rows)


def cramers_v(ct):
    chi2 = stats.chi2_contingency(ct)[0]
    n = ct.values.sum()
    return np.sqrt(chi2 / (n * (min(ct.shape) - 1)))


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    plt.rcParams.update({"figure.dpi": 150, "savefig.bbox": "tight", "font.size": 9})

    df = pd.read_csv(MASTER, sep="\t", low_memory=False)
    print(f"master: {len(df):,} plasmids x {df.shape[1]} columns")

    df["amr_pos"] = df.card_n_arg > 0
    df["compartment"] = df.apply(assign_compartment, axis=1)
    df["oh_domain"] = df.apply(assign_domain, axis=1)

    oh = df[df.compartment.notna()].copy()
    oh["compartment"] = pd.Categorical(oh.compartment, categories=ORDER, ordered=True)
    print(f"One Health analysis set: {len(oh):,} "
          f"({100*len(oh)/len(df):.1f}% of working set); "
          f"excluded {len(df)-len(oh):,} without a compartment")

    # --- 1. prevalence by compartment -------------------------------------------------
    comp_tab = prevalence_table(oh, "compartment", ORDER)
    comp_tab.to_csv(f"{PP}/onehealth_prevalence_by_compartment.tsv", sep="\t", index=False)
    print("\n== AMR prevalence by One Health compartment ==")
    print(comp_tab[["compartment", "n", "amr_pos", "amr_prevalence",
                    "multidrug_prevalence"]].to_string(index=False))

    ct = pd.crosstab(oh.compartment, oh.amr_pos)
    chi2, p, dof, _ = stats.chi2_contingency(ct)
    v = cramers_v(ct)
    print(f"chi2={chi2:.1f} dof={dof} p={p:.3e} Cramers_V={v:.3f}")

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
    print("\n== Fisher vs Environment-natural (BH-corrected) ==")
    print(pw.to_string(index=False))

    # --- 2. prevalence by mobility ----------------------------------------------------
    mob = oh[oh.mob_mobility.notna()].copy()
    mob_tab = prevalence_table(mob, "mob_mobility", MOB_ORDER)
    mob_tab.to_csv(f"{PP}/onehealth_prevalence_by_mobility.tsv", sep="\t", index=False)
    print("\n== AMR prevalence by mobility ==")
    print(mob_tab[["mob_mobility", "n", "amr_prevalence", "multidrug_prevalence"]]
          .to_string(index=False))

    # --- 3. compartment x mobility ----------------------------------------------------
    joint = (mob.groupby(["compartment", "mob_mobility"], observed=True)
                .agg(n=("amr_pos", "size"), amr_pos=("amr_pos", "sum")).reset_index())
    joint["amr_prevalence"] = joint.amr_pos / joint.n
    joint.to_csv(f"{PP}/onehealth_compartment_x_mobility.tsv", sep="\t", index=False)
    pivot = joint.pivot(index="compartment", columns="mob_mobility",
                        values="amr_prevalence").reindex(index=ORDER, columns=MOB_ORDER)
    npivot = joint.pivot(index="compartment", columns="mob_mobility",
                         values="n").reindex(index=ORDER, columns=MOB_ORDER)
    print("\n== AMR prevalence, compartment x mobility ==")
    print((pivot * 100).round(1).to_string())

    # --- 4. logistic regression -------------------------------------------------------
    reg = mob[mob.size_bp > 0].copy()
    reg["amr_pos"] = reg.amr_pos.astype(int)  # statsmodels reads bool endog as categorical
    reg["log_size"] = np.log10(reg.size_bp)
    reg["comp"] = reg.compartment.astype(str)
    reg["mobc"] = reg.mob_mobility.astype(str)
    model = smf.logit(
        f'amr_pos ~ C(comp, Treatment(reference="{REF}")) '
        f'+ C(mobc, Treatment(reference="non-mobilizable")) + log_size + gc_percent',
        data=reg).fit(disp=0)
    res = pd.DataFrame({"term": model.params.index, "coef": model.params.values,
                        "p": model.pvalues.values})
    ci = model.conf_int()
    res["or"] = np.exp(res.coef)
    res["or_lo"] = np.exp(ci[0].values)
    res["or_hi"] = np.exp(ci[1].values)
    res.to_csv(f"{PP}/onehealth_logit_odds_ratios.tsv", sep="\t", index=False)
    print(f"\n== Logistic regression (n={len(reg):,}, pseudo-R2={model.prsquared:.3f}) ==")
    print(res[["term", "or", "or_lo", "or_hi", "p"]].to_string(index=False))

    # --- 5. drug-class spectrum per compartment ---------------------------------------
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
            dc_rows.append({"compartment": g, "drug_class": c, "n_plasmids": k,
                            "frac_of_carriers": k / len(s)})
    dc = pd.DataFrame(dc_rows)
    dc.to_csv(f"{PP}/onehealth_drug_class_by_compartment.tsv", sep="\t", index=False)
    top_classes = (dc.groupby("drug_class").n_plasmids.sum()
                     .sort_values(ascending=False).head(15).index.tolist())
    dc_pivot = (dc[dc.drug_class.isin(top_classes)]
                .pivot(index="drug_class", columns="compartment", values="frac_of_carriers")
                .reindex(index=top_classes, columns=ORDER))

    # --- 6. cross-compartment ARG sharing (Jaccard over ARO sets) ---------------------
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
    print("\n== unique ARO richness per compartment ==")
    print(rich.to_string(index=False))
    print("\n== ARG-repertoire Jaccard (vs Human-clinical) ==")
    print(jac["Human - clinical"].round(3).to_string())

    # --- 7. multireplicon vs AMR ------------------------------------------------------
    mr = oh[oh.pf_n_inc.notna()].copy()
    mr["multireplicon"] = mr.pf_n_inc >= 2
    mr_rows = []
    for g in ORDER:
        s = mr[mr.compartment == g]
        a = s[s.multireplicon]
        b = s[~s.multireplicon]
        if len(a) == 0 or len(b) == 0:
            continue
        orv, pv = stats.fisher_exact(
            [[int(a.amr_pos.sum()), len(a) - int(a.amr_pos.sum())],
             [int(b.amr_pos.sum()), len(b) - int(b.amr_pos.sum())]])
        mr_rows.append({"compartment": g, "n_multireplicon": len(a),
                        "amr_prev_multireplicon": a.amr_pos.mean(),
                        "n_single_or_none": len(b),
                        "amr_prev_other": b.amr_pos.mean(),
                        "odds_ratio": orv, "p_raw": pv})
    mrt = pd.DataFrame(mr_rows)
    mrt["p_bh"] = multipletests(mrt.p_raw, method="fdr_bh")[1]
    mrt.to_csv(f"{PP}/onehealth_multireplicon_amr.tsv", sep="\t", index=False)
    print("\n== multireplicon (PlasmidFinder >=2 Inc) vs AMR, per compartment ==")
    print(mrt.to_string(index=False))

    # --- 8. robustness: one plasmid per MOB cluster ----------------------------------
    dedup = (oh[oh.mob_cluster.notna()]
             .groupby("mob_cluster", observed=True)
             .sample(n=1, random_state=SEED))
    ded_tab = prevalence_table(dedup, "compartment", ORDER)[
        ["compartment", "n", "amr_prevalence"]].rename(
        columns={"n": "n_dedup", "amr_prevalence": "amr_prevalence_dedup"})
    comp_cmp = comp_tab[["compartment", "n", "amr_prevalence"]].merge(
        ded_tab, on="compartment", how="left")
    comp_cmp.to_csv(f"{PP}/onehealth_dereplication_sensitivity.tsv", sep="\t", index=False)
    print(f"\n== dereplication robustness (1 per MOB cluster, n={len(dedup):,}) ==")
    print(comp_cmp.to_string(index=False))

    # ---------------------------------------------------------------- figures --------
    colors = plt.cm.tab10(np.linspace(0, 1, 10))

    # fig 1: prevalence by compartment
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
    plt.close(fig)

    # fig 2: compartment x mobility heatmap
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    im = ax.imshow(pivot.values * 100, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(MOB_ORDER)))
    ax.set_xticklabels(MOB_ORDER, rotation=20, ha="right")
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels(ORDER)
    for i in range(len(ORDER)):
        for j in range(len(MOB_ORDER)):
            val = pivot.values[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val*100:.0f}%\nn={int(npivot.values[i,j]):,}",
                        ha="center", va="center", fontsize=6.5,
                        color="white" if val > 0.45 else "black")
    fig.colorbar(im, ax=ax, label="AMR prevalence (%)")
    ax.set_title("AMR prevalence: compartment x mobility")
    fig.savefig(f"{FIGDIR}/amr_compartment_x_mobility.png")
    plt.close(fig)

    # fig 3: drug-class spectrum
    fig, ax = plt.subplots(figsize=(7.6, 5))
    im = ax.imshow(dc_pivot.values, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=7.5)
    ax.set_yticks(range(len(dc_pivot)))
    ax.set_yticklabels(dc_pivot.index, fontsize=7.5)
    fig.colorbar(im, ax=ax, label="fraction of AMR+ plasmids in compartment")
    ax.set_title("Drug-class spectrum of the plasmid resistome, by compartment")
    fig.savefig(f"{FIGDIR}/drug_class_by_compartment.png")
    plt.close(fig)

    # fig 4: ARG-repertoire Jaccard
    fig, ax = plt.subplots(figsize=(6.2, 5))
    im = ax.imshow(jac.values.astype(float), cmap="magma", vmin=0, vmax=1)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(ORDER, rotation=35, ha="right", fontsize=7.5)
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels(ORDER, fontsize=7.5)
    for i in range(len(ORDER)):
        for j in range(len(ORDER)):
            v_ = jac.values[i, j]
            if not np.isnan(v_):
                ax.text(j, i, f"{v_:.2f}", ha="center", va="center", fontsize=6.5,
                        color="white" if v_ < 0.6 else "black")
    fig.colorbar(im, ax=ax, label="Jaccard index (shared AROs)")
    ax.set_title("Cross-compartment sharing of ARG repertoires")
    fig.savefig(f"{FIGDIR}/arg_repertoire_jaccard.png")
    plt.close(fig)

    # fig 5: odds-ratio forest plot
    plot_terms = res[res.term.str.contains(r"C\(comp|C\(mobc")].copy()
    plot_terms["label"] = (plot_terms.term
                           .str.replace(r".*\[T\.", "", regex=True)
                           .str.replace(r"\]$", "", regex=True))
    fig, ax = plt.subplots(figsize=(6.4, 4))
    yy = np.arange(len(plot_terms))
    ax.errorbar(plot_terms["or"], yy,
                xerr=[plot_terms["or"] - plot_terms.or_lo,
                      plot_terms.or_hi - plot_terms["or"]],
                fmt="o", color="#22456b", ecolor="0.4", capsize=3)
    ax.axvline(1, color="crimson", ls="--", lw=1)
    ax.set_yticks(yy)
    ax.set_yticklabels(plot_terms.label, fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("adjusted odds ratio for carrying >=1 ARG (log scale)")
    ax.set_title(f"AMR odds, adjusted for plasmid size and GC\n"
                 f"(ref: {REF} / non-mobilizable)")
    fig.savefig(f"{FIGDIR}/amr_adjusted_odds_ratios.png")
    plt.close(fig)

    # fig 6: dereplication sensitivity
    fig, ax = plt.subplots(figsize=(7, 3.8))
    xx = np.arange(len(comp_cmp))
    ax.bar(xx - 0.2, comp_cmp.amr_prevalence * 100, 0.4, label="all plasmids")
    ax.bar(xx + 0.2, comp_cmp.amr_prevalence_dedup * 100, 0.4,
           label="1 per MOB cluster")
    ax.set_xticks(xx)
    ax.set_xticklabels(comp_cmp.compartment, rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel("AMR prevalence (%)")
    ax.legend(frameon=False)
    ax.set_title("Clonal-redundancy robustness check")
    fig.savefig(f"{FIGDIR}/dereplication_sensitivity.png")
    plt.close(fig)

    print(f"\nwrote 10 tables to {PP}/onehealth_*.tsv and 6 figures to {FIGDIR}/")


if __name__ == "__main__":
    main()
