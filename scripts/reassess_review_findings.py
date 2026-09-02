#!/usr/bin/env python3
"""Reassessment of the round-1 peer-review findings on manuscript/amr_onehealth_NAR.md.

Executes the computable items of the revision roadmap in `manuscript/review_round1.md`:

  R1 Host coverage        — can the host-taxon confounder (D1/DA3) be tested at all?
  R2 Collinearity         — VIF and mobility-coefficient stability with/without log size (M2)
  R3 Cluster-robust SEs   — refit with SEs clustered on MOB cluster (M1)
  R4 Host-taxon adjust    — compartment coefficients before/after host taxon (D1)
  R5 Headline contrast    — adjusted predicted probabilities at fixed size (M3/DA1)
  R6 Assignment bias      — AMR prevalence in compartment-assigned vs unassigned plasmids (P2)
  R7 Typing coverage      — PlasmidFinder coverage per compartment (D3, multireplicon validity)
  R8 is_clinical basis    — how much of the clinical flag rests on PLSDB disease tags (D4)

Host taxonomy is joined from PlasmidScope `all_metadata.tsv` via the provenance representative key,
following the join already used by `scripts/build_analysis_table.py`.

Inputs:  data/plasmidscope_primary/{plasmid_metadata_master.tsv, all_metadata.tsv,
         complete_provenance.tsv, environment_reconciled.tsv}
Outputs: data/plasmidscope_primary/reassess_*.tsv (6 tables)
         reports/figures/onehealth/reassess_compartment_or_shift.png
Run:     python3 scripts/reassess_review_findings.py
"""
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor

csv.field_size_limit(sys.maxsize)

PP = "data/plasmidscope_primary"
FIGDIR = "reports/figures/onehealth"

HUMAN_SUBS = {
    "Human: unspecified", "Human: gut/faeces", "Human: blood", "Human: urine",
    "Human: respiratory", "Human: skin/wound", "Human: other clinical", "Human: oral",
}
COMPARTMENT_MAP = {
    "Livestock: pig": "Animal - livestock & poultry",
    "Livestock: cattle": "Animal - livestock & poultry",
    "Livestock: other": "Animal - livestock & poultry",
    "Poultry/bird": "Animal - livestock & poultry",
    "Companion animal": "Animal - companion & aquaculture",
    "Fish/aquaculture": "Animal - companion & aquaculture",
    "Other mammal": "Animal - wildlife & other",
    "Rodent": "Animal - wildlife & other",
    "Insect/arthropod": "Animal - wildlife & other",
    "Other invertebrate": "Animal - wildlife & other",
    "Gut/faeces: unspecified host": "Animal - wildlife & other",
    "Wastewater/sewage": "Environment - engineered interface",
    "Bioreactor": "Environment - engineered interface",
    "Built environment": "Environment - engineered interface",
    "Solid waste/compost": "Environment - engineered interface",
    "Industrial/remediation": "Environment - engineered interface",
    "Terrestrial/soil": "Environment - natural",
    "Aquatic: freshwater/other": "Environment - natural",
    "Aquatic: marine": "Environment - natural",
    "Air": "Environment - natural",
    "Extreme/other": "Environment - natural",
    "Environmental: unspecified": "Environment - natural",
    "Plant": "Plant & food",
    "Food/fermentation": "Plant & food",
    "Algae": "Plant & food",
    "Fungi": "Plant & food",
}
ORDER = [
    "Human - clinical", "Human - community",
    "Animal - livestock & poultry", "Animal - companion & aquaculture",
    "Animal - wildlife & other",
    "Environment - engineered interface", "Environment - natural",
    "Plant & food",
]
REF = "Environment - natural"

# Enterobacteriaceae (sensu lato, as used for the detection-bias argument) genera.
ENTERO = {
    "Escherichia", "Klebsiella", "Salmonella", "Enterobacter", "Citrobacter", "Serratia",
    "Proteus", "Shigella", "Cronobacter", "Providencia", "Morganella", "Yersinia",
    "Kluyvera", "Raoultella", "Leclercia", "Edwardsiella", "Hafnia", "Pluralibacter",
    "Lelliottia", "Atlantibacter", "Phytobacter", "Kosakonia",
}
NULL_HOST = {"-", "", "synthetic construct", "unclassified", "uncultured bacterium"}


def compartment(row):
    if row["hab_sub"] in HUMAN_SUBS:
        return "Human - clinical" if row["is_clinical"] == 1 else "Human - community"
    return COMPARTMENT_MAP.get(row["hab_sub"])


def or_table(model, label):
    ci = model.conf_int()
    out = pd.DataFrame({
        "term": model.params.index, "or": np.exp(model.params.values),
        "or_lo": np.exp(ci[0].values), "or_hi": np.exp(ci[1].values),
        "p": model.pvalues.values,
    })
    out["model"] = label
    out["k"] = (out.term.str.replace(r".*\[T\.", "", regex=True)
                        .str.replace(r"\]$", "", regex=True))
    return out


def main():
    plt.rcParams.update({"figure.dpi": 150, "savefig.bbox": "tight", "font.size": 9})
    df = pd.read_csv(f"{PP}/plasmid_metadata_master.tsv", sep="\t", low_memory=False)
    df["amr_pos"] = (df.card_n_arg > 0).astype(int)
    df["compartment"] = df.apply(compartment, axis=1)

    # ---- join PlasmidScope host taxonomy via the provenance representative key --------
    rep2full = {}
    with open(f"{PP}/complete_provenance.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rep2full[r["plasmid_id"].split(",")[0]] = r["plasmid_id"]
    host = {}
    with open(f"{PP}/all_metadata.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            host[r["Plasmid_ID"]] = r.get("Host", "")
    df["host_raw"] = df.plasmid_id.map(lambda p: host.get(rep2full.get(p, ""), ""))
    df["host_raw"] = df.host_raw.fillna("")
    df["has_host"] = ~df.host_raw.str.strip().isin(NULL_HOST)
    df["host_genus"] = np.where(df.has_host, df.host_raw.str.split().str[0], None)
    df["is_entero"] = df.host_genus.isin(ENTERO).astype(int)

    oh = df[df.compartment.notna()].copy()
    print(f"working set {len(df):,}; One Health set {len(oh):,}")

    # =========== R1: host coverage ====================================================
    print("\n=== R1 host-taxonomy coverage by compartment ===")
    cov = (oh.groupby("compartment")
             .agg(n=("amr_pos", "size"), with_host=("has_host", "sum"),
                  entero_frac=("is_entero", "mean"))
             .reindex(ORDER))
    cov["host_coverage"] = cov.with_host / cov.n
    cov["entero_frac_of_hosted"] = (
        oh[oh.has_host].groupby("compartment").is_entero.mean().reindex(ORDER))
    cov.to_csv(f"{PP}/reassess_host_coverage.tsv", sep="\t")
    print(cov[["n", "with_host", "host_coverage", "entero_frac_of_hosted"]].round(3).to_string())

    # =========== R6: assignment bias ==================================================
    print("\n=== R6 assignment bias: AMR prevalence assigned vs unassigned ===")
    df["assigned"] = df.compartment.notna()
    ab = df.groupby("assigned").agg(n=("amr_pos", "size"), amr=("amr_pos", "sum"))
    ab["prevalence"] = ab.amr / ab.n
    # unassigned split by reason
    df["unassigned_kind"] = np.where(
        df.assigned, "assigned",
        np.where(df.hab_sub == "Simulated communities", "simulated",
                 np.where(df.hab_sub.isin(["Unlabelled", "Geography only (no habitat)",
                                           "Unresolved"]), "no label", "other excluded")))
    ak = df.groupby("unassigned_kind").agg(n=("amr_pos", "size"), amr=("amr_pos", "sum"))
    ak["prevalence"] = ak.amr / ak.n
    ak.to_csv(f"{PP}/reassess_assignment_bias.tsv", sep="\t")
    print(ak.round(4).to_string())

    # =========== R7: PlasmidFinder typing coverage ====================================
    print("\n=== R7 PlasmidFinder typing coverage by compartment ===")
    oh["typed"] = (oh.pf_n_inc.fillna(0) >= 1).astype(int)
    tc = (oh.groupby("compartment")
            .agg(n=("typed", "size"), typed=("typed", "sum")).reindex(ORDER))
    tc["typing_coverage"] = tc.typed / tc.n
    tc.to_csv(f"{PP}/reassess_typing_coverage.tsv", sep="\t")
    print(tc.round(3).to_string())

    # =========== R8: is_clinical basis ================================================
    print("\n=== R8 clinical-flag basis ===")
    clin = oh[oh.compartment == "Human - clinical"]
    body_site = clin.hab_sub.isin(["Human: blood", "Human: urine", "Human: respiratory",
                                   "Human: skin/wound", "Human: other clinical"])
    has_disease = clin.plsdb_disease_tags.notna() & (clin.plsdb_disease_tags != "")
    basis = pd.DataFrame({
        "basis": ["clinical body site", "no body site, has PLSDB disease tag",
                  "no body site, no disease tag (context only)"],
        "n": [int(body_site.sum()),
              int((~body_site & has_disease).sum()),
              int((~body_site & ~has_disease).sum())],
    })
    basis["frac"] = basis.n / len(clin)
    basis.to_csv(f"{PP}/reassess_clinical_basis.tsv", sep="\t", index=False)
    print(basis.round(3).to_string(index=False))

    # =========== modelling set =========================================================
    reg = oh[oh.mob_mobility.notna() & (oh.size_bp > 0)].copy()
    reg["log_size"] = np.log10(reg.size_bp)
    reg["comp"] = reg.compartment.astype(str)
    reg["mobc"] = reg.mob_mobility.astype(str)
    reg["cluster"] = reg.mob_cluster.fillna("NA")

    BASE = (f'amr_pos ~ C(comp, Treatment(reference="{REF}")) '
            f'+ C(mobc, Treatment(reference="non-mobilizable")) + log_size + gc_percent')

    # =========== R2: collinearity =====================================================
    print("\n=== R2 collinearity diagnostics ===")
    dm = pd.get_dummies(reg[["comp", "mobc"]], drop_first=True).astype(float)
    dm["log_size"] = reg.log_size.values
    dm["gc_percent"] = reg.gc_percent.values
    dm.insert(0, "const", 1.0)
    vif = pd.DataFrame({
        "term": dm.columns,
        "VIF": [variance_inflation_factor(dm.values, i) for i in range(dm.shape[1])],
    })
    vif = vif[vif.term != "const"]
    vif.to_csv(f"{PP}/reassess_vif.tsv", sep="\t", index=False)
    print(vif.sort_values("VIF", ascending=False).round(2).to_string(index=False))

    m_full = smf.logit(BASE, data=reg).fit(disp=0)
    m_nosize = smf.logit(BASE.replace(" + log_size", ""), data=reg).fit(disp=0)
    print("\nmobility OR with vs without log_size in the model:")
    for k in ["conjugative", "mobilizable"]:
        a = or_table(m_full, "with size").set_index("k").loc[k]
        b = or_table(m_nosize, "without size").set_index("k").loc[k]
        print(f"  {k:16s} with size OR={a['or']:.2f}  without size OR={b['or']:.2f}")

    # =========== R3: cluster-robust SEs ===============================================
    print("\n=== R3 cluster-robust standard errors (groups = MOB cluster) ===")
    m_clu = smf.logit(BASE, data=reg).fit(
        disp=0, cov_type="cluster", cov_kwds={"groups": reg["cluster"]})
    a = or_table(m_full, "nominal").set_index("k")
    b = or_table(m_clu, "clustered").set_index("k")
    cmp_rows = []
    for k in a.index:
        if k == "Intercept":
            continue
        cmp_rows.append({
            "term": k, "or": round(a.loc[k, "or"], 3),
            "nominal_lo": round(a.loc[k, "or_lo"], 3), "nominal_hi": round(a.loc[k, "or_hi"], 3),
            "clustered_lo": round(b.loc[k, "or_lo"], 3), "clustered_hi": round(b.loc[k, "or_hi"], 3),
            "nominal_p": a.loc[k, "p"], "clustered_p": b.loc[k, "p"],
            "ci_width_ratio": round(
                (b.loc[k, "or_hi"] - b.loc[k, "or_lo"]) / (a.loc[k, "or_hi"] - a.loc[k, "or_lo"]), 1),
            "still_sig_0.05": bool(b.loc[k, "p"] < 0.05),
        })
    clu = pd.DataFrame(cmp_rows)
    clu.to_csv(f"{PP}/reassess_cluster_robust.tsv", sep="\t", index=False)
    print(clu.to_string(index=False))
    print(f"n clusters = {reg['cluster'].nunique():,}")

    # =========== R4: host-taxon adjustment ============================================
    print("\n=== R4 compartment ORs before vs after host-taxon adjustment ===")
    hs = reg[reg.has_host].copy()
    top_genera = hs.host_genus.value_counts().head(20).index.tolist()
    hs["genus_grp"] = np.where(hs.host_genus.isin(top_genera), hs.host_genus, "other")
    m_h0 = smf.logit(BASE, data=hs).fit(disp=0)                      # same subset, no host
    m_h1 = smf.logit(BASE + " + is_entero", data=hs).fit(disp=0)     # + Entero indicator
    m_h2 = smf.logit(BASE + ' + C(genus_grp, Treatment(reference="other"))',
                     data=hs).fit(disp=0)                            # + top-20 genus
    t0 = or_table(m_h0, "no host").set_index("k")
    t1 = or_table(m_h1, "+entero").set_index("k")
    t2 = or_table(m_h2, "+genus").set_index("k")
    rows = []
    for k in [c for c in ORDER if c != REF] + ["conjugative", "mobilizable", "log_size"]:
        if k not in t0.index:
            continue
        rows.append({
            "term": k,
            "OR_no_host": round(t0.loc[k, "or"], 2),
            "OR_plus_entero": round(t1.loc[k, "or"], 2),
            "OR_plus_genus": round(t2.loc[k, "or"], 2),
            "pct_change_genus": round(100 * (t2.loc[k, "or"] - t0.loc[k, "or"]) / t0.loc[k, "or"], 1),
            "p_genus_model": t2.loc[k, "p"],
        })
    host_cmp = pd.DataFrame(rows)
    host_cmp.to_csv(f"{PP}/reassess_host_adjustment.tsv", sep="\t", index=False)
    print(f"(subset with host taxonomy: n={len(hs):,})")
    print(host_cmp.to_string(index=False))
    print(f"pseudo-R2: no host {m_h0.prsquared:.3f} | +entero {m_h1.prsquared:.3f} "
          f"| +genus {m_h2.prsquared:.3f}")

    # =========== R5: adjusted headline contrast =======================================
    print("\n=== R5 headline contrast, adjusted predicted probabilities ===")
    med_gc = float(reg.gc_percent.median())
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
    piv = sc.pivot_table(index="size_bp", columns=["comp", "mobc"], values="pred_prob")
    sc.to_csv(f"{PP}/reassess_headline_contrast.tsv", sep="\t", index=False)
    print("predicted P(AMR+) at matched plasmid size, median GC:")
    print((piv * 100).round(1).to_string())

    # observed median sizes of the two headline cells
    cell_a = reg[(reg.comp == "Environment - natural") & (reg.mobc == "conjugative")]
    cell_b = reg[(reg.comp == "Human - clinical") & (reg.mobc == "non-mobilizable")]
    print(f"\nobserved: env-natural/conjugative  n={len(cell_a):,} median={cell_a.size_bp.median():,.0f} bp "
          f"raw prevalence={cell_a.amr_pos.mean():.3f}")
    print(f"observed: clinical/non-mobilizable n={len(cell_b):,} median={cell_b.size_bp.median():,.0f} bp "
          f"raw prevalence={cell_b.amr_pos.mean():.3f}")

    # =========== figure: compartment OR shift under host adjustment ===================
    plot = host_cmp[host_cmp.term.isin([c for c in ORDER if c != REF])]
    fig, ax = plt.subplots(figsize=(7, 4))
    y = np.arange(len(plot))
    ax.barh(y - 0.2, plot.OR_no_host, 0.4, label="compartment only")
    ax.barh(y + 0.2, plot.OR_plus_genus, 0.4, label="+ host genus")
    ax.axvline(1, color="crimson", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(plot.term, fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("adjusted odds ratio for carrying >=1 ARG (log scale)")
    ax.set_title("Compartment effects before and after host-taxon adjustment")
    ax.legend(frameon=False)
    fig.savefig(f"{FIGDIR}/reassess_compartment_or_shift.png")
    plt.close(fig)

    print(f"\nwrote 6 tables to {PP}/reassess_*.tsv and 1 figure to {FIGDIR}/")


if __name__ == "__main__":
    main()
