import pandas as pd, numpy as np, sys
R="/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis/"
m=pd.read_csv(R+"data/plasmidscope_primary/plasmid_metadata_master.tsv",sep="\t",low_memory=False)
def P(k,got,claim=None):
    ok="" if claim is None else ("  OK" if str(got)==str(claim) else f"  <<< CLAIM {claim}")
    print(f"{k:<62} {got}{ok}")

print("### A. RESOURCE")
P("master rows",len(m),208248)
P("master cols",m.shape[1],67)
P("unique plasmid_id",m.plasmid_id.nunique(),208248)

print("\n-- hab_top --")
print(m.hab_top.value_counts(dropna=False).to_string())
sim=(m.hab_top=="Simulated-artifact").sum(); lab=(m.hab_top=="Lab-artifact").sum()
P("Simulated-artifact",sim,64658); P("Lab-artifact",lab,87)
ana=m[~m.hab_top.isin(["Simulated-artifact","Lab-artifact"])]
P("analysis set (excl sim+lab)",len(ana),143503)
P("excl sim ONLY",len(m)-sim,143590)
P("is_clinical (full set)",int(m.is_clinical.fillna(0).astype(float).sum()),19446)
P("is_clinical (analysis set)",int(ana.is_clinical.fillna(0).astype(float).sum()))

print("\n-- AMR (full working set) --")
arg=m.card_n_arg.fillna(0)
P("AMR+ card_n_arg>0",int((arg>0).sum()),29396)
P("AMR+ pct",round(100*(arg>0).mean(),2),14.1)
P("card_multidrug",int(m.card_multidrug.fillna(0).astype(float).sum()),20741)

print("\n-- typing --")
P("pf_n_inc>=1",int((m.pf_n_inc.fillna(0)>=1).sum()),41555)
P("pf_n_inc>=2 (multireplicon)",int((m.pf_n_inc.fillna(0)>=2).sum()),24583)
P("pf pct of 208248",round(100*(m.pf_n_inc.fillna(0)>=1).mean(),1),20)
P("mob_cluster non-null",int(m.mob_cluster.notna().sum()))
P("mob_cluster distinct",m.mob_cluster.nunique(),7054)
P("mob_host_range non-null pct",round(100*m.mob_host_range.notna().mean(),1),58)
anyrep = (m.pf_n_inc.fillna(0)>=1) | m.mob_rep_types.notna() | (m.plasann_n_replicons.fillna(0)>=1)
P("any-method replicon call",int(anyrep.sum()),95442)
P("any-method pct",round(100*anyrep.mean(),1),46)
