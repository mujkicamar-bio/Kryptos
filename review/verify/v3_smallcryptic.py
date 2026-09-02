import pandas as pd, numpy as np
R="/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis/"
m=pd.read_csv(R+"data/plasmidscope_primary/plasmid_metadata_master.tsv",sep="\t",low_memory=False)
def P(k,got,claim=None):
    ok="" if claim is None else ("  OK" if str(got)==str(claim) else f"  <<< CLAIM {claim}")
    print(f"{k:<58} {got}{ok}")
a=m[~m.hab_top.isin(["Simulated-artifact","Lab-artifact"])].copy()
P("analysis set",len(a),143503)

# mob-based multireplicon (standing rule)
def n_semi(s):
    if pd.isna(s) or str(s).strip()=="" : return 0
    return len(set(x.strip() for x in str(s).split(",") if x.strip()))
a["mob_nrep"]=a.mob_rep_types.map(n_semi)
print("\n-- mob_rep_types based replicon count (analysis set) --")
print(a.mob_nrep.value_counts().sort_index().head(8).to_string())
P("mob >=1 rep",int((a.mob_nrep>=1).sum()))
P("mob >=2 rep",int((a.mob_nrep>=2).sum()))
P("mob multi % of typed",round(100*(a.mob_nrep>=2).sum()/max(1,(a.mob_nrep>=1).sum()),1))
P("mob multi % of analysis set",round(100*(a.mob_nrep>=2).mean(),1),24.7)

print("\n### B. SMALL CRYPTIC (report: size<10kb)")
small=a[a.size_bp<10000]
P("small <10kb",len(small),70243)
P("small pct",round(100*len(small)/len(a),1),48.9)
z=lambda c: a[c].fillna(0)
cryptic_mask=(z("card_n_arg")==0)&(z("plasann_n_virulence")==0)&(z("plasann_n_metal_biocide")==0)&(z("plasann_n_conjugation")==0)
sc=a[(a.size_bp<10000)&cryptic_mask]
P("small-cryptic",len(sc),47031)
P("small-cryptic % of analysis",round(100*len(sc)/len(a),1),32.8)
P("median size",sc.size_bp.median(),4148)
P("median GC",round(sc.gc_percent.median(),1),43.6)
P("median n_cds",sc.plasann_n_cds.median(),5)

# typing blind spot
def anycall(df):
    return (df.pf_n_inc.fillna(0)>=1)|df.mob_rep_types.notna()|(df.plasann_n_replicons.fillna(0)>=1)
P("small-cryptic PF %",round(100*(sc.pf_n_inc.fillna(0)>=1).mean(),1),10.5)
P("small-cryptic mob_typer %",round(100*sc.mob_rep_types.notna().mean(),1),27.5)
P("small-cryptic PlasAnn %",round(100*(sc.plasann_n_replicons.fillna(0)>=1).mean(),1),17.6)
P("small-cryptic ANY %",round(100*anycall(sc).mean(),1),28.9)
P("small-cryptic untyped n",int((~anycall(sc)).sum()),33427)
P("untyped % of analysis set",round(100*(~anycall(sc)).sum()/len(a),1),23.3)
lc=a[(a.size_bp>=10000)&(a.mob_mobility=="conjugative")]
P("large-conjugative n",len(lc),25631)
P("large-conj ANY %",round(100*anycall(lc).mean(),1),92.2)

print("\n### C. DARK-ORF SET (report: size<20kb, payload-free w/o conjugation filter)")
pf_mask=(z("card_n_arg")==0)&(z("plasann_n_virulence")==0)&(z("plasann_n_metal_biocide")==0)
P("small <20kb",int((a.size_bp<20000).sum()),82261)
dk=a[(a.size_bp<20000)&pf_mask]
P("payload-free small (dark set)",len(dk),71414)
ids=set(open(R+"data/dark_orf_run/cryptic_small_ids.txt").read().split())
P("cryptic_small_ids.txt n",len(ids),71414)
P("recomputed set == id file",len(set(dk.plasmid_id)&ids)==len(ids) and len(dk)==len(ids))
P("562 leak: PlasAnn AMR CDS in dark set",int((dk.plasann_n_amr.fillna(0)>0).sum()),562)
print("\nOVERLAP small-cryptic(47031) vs dark set(71414):",len(set(sc.plasmid_id)&ids))
