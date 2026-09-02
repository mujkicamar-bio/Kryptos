import pandas as pd, numpy as np
R="/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis/"
m=pd.read_csv(R+"data/plasmidscope_primary/plasmid_metadata_master.tsv",sep="\t",low_memory=False,
              usecols=["plasmid_id","pf_inc_types","pf_inc_families","pf_n_inc","size_bp","mob_rep_types","plasann_n_replicons"])
t=m[m.pf_n_inc.fillna(0)>=1]
print("typed:",len(t))
print("\n-- pf_n_inc distribution --"); print(t.pf_n_inc.value_counts().sort_index().head(15).to_string())
print("\n-- examples of pf_n_inc>=2 --")
print(t[t.pf_n_inc>=2][["plasmid_id","size_bp","pf_inc_types","pf_inc_families","pf_n_inc"]].head(12).to_string())
# distinct FAMILIES rather than allele calls
def nfam(s):
    if pd.isna(s): return 0
    return len(set(x.strip() for x in str(s).split(";") if x.strip()))
t=t.copy()
t["n_fam"]=t.pf_inc_families.map(nfam)
t["n_typ"]=t.pf_inc_types.map(nfam)
print("\nmultireplicon by pf_n_inc>=2 :",(t.pf_n_inc>=2).sum())
print("multireplicon by distinct pf_inc_types>=2 :",(t.n_typ>=2).sum())
print("multireplicon by distinct pf_inc_FAMILIES>=2:",(t.n_fam>=2).sum())
print("\nmedian size, single vs multi (by family):")
print(t.groupby(t.n_fam>=2).size_bp.median().to_string())
