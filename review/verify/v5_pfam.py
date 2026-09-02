import pandas as pd, numpy as np
R="/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis/"
def P(k,got,claim=None):
    ok="" if claim is None else ("  OK" if str(got)==str(claim) else f"  <<< CLAIM {claim}")
    print(f"{k:<58} {got}{ok}")
ps=pd.read_csv(R+"data/pfam_run/pfam_per_seq.tsv",sep="\t",low_memory=False)
print("rows in pfam_per_seq:",len(ps))
hit=ps[ps.pfam_name.notna()]
print(hit.groupby("set").size().to_string())
P("NAMED with hit",int((hit.set=="NAMED").sum()),60124)
P("NAMED pct",round(100*(hit.set=="NAMED").sum()/78397,1),76.7)
P("DARKREP with hit",int((hit.set=="DARKREP").sum()),24357)
P("DARKREP pct",round(100*(hit.set=="DARKREP").sum()/92752,1),26.3)

# family-size-weighted
fs=pd.read_csv(R+"data/dark_orf_run/darkfam_stats.tsv",sep="\t")
fp=pd.read_csv(R+"data/pfam_run/darkfam_pfam.tsv",sep="\t")
P("darkfam_stats families",len(fs),92752)
P("sum proteins",int(fs.proteins.sum()),378552)
mg=fs.merge(fp[["rep","pfam_hit"]],on="rep",how="left")
P("weighted ORFs in Pfam-hit families",int(mg[mg.pfam_hit==True].proteins.sum()),149160)
P("weighted pct",round(100*mg[mg.pfam_hit==True].proteins.sum()/378552,1),39.4)
sub=mg[mg.proteins>=10]
P("families >=10 members",len(sub),5623)
P("of those, pct Pfam hit",round(100*(sub.pfam_hit==True).mean(),1),40.5)
# 4.4 deliverable
nov=sub[sub.pfam_hit!=True]
P("Pfam-negative families >=10",len(nov),3343)
P("their dark ORFs",int(nov.proteins.sum()),115811)
P("pct of dark proteome",round(100*nov.proteins.sum()/378552,1),30.6)
dn=pd.read_csv(R+"data/pfam_run/dark_novel_families.tsv",sep="\t")
P("dark_novel_families.tsv rows",len(dn),3343)
P("pct spanning >=2 lineages",round(100*(dn.lineages>=2).mean(),1),72.9)
P("median lineages",dn.lineages.median(),3); P("median habitats",dn.habitats.median(),4)

# 4.2 cross-tab
mix=pd.read_csv(R+"data/dark_orf_run/mix30_cluster.tsv",sep="\t",header=None,names=["rep","mem"])
mix["is_dark"]=mix.mem.str.endswith("|D")
g=mix.groupby("rep").is_dark.agg(['sum','size'])
donly=set(g[(g['sum']==g['size'])&(g['sum']>0)].index)
dark=mix[mix.is_dark].copy()
dark["dark_only"]=dark.rep.isin(donly)
dark["pid"]=dark.mem.str.rsplit("|",n=1).str[0]          # strip |D
cl=pd.read_csv(R+"data/dark_orf_run/dark30_cluster.tsv",sep="\t",header=None,names=["rep30","mem"])
famhit=dict(zip(fp.rep,fp.pfam_hit==True))
cl["pfam"]=cl.rep30.map(famhit)
m=dark.merge(cl[["mem","pfam"]],left_on="pid",right_on="mem",how="left",suffixes=("","_x"))
P("join coverage (should be 378552)",int(m.pfam.notna().sum()),378552)
ct=pd.crosstab(m.dark_only,m.pfam,normalize=True)*100
print("\ncross-tab (% of 378,552), rows=dark_only, cols=Pfam match")
print(ct.round(1).to_string())
print("CLAIM: dark-only/no-Pfam 57.4 | dark-only/Pfam 28.5 | named-homolog/Pfam 10.9 | named-homolog/no-Pfam 3.2")
