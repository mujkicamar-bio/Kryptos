import pandas as pd, numpy as np
R="/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis/"
R2=R+"data/dark_orf_run/recluster/"
claims={ # threshold -> (families, pct1plasmid, recurrent%, coreN)
 "30":(99362,62.9,83.3,61),"50":(117393,64.5,79.8,39),"70":(137149,66.7,75.7,19),"90":(160518,69.6,70.3,5)}
print(f"{'thr':<5}{'families':>10}{'claim':>10}{'recurr%':>10}{'claim':>8}{'1plasmid%':>12}{'claim':>8}{'core':>7}{'claim':>7}")
for t in ["30","50","70","90"]:
    cl=pd.read_csv(R2+f"dark{t}_cluster.tsv",sep="\t",header=None,names=["rep","mem"])
    fam=cl.groupby("rep").size()
    fs=pd.read_csv(R2+f"dark{t}_famstats.tsv",sep="\t")
    fp=pd.read_csv(R2+f"dark{t}_fampfam.tsv",sep="\t")
    mg=fs.merge(fp[["rep","pfam_hit"]],on="rep",how="left")
    core=mg[(mg.pfam_hit!=True)&(mg.plasmids>=100)&(mg.lineages>=10)]
    rec=100*fam[fam>=2].sum()/len(cl)
    p1=100*(fs.plasmids==1).mean()
    c=claims[t]
    print(f"{t:<5}{len(fam):>10}{c[0]:>10}{rec:>10.1f}{c[2]:>8}{p1:>12.1f}{c[1]:>8}{len(core):>7}{c[3]:>7}")
    if t=="30":
        print("   members sum:",len(cl),"unique:",cl.mem.nunique())
# original 30% core = 93
fs=pd.read_csv(R+"data/dark_orf_run/darkfam_stats.tsv",sep="\t")
fp=pd.read_csv(R+"data/pfam_run/darkfam_pfam.tsv",sep="\t")
mg=fs.merge(fp[["rep","pfam_hit"]],on="rep",how="left")
core=mg[(mg.pfam_hit!=True)&(mg.plasmids>=100)&(mg.lineages>=10)]
print("\noriginal 30% (no reassign) core families:",len(core),"CLAIM 93")
print("core % of dark ORFs:",round(100*core.proteins.sum()/378552,1),"CLAIM 5.9")
