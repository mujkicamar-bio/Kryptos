import pandas as pd, numpy as np, collections, gzip
R="/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis/"
D=R+"data/dark_orf_run/"
def P(k,got,claim=None):
    ok="" if claim is None else ("  OK" if str(got)==str(claim) else f"  <<< CLAIM {claim}")
    print(f"{k:<60} {got}{ok}")

# 1. FASTA counts
def nheaders(fp):
    n=0
    with open(fp) as fh:
        for l in fh:
            if l[0]=='>': n+=1
    return n
P("dark_orfs.faa seqs",nheaders(D+"dark_orfs.faa"),378552)
P("all_cds.faa seqs",nheaders(D+"all_cds.faa"),456949)

# 2. dark30 clustering
cl=pd.read_csv(D+"dark30_cluster.tsv",sep="\t",header=None,names=["rep","mem"])
P("dark30 member rows",len(cl),378552)
P("dark30 members unique",cl.mem.nunique(),378552)
fam=cl.groupby("rep").size()
P("dark30 families",len(fam),92752)
P("singletons",int((fam==1).sum()),57935)
P("ORFs in families >=2",int(fam[fam>=2].sum()),320617)
P("recurrence pct",round(100*fam[fam>=2].sum()/len(cl),1),84.7)
P("families >=10",int((fam>=10).sum()),5623)
P("ORFs in fams>=10",int(fam[fam>=10].sum()),222506)

# 3. mix30 co-clustering -> dark-only share
mix=pd.read_csv(D+"mix30_cluster.tsv",sep="\t",header=None,names=["rep","mem"])
P("mix30 member rows",len(mix),456949)
mix["is_dark"]=mix.mem.str.endswith("|D")
g=mix.groupby("rep").is_dark.agg(['sum','size'])
mixed=g[(g['sum']>0)&(g['sum']<g['size'])]
darkonly=g[(g['sum']==g['size'])&(g['sum']>0)]
P("mix30 families total",len(g),97500)
P("families w/ both dark+named",len(mixed),1674)
P("dark ORFs absorbed by mixed",int(mixed['sum'].sum()),53219)
P("dark-only families",len(darkonly),91159)
P("dark ORFs in dark-only fams",int(darkonly['sum'].sum()),325333)
P("dark-only pct of dark",round(100*darkonly['sum'].sum()/mix.is_dark.sum(),1),85.9)
