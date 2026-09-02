import _ctx  # noqa: F401
import csv, subprocess, tempfile, pathlib
from plasmidann.cascade import narrow

spec, faa = snakemake.params.spec, snakemake.input.faa
ids = [l[1:].split()[0] for l in open(faa) if l[0] == ">"]
hits = {}

if ids:
    tmp = tempfile.mkdtemp(dir=pathlib.Path(snakemake.output.hits).parent)
    if spec["method"] == "hmmer":
        raw = f"{tmp}/dom.tbl"
        subprocess.run(f"hmmsearch {spec['args']} --noali --cpu {snakemake.threads} "
                       f"--domtblout {raw} {spec['db']} {faa}",
                       shell=True, check=True, stdout=subprocess.DEVNULL)
        for line in open(raw):
            if line.startswith("#"):
                continue
            f = line.split()
            q, label, tlen, alen = f[0], f[3], int(f[2]), int(f[18]) - int(f[17]) + 1
            cov = alen / tlen if tlen else 0.0
            if q not in hits or cov > hits[q]["coverage"]:
                hits[q] = {"query": q, "label": label, "coverage": round(cov, 4),
                           "evalue": f[6]}
    else:
        raw = f"{tmp}/res.m8"
        subprocess.run(f"mmseqs easy-search {faa} {spec['db']} {raw} {tmp}/ms "
                       f"{spec['args']} --threads {snakemake.threads} --max-seqs 5 "
                       f"--format-output query,theader,qcov,evalue -v 1",
                       shell=True, check=True, stdout=subprocess.DEVNULL)
        for line in open(raw):
            q, header, qcov, ev = line.rstrip("\n").split("\t")
            label = header.split(None, 1)[1] if " " in header else header
            if q not in hits:
                hits[q] = {"query": q, "label": label, "coverage": float(qcov),
                           "evalue": ev}

cols = ["query", "label", "coverage", "evalue", "tier", "threshold"]
with open(snakemake.output.hits, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for h in hits.values():
        w.writerow({**h, "tier": spec["id"], "threshold": spec["args"]})

keep = set(narrow(ids, set(hits)))
with open(snakemake.output.unresolved, "w") as out:
    emit = False
    for line in open(faa):
        if line[0] == ">":
            emit = line[1:].split()[0] in keep
        if emit:
            out.write(line)
