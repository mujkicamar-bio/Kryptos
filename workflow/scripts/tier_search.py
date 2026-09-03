import _ctx  # noqa: F401
import csv, subprocess, tempfile, pathlib
from plasmidann.cascade import explained_fraction, narrow_by_explained

spec = snakemake.params.spec
faa = snakemake.input.faa
ids = [l[1:].split()[0] for l in open(faa) if l[0] == ">"]

# carry forward what earlier tiers already explained
spans, qlen = {}, {}
prev = snakemake.input.spans
for f in ([prev] if isinstance(prev, str) else list(prev)):
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            qlen[r["seq_id"]] = int(r["qlen"])
            spans[r["seq_id"]] = [tuple(map(int, iv.split("-")))
                                  for iv in filter(None, r["intervals"].split(";"))]

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
            q, label = f[0], f[3]
            tlen, hlen = int(f[2]), int(f[5])
            qcov = (int(f[18]) - int(f[17]) + 1) / tlen if tlen else 0.0
            tcov = (int(f[16]) - int(f[15]) + 1) / hlen if hlen else 0.0
            spans.setdefault(q, []).append((int(f[17]), int(f[18])))
            qlen[q] = tlen
            if q not in hits or qcov > hits[q]["coverage"]:
                hits[q] = {"query": q, "label": label, "coverage": round(qcov, 4),
                           "target_coverage": round(tcov, 4), "evalue": f[6]}
    else:
        raw = f"{tmp}/res.m8"
        subprocess.run(f"mmseqs easy-search {faa} {spec['db']} {raw} {tmp}/ms "
                       f"{spec['args']} --threads {snakemake.threads} --max-seqs 5 "
                       f"--format-output query,theader,qcov,tcov,evalue,qstart,qend,qlen -v 1",
                       shell=True, check=True, stdout=subprocess.DEVNULL)
        for line in open(raw):
            q, header, qcov, tcov, ev, qs, qe, ql = line.rstrip("\n").split("\t")
            label = header.split(None, 1)[1] if " " in header else header
            spans.setdefault(q, []).append((int(qs), int(qe)))
            qlen[q] = int(ql)
            if q not in hits:
                hits[q] = {"query": q, "label": label, "coverage": float(qcov),
                           "target_coverage": float(tcov), "evalue": ev}

cols = ["query", "label", "coverage", "target_coverage", "evalue", "tier", "threshold"]
with open(snakemake.output.hits, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for h in hits.values():
        w.writerow({**h, "tier": spec["id"], "threshold": spec["args"]})

# cumulative explained fraction decides who keeps going
explained = {q: explained_fraction(qlen.get(q, 0), spans.get(q, [])) for q in spans}
with open(snakemake.output.spans, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["seq_id", "qlen", "intervals", "explained_fraction"])
    for q, iv in spans.items():
        w.writerow([q, qlen.get(q, 0), ";".join(f"{a}-{b}" for a, b in iv), explained[q]])

keep = set(narrow_by_explained(ids, {k: v for k, v in explained.items() if k in set(ids)},
                               snakemake.params.min_explained))
with open(snakemake.output.unresolved, "w") as out:
    emit = False
    for line in open(faa):
        if line[0] == ">":
            emit = line[1:].split()[0] in keep
        if emit:
            out.write(line)
