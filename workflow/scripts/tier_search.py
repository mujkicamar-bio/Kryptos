import _ctx  # noqa: F401
import csv, subprocess, tempfile, pathlib
from plasmidann.cascade import explained_fraction, narrow_by_explained, is_informative

spec = snakemake.params.spec
faa = snakemake.input.faa
ids = [l[1:].split()[0] for l in open(faa) if l[0] == ">"]

# spans carried forward hold ONLY informative alignments: a hypothetical hit explains
# nothing and must never satisfy the narrowing threshold.
spans, qlen = {}, {}
prev = snakemake.input.spans
for f in ([prev] if isinstance(prev, str) else list(prev)):
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            qlen[r["seq_id"]] = int(r["qlen"])
            spans[r["seq_id"]] = [tuple(map(int, iv.split("-")))
                                  for iv in filter(None, r["intervals"].split(";"))]

best = {}     # best informative hit per query
unnamed = {}  # every uninformative hit, kept as evidence someone else has seen it


def record(q, label, qcov, tcov, ev, start, end, tlen):
    qlen[q] = tlen
    if is_informative(label):
        spans.setdefault(q, []).append((start, end))
        if q not in best or qcov > best[q]["coverage"]:
            best[q] = {"query": q, "label": label, "coverage": round(qcov, 4),
                       "target_coverage": round(tcov, 4), "evalue": ev,
                       "informative": True}
    else:
        unnamed.setdefault(q, []).append(
            {"query": q, "label": label, "coverage": round(qcov, 4),
             "target_coverage": round(tcov, 4), "evalue": ev, "informative": False})


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
            tlen, hlen = int(f[2]), int(f[5])
            a, b = int(f[17]), int(f[18])
            record(f[0], f[3], (b - a + 1) / tlen if tlen else 0.0,
                   (int(f[16]) - int(f[15]) + 1) / hlen if hlen else 0.0, f[6], a, b, tlen)
    else:
        raw = f"{tmp}/res.m8"
        subprocess.run(f"diamond blastp -q {faa} -d {spec['db']} -o {raw} "
                       f"{spec['args']} --threads {snakemake.threads} --max-target-seqs 5 "
                       f"--outfmt 6 qseqid stitle qcovhsp scovhsp evalue qstart qend qlen "
                       f"--quiet", shell=True, check=True)
        for line in open(raw):
            q, title, qc, tc, ev, qs, qe, ql = line.rstrip("\n").split("\t")
            record(q, title, float(qc) / 100, float(tc) / 100, ev,
                   int(qs), int(qe), int(ql))

cols = ["query", "label", "coverage", "target_coverage", "evalue", "informative",
        "tier", "threshold"]
with open(snakemake.output.hits, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for h in list(best.values()) + [x for v in unnamed.values() for x in v]:
        w.writerow({**h, "tier": spec["id"], "threshold": spec["args"]})

explained = {q: explained_fraction(qlen.get(q, 0), spans.get(q, [])) for q in qlen}
with open(snakemake.output.spans, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["seq_id", "qlen", "intervals", "explained_fraction"])
    for q in qlen:
        w.writerow([q, qlen[q], ";".join(f"{a}-{b}" for a, b in spans.get(q, [])),
                    explained.get(q, 0.0)])

keep = set(narrow_by_explained(ids, {k: v for k, v in explained.items() if k in set(ids)},
                               snakemake.params.min_explained))
with open(snakemake.output.unresolved, "w") as out:
    emit = False
    for line in open(faa):
        if line[0] == ">":
            emit = line[1:].split()[0] in keep
        if emit:
            out.write(line)
