import _ctx  # noqa: F401
import csv
import pyrodigal

aa = snakemake.params.min_orf_aa
gf = pyrodigal.GeneFinder(meta=True, min_gene=aa * 3,
                          min_edge_gene=min(aa * 3, 60))
with open(snakemake.output.faa, "w") as faa, open(snakemake.output.tsv, "w") as tsv:
    w = csv.writer(tsv, delimiter="\t")
    w.writerow(["plasmid_id", "start", "end", "strand", "partial", "seq"])
    pid, seq = None, []
    def flush():
        if pid is None:
            return
        for g in gf.find_genes("".join(seq)):
            prot = g.translate().rstrip("*")
            w.writerow([pid, g.begin, g.end, g.strand,
                        int(g.partial_begin or g.partial_end), prot])
    for line in open(snakemake.input[0]):
        if line[0] == ">":
            flush(); pid, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    flush()
