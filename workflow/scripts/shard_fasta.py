import _ctx  # noqa: F401
import gzip

keep = set(open(snakemake.input.ids).read().split())
handles = [open(f, "w") for f in snakemake.output]

n, write = 0, False
with gzip.open(snakemake.input.fasta, "rt") as fh:
    for line in fh:
        if line[0] == ">":
            write = line[1:].split()[0] in keep
            if write:
                n += 1
        if write:
            handles[n % len(handles)].write(line)
for h in handles:
    h.close()
