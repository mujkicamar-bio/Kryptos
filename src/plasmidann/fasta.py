"""Streaming a plasmid FASTA, plain or gzipped."""
import gzip


def iter_fasta(paths):
    """Yield (plasmid_id, sequence) for every record across `paths`, in file order.

    Plain or gzipped, decided per file, because the working set is delivered compressed
    and a smoke set is written plain. Handling only one would make the other yield no
    sequence at all - a silent skip.

    The identifier is the header up to the first whitespace. Everything after it is
    description, and plasmid_id is the join key for every table in the run.
    """
    for path in paths:
        opener = gzip.open if str(path).endswith(".gz") else open
        name, chunks = None, []
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if name is not None:
                        yield name, "".join(chunks)
                    name, chunks = line[1:].split()[0], []
                else:
                    chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks)


def split_fasta(path, n, outdir):
    """Deal the records of one FASTA into at most `n` files of similar total length.

    For tools that walk replicons one at a time (IntegronFinder, ISEScan): n single-thread
    runs over n chunks use n cores, where one run with n threads kept one busy. Streaming,
    each record goes to the chunk with the fewest bases so far, so nothing is held in
    memory and no chunk is much longer than the others. Returns the chunk paths that
    received at least one record, in order.
    """
    import pathlib
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = [outdir / f"chunk_{i:03d}.fna" for i in range(n)]
    handles = [open(p, "w") for p in paths]
    bases = [0] * n
    try:
        for name, seq in iter_fasta([path]):
            i = bases.index(min(bases))
            handles[i].write(f">{name}\n{seq}\n")
            bases[i] += len(seq)
    finally:
        for h in handles:
            h.close()
    return [p for p, b in zip(paths, bases) if b]
