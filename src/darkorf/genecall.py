"""Calling the genes of one plasmid, as a unit of work for a process pool.

Pyrodigal is single-threaded and meta mode calls every plasmid alone, so rule orf_call deals
the records over a pool. The worker lives here rather than in the script because a pool
pickles its function by import path, and a function defined inside a Snakemake script has
none.
"""
import pyrodigal

from darkorf.circular import is_circular, overlap_for, resolve_origin_genes

_finder = None


def configure(min_gene_nt):
    """Set the gene caller's floor for this process. Called once per pool worker."""
    global _finder
    _finder = pyrodigal.GeneFinder(
        meta=True,                              # no training set: every plasmid is called alone
        min_gene=min_gene_nt,
        min_edge_gene=min(min_gene_nt, 60),     # Prodigal's default floor for edge genes
    )


def call_genes(plasmid_id, sequence, topology):
    """Call one record; circular topologies get origin repair.

    Returns (plasmid_id, genes, redundant_dropped).
    """
    length = len(sequence)
    if not length:
        return plasmid_id, [], 0

    circular = is_circular(topology)
    search_seq = sequence + sequence[:overlap_for(length)] if circular else sequence

    raw = [{"start": g.begin, "end": g.end, "strand": g.strand,
            # The call runs off the left (begin) or right (end) edge of the sequence.
            "partial": int(g.partial_begin or g.partial_end),
            "partial_begin": int(g.partial_begin), "partial_end": int(g.partial_end),
            # Meta mode picks one model per call, and some models use table 4 (TGA = Trp).
            "translation_table": g.translation_table,
            "seq": g.translate().rstrip("*")}
           for g in _finder.find_genes(search_seq)]

    if circular:
        genes = resolve_origin_genes(raw, length, len(search_seq))
        return plasmid_id, genes, len(raw) - len(genes)
    return plasmid_id, [dict(g, origin_spanning=False) for g in raw], 0


def call_record(args):
    """Pool entry point: ((plasmid_id, sequence), topology) -> call_genes(...)."""
    (plasmid_id, sequence), topology = args
    return call_genes(plasmid_id, sequence, topology)
