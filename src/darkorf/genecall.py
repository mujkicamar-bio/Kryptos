"""Calling the genes of one plasmid, as a unit of work for a process pool.

Pyrodigal is single-threaded and meta mode calls every plasmid alone, so S1 deals the
records over a pool. The worker lives here rather than in the script because a pool
pickles its function by import path, and a function defined inside a Snakemake script has
none.
"""
import pyrodigal

from darkorf.circular import is_circular, overlap_for, resolve_origin_genes

_finder = None
_min_gene_nt = None


def configure(min_gene_nt):
    """Set the gene caller's floor for this process. Called once per pool worker."""
    global _finder, _min_gene_nt
    _min_gene_nt = min_gene_nt
    _finder = pyrodigal.GeneFinder(
        meta=True,                              # no training set: every plasmid is called alone
        min_gene=min_gene_nt,
        min_edge_gene=min(min_gene_nt, 60),     # edge genes on genuinely linear molecules
    )


def call_genes(plasmid_id, sequence, topology):
    """Call one record: (plasmid_id, genes, redundant_dropped).

    Topology decides whether the record gets origin repair. Unknown topology is treated as
    linear: extending a genuinely linear molecule would fabricate a junction that does not
    exist and could invent a chimeric gene across the two ends.
    """
    length = len(sequence)
    if not length:
        return plasmid_id, [], 0

    circular = is_circular(topology)
    if circular:
        # Append the head to the tail so a gene straddling the cut becomes contiguous.
        search_seq = sequence + sequence[:overlap_for(length)]
    else:
        search_seq = sequence

    raw = [{"start": g.begin, "end": g.end, "strand": g.strand,
            "partial": int(g.partial_begin or g.partial_end),
            # OPEN ISSUE (spec section 8, translation table): meta mode picks a model per
            # plasmid and some models use table 4 (TGA = Trp). Recorded so the affected
            # genes can be found; how to treat them is not decided yet.
            "translation_table": g.translation_table,
            "seq": g.translate().rstrip("*")}
           for g in _finder.find_genes(search_seq)]

    if circular:
        genes = resolve_origin_genes(raw, original_length=length)
        return plasmid_id, genes, len(raw) - len(genes)
    return plasmid_id, [dict(g, origin_spanning=False) for g in raw], 0


def call_record(args):
    """Pool entry point: ((plasmid_id, sequence), topology) -> call_genes(...)."""
    (plasmid_id, sequence), topology = args
    return call_genes(plasmid_id, sequence, topology)
