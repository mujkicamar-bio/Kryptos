"""DefenseFinder's two phases, each on the representation it needs.

Phase 1 is an HMM search with the profiles of the DefenseFinder, RM and CasFinder models
that asks per protein whether it looks like a defence component; it does not depend on gene
order, so it runs once on the dereplicated proteins. Phase 2 is MacSyFinder system calling
with the 554 models of the same three families, whose quorum and co-localisation rules
(e.g. inter_gene_max_space="3") count intervening genes, so it runs on every ORF of a
plasmid in genomic order. Between the two, component hits are propagated to every ORF
sharing the sequence, and plasmids without a component are left out of phase 2,
since they cannot satisfy any model's quorum.

Phase 2 and CONJScan read MacSyFinder's gembase format: one FASTA holding many replicons,
each ORF named <replicon>_<position>.
"""
import itertools


def propagate_components(component_seqs, orf_to_seq):
    """The ORFs whose protein is a component hit: an identical sequence has the same hits.

    `component_seqs` is the set of seq_ids with a phase 1 hit; `orf_to_seq` maps
    orf_id -> seq_id.
    """
    return {orf for orf, seq in orf_to_seq.items() if seq in component_seqs}


def candidate_plasmids(component_orfs):
    """The plasmids carrying at least one component ORF (orf_id is <plasmid_id>|<n>)."""
    return {orf.rsplit("|", 1)[0] for orf in component_orfs}


def order_orfs(orfs):
    """One plasmid's ORFs in genomic order.

    An origin-spanning gene runs start..length then 1..end, so its start is greater than its
    end. It is placed first, so that positions follow the molecule from its origin.
    """
    def key(o):
        spans = str(o.get("spans_origin", "0")) == "1"
        return (0 if spans else 1, int(o["start"]))
    return sorted(orfs, key=key)


def gembase_id(plasmid_id, position):
    """A MacSyFinder gembase identifier: <replicon>_<zero-padded position>.

    MacSyFinder takes the replicon to be everything before the LAST underscore, so the
    underscores of the plasmid id (COMPASS_AB007909.1) become hyphens. The position is
    zero-padded so that lexical order matches genomic order.
    """
    return f"{plasmid_id.replace('_', '-')}_{position:05d}"


def gembase_records(index_rows, keep=None):
    """(gembase_id, plasmid_id, row) for every ORF of an ORF index, in genomic order.

    `index_rows` are orf_index.tsv rows, which are sorted by plasmid; they are read one
    plasmid at a time, so only one plasmid is held in memory. `keep`, when given, is the set
    of plasmids to include. Raises ValueError when a replicon name appears twice: either the
    rows are not grouped by plasmid, or two plasmid ids differ only in '_' versus '-', and
    either would merge ORFs of different molecules in one replicon.
    """
    seen = set()
    for plasmid, rows in itertools.groupby(index_rows, key=lambda r: r["plasmid_id"]):
        if keep is not None and plasmid not in keep:
            continue
        replicon = plasmid.replace("_", "-")
        if replicon in seen:
            raise ValueError(f"gembase replicon {replicon} (plasmid {plasmid}) occurs twice: "
                             "the ORF index is not grouped by plasmid, or two plasmid ids "
                             "map to the same replicon name")
        seen.add(replicon)
        for position, orf in enumerate(order_orfs(list(rows)), start=1):
            yield gembase_id(plasmid, position), plasmid, orf


def parse_all_systems(paths):
    """Component hits from MacSyFinder's all_systems.tsv files (phase 1).

    MacSyFinder writes no best_solution.tsv under `--db-type unordered`; all_systems.tsv is
    its own output and holds the hit (hit_id), the component (gene_name), the model
    (model_fqn) and the independent e-value (hit_i_eval). The file opens with '#' lines and
    a blank line before its header, and a model family that matched nothing writes comments
    alone.
    """
    rows = []
    for path in paths:
        with open(path, newline="") as fh:
            lines = [l for l in fh if l.strip() and not l.startswith("#")]
        if not lines:
            continue
        header = lines[0].rstrip("\n").split("\t")
        for line in lines[1:]:
            fields = dict(zip(header, line.rstrip("\n").split("\t")))
            seq_id = fields.get("hit_id", "")
            if not seq_id:
                continue
            rows.append({
                "seq_id": seq_id,
                "component": fields.get("gene_name", ""),
                "model": fields.get("model_fqn", ""),
                "hit_evalue": fields.get("hit_i_eval", ""),
            })
    return rows
