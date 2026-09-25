"""IS elements from ISEScan (Xie & Tang, Bioinformatics 2017, 33:3340).

ISEScan finds insertion sequences de novo: profile HMMs of IS-family transposases locate
the transposase, and a search for terminal inverted repeats (TIRs) around it sets the
element's boundaries. It reports the IS family, the boundaries, the TIRs and whether the
element is complete ('c') or partial ('p': no TIR pair, or too short).

Only its `<fasta>.tsv` table is read. Columns used: seqID, family, cluster, isBegin,
isEnd, strand, E-value, type, tir. ISEScan leaves strand empty for some elements and
writes '-:-' as the TIR of an element without one; both are kept as absence ('').
"""
import csv

COLUMNS = ["plasmid_id", "is_id", "family", "cluster", "start", "end", "strand",
           "complete", "evalue", "tir"]


def parse_isescan(path):
    """One row per IS element, in file order; is_id is '<plasmid_id>|IS<n>' per plasmid."""
    rows = []
    n_by_plasmid = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pid = r["seqID"]
            n_by_plasmid[pid] = n_by_plasmid.get(pid, 0) + 1
            rows.append({
                "plasmid_id": pid,
                "is_id": f"{pid}|IS{n_by_plasmid[pid]}",
                "family": r["family"],
                "cluster": r["cluster"],
                "start": int(r["isBegin"]),
                "end": int(r["isEnd"]),
                "strand": r["strand"].strip(),
                "complete": int(r["type"] == "c"),
                "evalue": r["E-value"],
                "tir": "" if r["tir"] == "-:-" else r["tir"],
            })
    return rows
