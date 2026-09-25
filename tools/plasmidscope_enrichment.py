"""Reduce PlasmidScope's ALL.protein_list to the analysis set, keyed like our proteins.

PlasmidScope (Li et al., NAR 2025, 53:D179) ran eggNOG-mapper 2.1.12 on every protein.
The ORFs themselves are NOT uniform: for COMPASS, IMG-PR, PLSDB, mMGE and Kraken2 they are
Prodigal 2.6 calls with Prokka products, while for GenBank, RefSeq, DDBJ and EMBL they are
the submitters' deposited CDS and products (mostly NCBI PGAP). The column
`orf_source` keeps that distinction.

The protein table carries the sequence, so each row gets the identity our pipeline gives
a protein (plasmidann.dereplicate: sha256 of the sequence, first 32 hex characters) and
joins to our results without relying on coordinates.

ps_class is plasmidann.plasmidscope.ps_class (eggNOG fields only), the rule the pipeline
itself uses at S2p.
product_named is 1 when the product is anything but "hypothetical protein"; what that
means depends on orf_source.

Per unique protein, ps_class is the most informative class over its occurrences and
product_named is 1 if any occurrence is named.

Usage:
  python tools/plasmidscope_enrichment.py \
      data/PlasmidScope/annotation/ALL.protein_list.tsv.gz \
      data/plasmidscope_primary/analysis_set.tsv \
      data/PlasmidScope/annotation/analysis_set
Writes <prefix>_orfs.tsv.gz and <prefix>_proteins.tsv.gz, and prints a summary.
"""
import collections
import csv
import gzip
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from plasmidann.dereplicate import _seq_id  # noqa: E402  (the identity our proteins carry)
from plasmidann.plasmidscope import CLASSES, RANK, ps_class  # noqa: E402

ps_table, analysis_set, prefix = sys.argv[1:4]
csv.field_size_limit(sys.maxsize)

with open(analysis_set, newline="") as fh:
    wanted = {row["plasmid_id"] for row in csv.DictReader(fh, delimiter="\t")}


orf_cols = ["plasmid_id", "ps_protein_id", "orf_source", "start", "end", "strand", "seq_id",
            "ps_class", "product_named", "product", "cog_category", "cog_id", "kegg_ko",
            "pfams", "ec"]
proteins = {}                   # seq_id -> [length, ps_class, product_named, n_orfs, seq]
n_rows = n_kept = n_class_conflict = 0
orf_classes = collections.Counter()
plasmids_seen = set()

with gzip.open(ps_table, "rt", newline="") as fh, \
        gzip.open(f"{prefix}_orfs.tsv.gz", "wt", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(orf_cols)
    for r in csv.DictReader(fh, delimiter="\t"):
        n_rows += 1
        if r["Plasmid_ID"] not in wanted:
            continue
        n_kept += 1
        plasmids_seen.add(r["Plasmid_ID"])
        seq = r["Sequence"].rstrip("*")     # our ORFs are stored without the stop, too
        sid = _seq_id(seq)
        cls = ps_class(r)
        named = int(not r["Product"].startswith("hypothetical protein"))
        orf_classes[cls] += 1
        p = proteins.get(sid)
        if p is None:
            proteins[sid] = [len(seq), cls, named, 1, seq]
        else:
            n_class_conflict += p[1] != cls
            if RANK[cls] < RANK[p[1]]:
                p[1] = cls
            p[2] |= named
            p[3] += 1
        w.writerow([r["Plasmid_ID"], r["Protein_ID"], r["Orf Prediction Source"], r["Start"],
                    r["End"], r["Strand"], sid, cls, named, r["Product"], r["COG_category"],
                    r["COG_id"], r["KEGG_ko"], r["PFAMs"], r["EC_number"]])

with gzip.open(f"{prefix}_proteins.tsv.gz", "wt", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["seq_id", "length_aa", "ps_class", "product_named", "n_orfs", "seq"])
    for sid, (length, cls, named, n, seq) in proteins.items():
        w.writerow([sid, length, cls, named, n, seq])

prot_classes = collections.Counter(v[1] for v in proteins.values())
print(f"PlasmidScope rows read            {n_rows:,}")
print(f"rows on analysis-set plasmids     {n_kept:,}")
print(f"analysis-set plasmids covered     {len(plasmids_seen):,} / {len(wanted):,}")
print(f"unique proteins                   {len(proteins):,}")
print(f"identical sequence, other class   {n_class_conflict:,} occurrences")
for name, counts, total in (("ORFs", orf_classes, n_kept),
                            ("unique proteins", prot_classes, len(proteins))):
    print(f"ps_class over {name}:")
    for cls in CLASSES:
        print(f"  {cls:<17} {counts[cls]:>11,}  {counts[cls] / max(total, 1):6.1%}")
