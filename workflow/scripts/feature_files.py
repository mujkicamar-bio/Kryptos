"""S4: the feature files - GFF3 and GenBank.

The design names these as S4 deliverables alongside the TSV, and records that v1 declared
them and never wrote them. A TSV is what this project reasons over; a feature file is what
every genome browser, every downstream tool and every collaborator expects.

WHAT MAKES THIS MORE THAN FORMATTING

160,375 ORFs were reconstructed by S1 across the origin of a circular plasmid, and they
carry start > end. GFF3 forbids that outright and GenBank has dedicated syntax for it, so
the same gene has to be written two different ways. Get it wrong and the file still loads -
it just puts the gene somewhere else on the molecule. plasmidann.features owns both
conventions and is unit-tested against them.

NOTHING IS FILTERED HERE

Every ORF in the annotation table becomes a feature, including artefact-flagged ones. The
flags travel as attributes so a reader can act on them; the record stays complete.

ONE PASS OVER EACH INPUT

The annotation table is held indexed by plasmid, and the FASTA is streamed once. Reading
the annotation per plasmid instead would be a scan of a 9.3M-row file per record, which is
the quadratic pattern this pipeline has already had to fix twice.

THE SCOPE IS THE ANALYSIS SET

This stage used to stream the whole working-set FASTA and write a record for every sequence
in it. That is the corpus, not the run: on the 100-plasmid test configuration it produced
208,245 GenBank records and 11.9 GB. The TSV beside it was correctly scoped, so every count
a reader would think to check looked right. It now reads the analysis-set FASTA from S0.
"""
import collections
import csv

import _ctx  # noqa: F401

from darkorf.circular import is_circular
from plasmidann.fasta import iter_fasta
from plasmidann.features import genbank_location, gff3_features

# ------------------------------------------------------------------------------------
# Topology, from the master table. Length decides where a join() wraps, so it comes from
# the record itself (write_record), never from the table or the coordinates.
# ------------------------------------------------------------------------------------
topology = {}
with open(snakemake.input.master, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        topology[r["plasmid_id"]] = r.get("topology", "") or "linear"

# ------------------------------------------------------------------------------------
# Every ORF, grouped by plasmid.
# ------------------------------------------------------------------------------------
genes = collections.defaultdict(list)
with open(snakemake.input.annotation, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        genes[r["plasmid_id"]].append(r)
for rows in genes.values():
    rows.sort(key=lambda r: (int(r["start"]), int(r["end"])))


def attributes_for(r):
    """GFF3 attributes for one ORF.

    An empty value is omitted rather than written blank: a dark ORF has no product, and
    `product=` asserts that it has one which happens to be empty. Writing
    "hypothetical protein" would be worse still - it fabricates an annotation the cascade
    did not make, in the file most likely to be read by someone who did not run it.
    """
    attrs = {"ID": r["orf_id"]}
    if r.get("annot_label"):
        attrs["product"] = r["annot_label"]
    for key, column in (("functional_class", "functional_class"),
                        ("annot_tier", "annot_tier"),
                        ("dark_evidence", "dark_evidence")):
        if r.get(column):
            attrs[key] = r[column]
    if r.get("partial") == "1":
        attrs["partial"] = "true"
    if r.get("spans_origin") == "1":
        attrs["spans_origin"] = "true"
    if r.get("artefact_flag") == "1":
        attrs["artefact"] = r.get("artefact_reason") or "true"
    return attrs


def wrap(seq, width=60):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


n_features = n_records = 0
with open(snakemake.output.gff3, "w") as gff, open(snakemake.output.genbank, "w") as gbk:
    gff.write("##gff-version 3\n")

    def write_record(name, seq):
        """Write one plasmid's GFF3 and GenBank records."""
        global n_features, n_records
        # The record's own length, as S0 wrote it: S1 wrapped the coordinates on this
        # sequence, which is shorter than size_bp wherever S0 removed a terminal repeat.
        L = len(seq)
        rows = genes.get(name, [])
        n_records += 1

        gff.write(f"##sequence-region {name} 1 {L}\n")
        for r in rows:
            gene = {"plasmid_id": name, "orf_id": r["orf_id"], "start": int(r["start"]),
                    "end": int(r["end"]), "strand": r["strand"]}
            for feature in gff3_features(gene, length=L, attributes=attributes_for(r)):
                gff.write("\t".join(str(x) for x in feature) + "\n")
                n_features += 1

        circular = "circular" if is_circular(topology.get(name)) else "linear"
        gbk.write(f"LOCUS       {name:<20}{L} bp    DNA     {circular}  UNK\n")
        gbk.write(f"DEFINITION  {name} annotated by plasmidann.\n")
        gbk.write("FEATURES             Location/Qualifiers\n")
        gbk.write(f"     source          1..{L}\n")
        gbk.write('                     /mol_type="genomic DNA"\n')
        for r in rows:
            loc = genbank_location(int(r["start"]), int(r["end"]), r["strand"], L)
            gbk.write(f"     CDS             {loc}\n")
            gbk.write(f'                     /locus_tag="{r["orf_id"]}"\n')
            if r.get("annot_label"):
                gbk.write(f'                     /product="{r["annot_label"]}"\n')
            gbk.write(f'                     /note="functional_class='
                      f'{r.get("functional_class", "")}"\n')
        gbk.write("ORIGIN\n")
        for i, line in enumerate(wrap(seq.lower())):
            blocks = " ".join(line[j:j + 10] for j in range(0, len(line), 10))
            gbk.write(f"{i * 60 + 1:>9} {blocks}\n")
        gbk.write("//\n")

    for name, sequence in iter_fasta([snakemake.input.fasta]):
        write_record(name, sequence)

print(f"feature files: {n_records} records, {n_features} GFF3 features")
