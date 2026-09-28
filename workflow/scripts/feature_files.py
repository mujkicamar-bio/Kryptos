"""Write the annotation table as GFF3 and GenBank feature files.

Inputs: the per-ORF annotation table, the analysis-set FASTA (one record per sequence in
it, with or without ORFs) and the master table for each plasmid's topology. Every ORF
becomes a feature, artefact-flagged ones included; the flags travel as attributes.
Origin-spanning ORFs are encoded by plasmidann.features. The annotation table is held in
memory, indexed by plasmid, and the FASTA is streamed once.
"""
import collections
import csv
import datetime

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
        topology[r["plasmid_id"]] = r.get("topology", "")

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


def quoted(value):
    """A GenBank qualifier value in double quotes, with inner double quotes doubled."""
    return '"' + str(value).replace('"', '""') + '"'


def wrap(seq, width=60):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


# The LOCUS date is the date the file was written.
today = datetime.date.today().strftime("%d-%b-%Y").upper()
n_features = n_records = 0
with open(snakemake.output.gff3, "w") as gff, open(snakemake.output.genbank, "w") as gbk:
    gff.write("##gff-version 3\n")

    def write_record(name, seq):
        """Write one plasmid's GFF3 and GenBank records."""
        global n_features, n_records
        # The record's own length: ORF coordinates wrap on this sequence, which is shorter
        # than the master table's size_bp wherever the analysis set trimmed a terminal repeat.
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
        # INSDC column layout: name from column 13, length right-aligned in columns 30-40,
        # topology from column 56. A name longer than 16 characters shifts the rest right,
        # and the space after it keeps name and length apart for a parser.
        gbk.write(f"LOCUS       {name:<16} {L:>11} bp    DNA     {circular:<8} UNK "
                  f"{today}\n")
        gbk.write(f"DEFINITION  {name} annotated by plasmidann.\n")
        gbk.write("FEATURES             Location/Qualifiers\n")
        gbk.write(f"     source          1..{L}\n")
        gbk.write('                     /mol_type="genomic DNA"\n')
        for r in rows:
            loc = genbank_location(int(r["start"]), int(r["end"]), r["strand"], L,
                                   r["partial_begin"] == "1", r["partial_end"] == "1")
            gbk.write(f"     CDS             {loc}\n")
            gbk.write(f'                     /locus_tag={quoted(r["orf_id"])}\n')
            if r.get("annot_label"):
                gbk.write(f'                     /product={quoted(r["annot_label"])}\n')
            note = "functional_class=" + r.get("functional_class", "")
            gbk.write(f'                     /note={quoted(note)}\n')
        gbk.write("ORIGIN\n")
        for i, line in enumerate(wrap(seq.lower())):
            blocks = " ".join(line[j:j + 10] for j in range(0, len(line), 10))
            gbk.write(f"{i * 60 + 1:>9} {blocks}\n")
        gbk.write("//\n")

    for name, sequence in iter_fasta([snakemake.input.fasta]):
        write_record(name, sequence)

print(f"feature files: {n_records} records, {n_features} GFF3 features")
