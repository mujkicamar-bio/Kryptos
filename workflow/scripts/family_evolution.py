"""S7b: evolutionary evidence per dark family.

Everything upstream reports what is ABSENT - no Pfam hit, no Swiss-Prot hit, no nr hit.
Absence is exactly what a spurious ORF also produces, which is why FESNov discarded 94.3%
of its 7,052,473 novel clusters.

This stage reports something POSITIVE:

  dN/dS < 0.5        selection is suppressing replacement changes, so there IS a protein
  RNAcode P < 0.05   the region carries coding signal independent of the gene caller
  lineage breadth    it is not a single-lineage accident

WHY RNAcode SCORES BOTH STRANDS HERE

RNAcode evaluates all six reading frames and reports the sense and antisense signal
separately. For this project that is the point, not a detail. A shadow ORF is the reverse
complement of a real gene, so the coding signal on its ANTISENSE strand should be the
stronger of the two - and a shadow ORF is exactly the artefact class that survives every
absence-based test in this pipeline. Recording only the sense P would discard the one
number that separates them.

The two are reported, never combined into a verdict: `rnacode_p`, `rnacode_p_antisense`,
and `coding_signal` as the declared threshold applied to the sense P.

dN/dS is also the only filter that catches the artefact class nothing else does. A shadow
ORF on the reverse-complement strand of a real gene is conserved, multi-species, and passes
every absence-based test - but the selection acting on that locus is acting on the gene on
the OTHER strand.

The consensus re-check matters because a protein can miss every per-sequence threshold
while its family is collectively recognisable; this removed 6.5% of clusters in
Pavlopoulos et al.
"""
import _ctx  # noqa: F401
import collections
import csv
import pathlib
import statistics
import subprocess

from plasmidann import scratch
from plasmidann.evolution import dnds_detail, back_translate, consensus

cfg = snakemake.params.evolution

def read_fasta(path):
    out, name, buf = {}, None, []
    for line in open(path):
        if line[0] == ">":
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out

proteins = read_fasta(snakemake.input.faa)
cds = read_fasta(snakemake.input.cds)

families = []
with open(snakemake.input.families, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        families.append(r)

tmpdir = str(scratch.scratch_dir(pathlib.Path(snakemake.output[0]).parent))
# The consensus carries the family's shared signal and is re-searched at S7c.
consensus_out = open(snakemake.output.consensus, "w")
cols = ["family_id", "n_aligned", "dnds_median", "dnds_min", "n_pairs",
        "under_purifying_selection", "dnds_status", "rnacode_p", "rnacode_p_antisense",
        "rnacode_status", "coding_signal", "evidence_note"]


def clustal(alignment, path):
    """Write a codon alignment as Clustal W.

    RNAcode reads Clustal W or MAF and nothing else - handed FASTA it prints
    "ERROR: Unknown alignment file format" and EXITS 0, so the caller has to look at what
    it wrote rather than at its status.
    """
    names = list(alignment)
    width = 60
    with open(path, "w") as fh:
        fh.write("CLUSTAL W (1.81) multiple sequence alignment\n\n\n")
        length = len(alignment[names[0]])
        for i in range(0, length, width):
            for n in names:
                fh.write(f"{n[:15]:<16}{alignment[n][i:i + width]}\n")
            fh.write("\n\n")


def rnacode(alignment, path, threads):
    """Best sense and antisense P from RNAcode over a codon alignment.

    Returns (p_sense, p_antisense, status). Tabular columns are
    HSS, strand, frame, length, from, to, name, start, end, score, P.

    A tool that reports failure on stdout and exits 0 cannot be trusted to its return code,
    so an unparsable output is NO_OUTPUT rather than "no coding signal". The distinction
    matters: no signal is evidence against a family, and a tool that did not run is not.
    """
    clustal(alignment, path)
    proc = subprocess.run(f"RNAcode -t {path}", shell=True, capture_output=True, text=True)
    best = {"+": None, "-": None}
    for line in proc.stdout.splitlines():
        f = line.split()
        if len(f) < 11 or f[1] not in best:
            continue
        try:
            p = float(f[10])
        except ValueError:
            continue
        if best[f[1]] is None or p < best[f[1]]:
            best[f[1]] = p
    if best["+"] is None and best["-"] is None:
        note = proc.stdout.strip().splitlines()
        return None, None, ("NO_OUTPUT" if any(l.startswith("ERROR") for l in note)
                            else "NO_SIGNAL")
    return best["+"], best["-"], "MEASURED"

n_tested = n_purifying = n_coding = 0
with open(snakemake.output.tsv, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()

    for fam in families:
        members = fam["members"].split(",")
        row = {"family_id": fam["family_id"], "n_aligned": 0, "dnds_median": "",
               "dnds_min": "", "n_pairs": 0, "under_purifying_selection": "",
               "dnds_status": "", "rnacode_p": "", "rnacode_p_antisense": "",
               "rnacode_status": "", "coding_signal": "", "evidence_note": ""}

        usable = [m for m in members if m in cds and m in proteins]
        if len(usable) < cfg["min_members_for_dnds"]:
            # Not an error: most families are small. Reported so the shortfall is visible
            # rather than looking like a failed test.
            row["evidence_note"] = "too_few_members"
            row["dnds_status"] = "TOO_FEW_MEMBERS"
            row["rnacode_status"] = "TOO_FEW_MEMBERS"
            w.writerow(row)
            continue

        # Alignment is superlinear and the marginal information from the 200th member is
        # negligible, so cap it. Members are taken in file order, which is deterministic.
        usable = usable[:cfg["max_members_aligned"]]
        row["n_aligned"] = len(usable)

        pf = f"{tmpdir}/{fam['family_id']}.faa"
        with open(pf, "w") as fh:
            for m in usable:
                fh.write(f">{m}\n{proteins[m]}\n")

        # Proteins are aligned, not nucleotides: protein alignment is far more reliable at
        # the 30% identities these families show. The alignment is then projected onto
        # codons, so every protein gap becomes exactly three nucleotide gaps.
        aln = subprocess.run(f"mafft --auto --quiet --thread 1 {pf}",
                             shell=True, capture_output=True, text=True)
        if aln.returncode != 0:
            row["evidence_note"] = "alignment_failed"
            row["dnds_status"] = "ALIGNMENT_FAILED"
            row["rnacode_status"] = "ALIGNMENT_FAILED"
            w.writerow(row)
            continue

        aligned = {}
        name, buf = None, []
        for line in aln.stdout.splitlines():
            if line.startswith(">"):
                if name:
                    aligned[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if name:
            aligned[name] = "".join(buf)

        # One consensus per family, from the PROTEIN alignment. Written for every family
        # that could be aligned at all, including those with no usable codon pairs: the
        # re-check asks whether the family is collectively recognisable, which is a
        # different question from whether its divergence can be measured.
        family_consensus = consensus(aligned)
        if family_consensus:
            consensus_out.write(f">{fam['family_id']}\n{family_consensus}\n")

        codon_aln = {}
        for m, ap in aligned.items():
            try:
                codon_aln[m] = back_translate(ap, cds[m])
            except ValueError:
                continue

        # Statuses are counted, not just values. A family of identical sequences
        # (NO_DIVERGENCE) is a different thing from one measured and found neutral, and S9
        # scores them differently - the first withholds judgement, the second is evidence
        # against. Collapsing both into a bare None penalised the most conserved families.
        ratios, statuses = [], collections.Counter()
        ids = sorted(codon_aln)
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                # min_codons comes from config and is applied HERE, per pair. It was
                # declared and never read: the only floor in force was the arithmetic
                # minimum of three, so eight-codon fragments produced dN/dS values that
                # then fired purifying_selection, the strongest of the four reality tests.
                r, status = dnds_detail(codon_aln[ids[i]], codon_aln[ids[j]],
                                        min_codons=cfg["min_codons"])
                statuses[status] += 1
                if r is not None and r != float("inf"):
                    ratios.append(r)

        # Coding potential, independent of the gene caller and of dN/dS. Both strands: for
        # a shadow ORF the antisense signal is expected to be the stronger one.
        if len(codon_aln) >= cfg["min_members_for_dnds"]:
            p_sense, p_anti, rc_status = rnacode(
                codon_aln, f"{tmpdir}/{fam['family_id']}.aln", snakemake.threads)
            row["rnacode_status"] = rc_status
            if p_sense is not None:
                row["rnacode_p"] = p_sense
                row["coding_signal"] = int(p_sense < cfg["rnacode_max_p"])
                n_coding += int(p_sense < cfg["rnacode_max_p"])
            if p_anti is not None:
                row["rnacode_p_antisense"] = p_anti
        else:
            row["rnacode_status"] = "TOO_FEW_MEMBERS"

        row["n_pairs"] = len(ratios)
        if ratios:
            n_tested += 1
            row["dnds_median"] = round(statistics.median(ratios), 4)
            row["dnds_min"] = round(min(ratios), 4)
            purifying = int(row["dnds_median"] < cfg["dnds_purifying_max"])
            row["under_purifying_selection"] = purifying
            n_purifying += purifying
            row["dnds_status"] = "MEASURED"
        else:
            row["evidence_note"] = "no_informative_pairs"
            # Which kind of absence? The commonest status across the pairs is the honest
            # summary, and it is what S9 reads.
            row["dnds_status"] = (statuses.most_common(1)[0][0] if statuses
                                  else "NO_INFORMATIVE_PAIRS")

        w.writerow(row)

consensus_out.close()

print(f"families={len(families)} with_dnds={n_tested} purifying={n_purifying} "
      f"coding_signal={n_coding}")

# The scratch directory is removed only here, on the ordinary path. A script that raised
# never reaches this line, and its intermediates are what the failure is diagnosed from.
scratch.release(tmpdir)

