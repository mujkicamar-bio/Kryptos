"""S8d: structural homology for EVERY dark protein, via Foldseek.

Sequence search has already failed on these proteins by definition. Structure reaches
further back in evolutionary time, so a family with no sequence homolog can still have a
recognisable fold - and a fold implies a mechanism, which is a testable hypothesis.

TWO OUTCOMES, BOTH INFORMATIVE
  a significant match to a characterised fold  -> a functional hypothesis
  a confident structure with NO match          -> a candidate novel fold, the highest-risk
                                                  and highest-reward stratum

WHY ProstT5 RATHER THAN FOLDING EVERY FAMILY
Foldseek can derive the 3Di structural alphabet directly from sequence with its ProstT5
model, which avoids running ColabFold over hundreds of thousands of families - days of GPU
time - before knowing which are worth folding. Proteins that Foldseek flags as promising
are then worth folding properly for confirmation.

This is a screening step, not a structure determination.

WHY EVERY DARK PROTEIN AND NOT ONE REPRESENTATIVE PER FAMILY

The input is results/s6/dark_proteins.faa, the whole dark set. A cluster representative is
chosen by MMseqs2 on sequence criteria that have nothing to do with which member is most
structurally informative, and at 30% identity - FESNov's deep-homology setting - members of
one family can differ enough that only some of them reach a recognisable fold. Searching
only representatives would therefore lose real structural evidence for reasons unrelated to
structure. `folds` is one of the four reality tests and the only route to a fourth line of
evidence, so losing it silently costs candidates their eligibility.

WHY THE DESCRIPTION IS CARRIED, NOT JUST THE ACCESSION

Foldseek's `target` is a PDB accession such as `12as-assembly1_A`, which names nothing. The
`theader` field carries the real description ("... ASPARAGINE SYNTHETASE ..."), and two
things downstream need it: the nucleic_acid_binding stratum, which is assigned from what
the fold IS, and the hypothesis written into the synthesis order, which a bench scientist
has to be able to read.
"""
import _ctx  # noqa: F401
import csv
import pathlib
import subprocess

cfg = snakemake.params.structure
target_db = snakemake.params.target_db
tmp = pathlib.Path(snakemake.output[0]).parent / "foldseek_tmp"

subprocess.run(
    f"foldseek easy-search {snakemake.input.faa} {target_db} "
    f"{snakemake.output[0]}.raw {tmp} --prostt5-model {snakemake.params.prostt5} "
    f"-e {cfg['max_evalue']} --threads {snakemake.threads} "
    # `prob` is a valid output field but is derived from Calpha coordinates, which a
    # ProstT5 query database does not carry - requesting it makes foldseek exit 1 and the
    # whole stage returns nothing. Verified in review: dropping this one field gives
    # exit=0 with rows returned.
    #
    # `theader` is the target's description line. Without it the only thing recorded about
    # a structural match is its PDB accession, which cannot be read by a person and cannot
    # be classified into a stratum.
    f"--format-output query,target,theader,fident,alnlen,evalue,bits "
    f"--max-seqs 10 -v 1",
    shell=True, check=True)

rows = []
raw = pathlib.Path(f"{snakemake.output[0]}.raw")
if raw.exists():
    best = {}
    for line in open(raw):
        q, t, theader, fident, alnlen, ev, bits = line.rstrip("\n").split("\t")
        if q not in best or float(ev) < float(best[q]["evalue"]):
            best[q] = {"seq_id": q, "target": t, "target_description": theader,
                       "fident": fident, "alnlen": alnlen, "evalue": ev, "bits": bits}
    rows = list(best.values())
else:
    # Reaching here means foldseek exited 0 but produced no file at all, which is a tool
    # contract violation rather than a legitimate "no hits" result.
    raise SystemExit(
        f"foldseek reported success but wrote no output to {raw}. Structure evidence would "
        "be silently empty, which would quietly empty the novel_fold stratum.")

with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=["seq_id", "target", "target_description",
                                        "fident", "alnlen", "evalue", "bits"],
                       delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"structural matches={len(rows)}")
