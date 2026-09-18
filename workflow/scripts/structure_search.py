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

SCOPE: REPRESENTATIVES BY DEFAULT

Spec section 49 sets the discovery-scale strategy - "analyze dark-family representatives
where practical" - and section 79 makes it a success criterion: "structure is performed at
representative scale in the production run". So `structure.scope: representatives` is the
default and searches one sequence per dark family.

The cost is why. ProstT5 predicts the 3Di alphabet for every query, and it is a transformer:
on the full collection the difference between all dark proteins and one per family is the
difference between the largest job in the pipeline and a modest one.

`scope: all` searches every dark protein and is a real option, because a representative is
chosen by MMseqs2 on sequence criteria that have nothing to do with which member is most
structurally informative - at 30% identity, members of one family can differ enough that
only some reach a recognisable fold. Spec section 49 calls that "candidate scale" and puts
it downstream of discovery. Setting it here is supported and expensive; the default is what
the specification asks for.

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
import sys

from plasmidann import scratch

cfg = snakemake.params.structure
target_db = snakemake.params.target_db
tmp = scratch.scratch_dir(pathlib.Path(snakemake.output[0]).parent, "foldseek_tmp")

# ---- the query set: one sequence per dark family, or every dark protein ----------------
scope = cfg.get("scope", "representatives")
if scope not in ("representatives", "all"):
    sys.exit(f"structure.scope must be 'representatives' or 'all', not {scope!r}")

query_faa = snakemake.input.faa
if scope == "representatives":
    wanted = set()
    with open(snakemake.input.families, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("representative"):
                wanted.add(row["representative"])

    query_faa = str(pathlib.Path(snakemake.output[0]).parent / "structure_query.faa")
    n_written = 0
    with open(query_faa, "w") as out:
        emit = False
        for line in open(snakemake.input.faa):
            if line[0] == ">":
                emit = line[1:].split()[0] in wanted
                n_written += emit
            if emit:
                out.write(line)

    # A representative named in the families table but absent from the dark FASTA means the
    # two disagree about what the dark set is, and every structural count would inherit it.
    if n_written != len(wanted):
        sys.exit(f"S8d: {len(wanted)} family representatives declared but {n_written} found "
                 f"in {snakemake.input.faa} - the families table and the dark set disagree.")
    print(f"S8d: scope=representatives, searching {n_written} of "
          f"{sum(1 for l in open(snakemake.input.faa) if l[0] == '>')} dark proteins")
else:
    print("S8d: scope=all, searching every dark protein")

# GPU, when one is allocated.
#
# ProstT5 is a transformer, and predicting 3Di for every dark protein is the cost of this
# stage - not the Foldseek search that follows it. On CPU that is the difference between
# minutes and hours, so the flag is worth having; foldseek 10 takes --gpu on both createdb
# (the ProstT5 step) and search.
#
# It is CONFIGURED rather than detected. Auto-detecting a GPU would make the stage behave
# differently depending on which node it landed on, with nothing in the output saying
# which - and this pipeline records the settings that produced every row precisely so that
# a result can be traced. A run that asked for a GPU and did not get one should fail
# visibly, not silently take ten times longer.
gpu_flag = " --gpu 1" if cfg.get("gpu", False) else ""

subprocess.run(
    f"foldseek easy-search {query_faa} {target_db} "
    f"{snakemake.output[0]}.raw {tmp} --prostt5-model {snakemake.params.prostt5} "
    f"-e {cfg['max_evalue']} --threads {snakemake.threads}{gpu_flag} "
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

# The scratch directory is removed only here, on the ordinary path. A script that raised
# never reaches this line, and its intermediates are what the failure is diagnosed from.
scratch.release(tmp)

