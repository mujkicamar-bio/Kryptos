#!/usr/bin/env python3
"""Phase A.5 (local part + fetch-list assembly). For every deduplicated complete plasmid, parse its
per-source original IDs (the comma-joined `Plasmid_ID`), recover the IMG/PR GOLD ecosystem locally,
and assemble the exact list of NCBI accessions / SRA runs still needing an online fetch (reusing the
env_raw.tsv we already have so we don't refetch). Data-only; the GOLD part needs no network.

Environment channels (kept RAW, all labels retained per decision #3; blanks -> unknown, never dropped):
  IMG-PR                        -> GOLD `ecosystem` string (local IMGPR_plasmid_data.tsv)
  RefSeq/GenBank/PLSDB/COMPASS/
  ENA(EMBL)/DDBJ/Kraken2/TPA    -> INSDC nuccore accession -> NCBI BioSample/source (fetched later)
  mMGE                          -> SRA run accession        -> NCBI BioSample (fetched later)

Outputs (data/plasmidscope_primary/):
  env_local.tsv                 rep_id, sources, imgpr_ecosystem_raw, insdc_accessions, sra_runs
  fetch_insdc_accessions.txt    unique INSDC accessions still needing a fetch (env_raw reused)
  fetch_sra_runs.txt            unique SRA runs (mMGE) needing a fetch
Run:  python3 scripts/recover_env_local.py
"""
import csv
import os
import re
import sys

PROV = "data/plasmidscope_primary/complete_provenance.tsv"
IMGPR = "data/Cus_PR/IMG_VR_2023-08-08_1/IMGPR_plasmid_data.tsv"
ENV_RAW = "data/processed/plasmidscope/biosample/env_raw.tsv"
OUTDIR = "data/plasmidscope_primary"
csv.field_size_limit(sys.maxsize)

INSDC_PREFIXES = ["RefSeq_", "GenBank_", "PLSDB_", "PLDSB_", "COMPASS_", "EMBL_", "DDBJ_",
                  "Kraken2_", "TPA_"]


def strip_ver(acc):
    return acc.rsplit(".", 1)[0]


def main():
    # IMG/PR GOLD ecosystem lookup (id -> raw ecosystem string), non-blank only
    gold = {}
    with open(IMGPR, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            eco = (r.get("ecosystem") or "").strip()
            if eco and eco.lower() != "unclassified":
                gold[r["plasmid_id"]] = eco
    print(f"GOLD ecosystem lookup: {len(gold):,} IMG/PR plasmids with a label", flush=True)

    # accessions already fetched (reuse) -> version-stripped set
    done = set()
    if os.path.exists(ENV_RAW):
        with open(ENV_RAW, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                done.add(strip_ver(r["accession"]))

    insdc_needed, sra_needed = set(), set()
    n_gold = n_any_insdc = n_any_sra = 0
    os.makedirs(OUTDIR, exist_ok=True)
    with open(PROV, newline="") as fh, open(os.path.join(OUTDIR, "env_local.tsv"), "w", newline="") as out:
        r = csv.DictReader(fh, delimiter="\t")
        w = csv.writer(out, delimiter="\t")
        w.writerow(["rep_id", "sources", "imgpr_ecosystem_raw", "insdc_accessions", "sra_runs"])
        for row in r:
            ids = row["plasmid_id"].split(",")
            rep = ids[0]
            ecos, insdc, sra = [], [], []
            for pid in ids:
                if pid.startswith("IMGPR_plasmid_"):
                    e = gold.get(pid)
                    if e and e not in ecos:
                        ecos.append(e)
                elif pid.startswith("mMGEs_"):
                    parts = pid.split("_")
                    # only real SRA run accessions are fetchable; many mMGE ids instead encode an
                    # assembly-contig prefix (mMGEs_k119_..., mMGEs_k141_...) with NO run -> env not
                    # recoverable from the id (would need mMGE's own contig->sample table).
                    if len(parts) > 1 and re.match(r"^[SDE]RR\d+$", parts[1]):
                        sra.append(parts[1])
                else:
                    for p in INSDC_PREFIXES:
                        if pid.startswith(p):
                            insdc.append(pid[len(p):])  # bare accession (with version)
                            break
            insdc = sorted(set(insdc)); sra = sorted(set(sra))
            if ecos:
                n_gold += 1
            if insdc:
                n_any_insdc += 1
            if sra:
                n_any_sra += 1
            for a in insdc:
                if strip_ver(a) not in done:
                    insdc_needed.add(a)
            sra_needed.update(sra)
            w.writerow([rep, row["sources"], " | ".join(ecos), ";".join(insdc), ";".join(sra)])

    with open(os.path.join(OUTDIR, "fetch_insdc_accessions.txt"), "w") as fh:
        fh.write("\n".join(sorted(insdc_needed)) + ("\n" if insdc_needed else ""))
    with open(os.path.join(OUTDIR, "fetch_sra_runs.txt"), "w") as fh:
        fh.write("\n".join(sorted(sra_needed)) + ("\n" if sra_needed else ""))

    print(f"plasmids with a GOLD ecosystem (local):   {n_gold:,}")
    print(f"plasmids with >=1 INSDC accession:        {n_any_insdc:,}")
    print(f"plasmids with >=1 mMGE SRA run:           {n_any_sra:,}")
    print(f"INSDC accessions still to fetch (gap):    {len(insdc_needed):,}  (env_raw reused: {len(done):,})")
    print(f"SRA runs to fetch (mMGE):                 {len(sra_needed):,}")


if __name__ == "__main__":
    main()
