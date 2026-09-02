#!/usr/bin/env python3
"""Phase A.1-A.3 of PLASMIDSCOPE_PRIMARY_PLAN.md: characterize PlasmidScope's deduplicated `ALL`
master, apply the completeness gate with a full ledger, and reconstruct multi-source provenance for
the complete set. Data-only, stdlib only, no tools.

PlasmidScope's `ALL` table IS the non-redundant curated database: cross-source (and IMG/PR-internal)
duplicates were collapsed with MMseqs2 (--cov-mode 0 -c 1.0 --min-seq-id 1.0, i.e. 100% identity &
coverage), each surviving plasmid tagged in `Data_Source` with the comma-separated set of source DBs
it was consolidated from. We therefore never see cross-dataset duplicates.

Inputs:  data/plasmidscope_primary/all_metadata.tsv   (ALL.plasmid_list, downloaded)
Outputs:
  data/plasmidscope_primary/complete_provenance.tsv   (the deduplicated COMPLETE set + provenance)
  reports/plasmidscope_primary_characterization.md
Run:  python3 scripts/characterize_plasmidscope_all.py
"""
import csv
import os
from collections import Counter

ALL = "data/plasmidscope_primary/all_metadata.tsv"
OUT = "data/plasmidscope_primary/complete_provenance.tsv"
REPORT = "reports/plasmidscope_primary_characterization.md"

META_SOURCES = {"IMG-PR", "mMGE"}                                    # metagenome-derived
ISO_SOURCES = {"GenBank", "RefSeq", "PLSDB", "COMPASS", "ENA", "DDBJ", "Kraken2", "TPA"}  # isolate/INSDC


def norm(s):
    return "PLSDB" if s == "PLDSB" else s  # PlasmidScope misspells PLSDB as PLDSB throughout


def lifestyle_of(srcset):
    has_m, has_i = bool(srcset & META_SOURCES), bool(srcset & ISO_SOURCES)
    if has_m and has_i:
        return "mixed"
    return "metagenomic" if has_m else "isolate" if has_i else "other"


def main():
    completeness = Counter()
    member = Counter()
    lifestyle = Counter()
    nsrc = Counter()
    n_complete = raw_entries = 0

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(ALL, newline="") as fh, open(OUT, "w", newline="") as out:
        r = csv.DictReader(fh, delimiter="\t")
        w = csv.writer(out, delimiter="\t")
        w.writerow(["plasmid_id", "sources", "n_sources", "n_raw_entries", "lifestyle",
                    "topology", "size_bp", "host", "predicted_mobility"])
        for row in r:
            completeness[row["Completeness"] or "-"] += 1
            if row["Completeness"] != "complete":
                continue
            n_complete += 1
            parts = [norm(x) for x in row["Data_Source"].split(",") if x]
            raw_entries += len(parts)
            srcset = set(parts)
            nsrc[len(srcset)] += 1
            for db in srcset:
                member[db] += 1
            life = lifestyle_of(srcset)
            lifestyle[life] += 1
            w.writerow([row["Plasmid_ID"], ";".join(sorted(srcset)), len(srcset), len(parts),
                        life, row["Topology"], row["Size (bp)"], row["Host"],
                        row["Predicted_Mobility"]])

    total = sum(completeness.values())
    L, w = [], lambda s: L.append(s)
    w("# PlasmidScope-primary — data characterization (Phase A.1-A.3)\n")
    w("**Pipeline:** `scripts/characterize_plasmidscope_all.py` reading "
      "`data/plasmidscope_primary/all_metadata.tsv` (PlasmidScope's deduplicated `ALL` table). "
      "Regenerate with:\n```\npython3 scripts/characterize_plasmidscope_all.py\n```\n")
    w("**Deduplication (PlasmidScope's own, quoted from the paper):** duplicates detected with "
      "MMseqs2 v15.6f452 at 100% identity & 100% coverage "
      "(`--cov-mode 0 -c 1.0 --min-seq-id 1.0`); the surviving non-redundant plasmids are the "
      "`ALL` table, each row tagged with the set of source DBs it was consolidated from. So this "
      "set contains **no duplicate plasmids across (or within) source databases**.\n")

    w("## Completeness gate (ledger)\n")
    w("| Completeness | n | disposition |")
    w("|---|---:|---|")
    w(f"| complete (closed) | {completeness.get('complete',0):,} | **KEEP** |")
    w(f"| incomplete | {completeness.get('incomplete',0):,} | discard (fragment/contig) |")
    w(f"| unknown ('-') | {completeness.get('-',0):,} | discard (completeness not asserted) |")
    w(f"| **total in ALL** | {total:,} | |\n")
    w(f"Working universe after the completeness gate: **{n_complete:,} complete plasmids** "
      f"({raw_entries:,} raw source-entries collapsed into them, avg "
      f"{raw_entries/n_complete:.2f}/plasmid — dedup consolidates, it does not drop biology).\n")

    w("## Source-database membership (complete set; a plasmid can belong to several)\n")
    w("| source | member plasmids | type |")
    w("|---|---:|---|")
    for db, c in member.most_common():
        w(f"| {db} | {c:,} | {'metagenome' if db in META_SOURCES else 'isolate'} |")
    w("")
    w("## Lifestyle composition (complete set)\n")
    w("| lifestyle | n | % |")
    w("|---|---:|---:|")
    for k in ("metagenomic", "isolate", "mixed", "other"):
        if lifestyle.get(k):
            w(f"| {k} | {lifestyle[k]:,} | {100*lifestyle[k]/n_complete:.1f}% |")
    w("")
    w("## Source multiplicity (how many DBs each complete plasmid came from)\n")
    w("| distinct source DBs | plasmids |")
    w("|---:|---:|")
    for k in sorted(nsrc):
        w(f"| {k} | {nsrc[k]:,} |")
    multi = sum(v for k, v in nsrc.items() if k >= 2)
    w(f"\n{multi:,} complete plasmids are shared across >=2 source databases.\n")
    w("## Next in Phase A\n")
    w("- A.4 lab-made/synthetic flagging (conservative, logged).\n")
    w("- A.5 raw environment recovery per plasmid (IMG/PR GOLD ecosystem + NCBI BioSample); blanks "
      "kept as `unknown`.\n")

    with open(REPORT, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"wrote {OUT} ({n_complete:,} complete plasmids) and {REPORT}")


if __name__ == "__main__":
    main()
