# PlasmidScope-primary — data characterization (Phase A.1-A.3)

**Pipeline:** `scripts/characterize_plasmidscope_all.py` reading `data/plasmidscope_primary/all_metadata.tsv` (PlasmidScope's deduplicated `ALL` table). Regenerate with:
```
python3 scripts/characterize_plasmidscope_all.py
```

**Deduplication (PlasmidScope's own, quoted from the paper):** duplicates detected with MMseqs2 v15.6f452 at 100% identity & 100% coverage (`--cov-mode 0 -c 1.0 --min-seq-id 1.0`); the surviving non-redundant plasmids are the `ALL` table, each row tagged with the set of source DBs it was consolidated from. So this set contains **no duplicate plasmids across (or within) source databases**.

## Completeness gate (ledger)

| Completeness | n | disposition |
|---|---:|---|
| complete (closed) | 208,360 | **KEEP** |
| incomplete | 530,547 | discard (fragment/contig) |
| unknown ('-') | 113,693 | discard (completeness not asserted) |
| **total in ALL** | 852,600 | |

Working universe after the completeness gate: **208,360 complete plasmids** (348,617 raw source-entries collapsed into them, avg 1.67/plasmid — dedup consolidates, it does not drop biology).

## Source-database membership (complete set; a plasmid can belong to several)

| source | member plasmids | type |
|---|---:|---|
| IMG-PR | 136,318 | metagenome |
| RefSeq | 56,042 | isolate |
| GenBank | 55,428 | isolate |
| PLSDB | 47,215 | isolate |
| COMPASS | 12,084 | isolate |
| mMGE | 7,207 | metagenome |
| DDBJ | 4,331 | isolate |
| ENA | 4,049 | isolate |
| Kraken2 | 431 | isolate |
| TPA | 7 | isolate |

## Lifestyle composition (complete set)

| lifestyle | n | % |
|---|---:|---:|
| metagenomic | 136,887 | 65.7% |
| isolate | 64,835 | 31.1% |
| mixed | 6,638 | 3.2% |

## Source multiplicity (how many DBs each complete plasmid came from)

| distinct source DBs | plasmids |
|---:|---:|
| 1 | 147,146 |
| 2 | 20,083 |
| 3 | 31,224 |
| 4 | 7,467 |
| 5 | 2,381 |
| 6 | 58 |
| 7 | 1 |

61,214 complete plasmids are shared across >=2 source databases.

## Next in Phase A

- A.4 lab-made/synthetic flagging (conservative, logged).

- A.5 raw environment recovery per plasmid (IMG/PR GOLD ecosystem + NCBI BioSample); blanks kept as `unknown`.

