# PlasmidScope-primary — working set & environment (Phase A.4 + A.6)

**Pipeline:** `scripts/assemble_working_set.py` merging `complete_provenance.tsv`, `env_local.tsv` (IMG/PR GOLD, `scripts/recover_env_local.py`), the NCBI BioSample table `env_raw.tsv` (`scripts/fetch_biosample_env.py`), and the mMGE SRA table `sra_env_raw.tsv` (`scripts/fetch_sra_env.py`). Regenerate:
```
python3 scripts/assemble_working_set.py
```

## Inclusion funnel (full accounting)

| stage | n | note |
|---|---:|---|
| PlasmidScope `ALL` (deduplicated) | 852,600 | MMseqs2 100% id/cov |
| complete / closed | 208,360 | −644,240 incomplete/unknown (Phase A.2) |
| − lab-made / synthetic | −112 | flagged from organism/host (see below) |
| **= WORKING SET** | **208,248** | everything else kept, incl. `unknown` environment |

## Lab-made discards (conservative, logged)

**112** plasmids discarded as lab-made, each logged in `discarded_labmade.tsv` with the matching signal. Signals used: `synthetic construct`, `cloning vector`, `expression vector`, `shuttle vector`.

| signal | n |
|---|---:|
| cloning vector | 88 |
| expression vector | 12 |
| synthetic construct | 8 |
| shuttle vector | 4 |

## Environment coverage of the working set (raw, un-reconciled)

- Any environment label recovered: **178,972 / 208,248 (85.9%)**; the rest are `unknown` (kept, decided later).
| channel | plasmids with a label from it |
|---|---:|
| IMG/PR GOLD ecosystem | 111,708 |
| NCBI BioSample (INSDC) | 63,762 |
| mMGE SRA sample | 3,517 |

(A plasmid can draw from several channels — provenance and environments are kept as sets.)

## Working set by lifestyle

| lifestyle | n | % |
|---|---:|---:|
| metagenomic | 136,887 | 65.7% |
| isolate | 64,725 | 31.1% |
| mixed | 6,636 | 3.2% |

**Environment reconciliation (raw → harmonised habitat taxonomy) is deferred**, per the locked plan — labels are kept verbatim here for inspection first.

