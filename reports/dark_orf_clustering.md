# The dark plasmid proteome — clustering methodology and results

*Written 2026-08-28; numbers revised the same day after two pipeline defects were fixed (§Corrections).
Every number in this report is produced by a named script; nothing here is asserted without a
pipeline behind it. All scripts were re-run from the project root, and reproduce the figures below
exactly.*

> **Read this with `reports/pfam_dark_validation.md`.** That report searches these families against
> Pfam-A 38.2 and shows the "dark" set is **not one thing**: ~29% of the dark proteome is a PlasAnn
> annotation failure that Pfam corrects, and ~57% resists both. The R3 result below (85.9% have no
> PlasAnn-named homolog) is arithmetically unchanged but **must not be read as evidence of novelty**
> on its own.

## Question

`notebooks_my/cryptic_plasmids.ipynb` (cells 16–22) established that **71.5% of the coding capacity
of the payload-free small plasmidome is unannotatable**, and that **43.5% of those plasmids are dark
end to end** — every CDS labelled `ORF` by PlasAnn. That result is compatible with two opposite
readings, and the notebook cannot separate them:

- **Idiosyncratic junk** — nearly every dark ORF unique, consistent with non-adaptive persistence.
- **Recurrent dark families** — a bounded set of uncharacterized protein families appearing over and
  over, i.e. an annotation gap rather than an absence of function.

This report separates them by clustering the dark proteome itself.

## Provenance

| what | produced by |
|---|---|
| Plasmid ID list (71,414) | inline snippet in this run → `data/dark_orf_run/cryptic_small_ids.txt` |
| Protein extraction | `scripts/extract_plasann_proteins.py` → `data/dark_orf_run/faa_dark/`, `faa_all/` |
| Clustering | `scripts/cluster_dark_orfs.sh` → `data/dark_orf_run/dark30_*`, `mix30_*` |
| Result tables + summary | `scripts/summarize_dark_orf_clusters.py` → `darkfam_stats.tsv`, `darkfam_top200.tsv`, `dark_orf_summary.txt` |

Inputs, all pre-existing:

| input | produced by | documented in |
|---|---|---|
| `data/plasmidscope_primary/plasmid_metadata_master.tsv` (67 cols × 208,248) | `scripts/build_metadata_master.py` | `PROJECT_OVERVIEW.md` §5 |
| `data/plasann_run/gbk/{shard,mshard}_*.gbk.tar.gz` (1,073 tarballs, 6.1 GB) | `scripts/plasann_harvest_shard.py` | `PROJECT_OVERVIEW.md` §4 |
| `mob_cluster`, `hab_sub` columns | `scripts/harvest_typing.py`, `scripts/reconcile_environment.py` | `reports/typing_methodology.md`, `reports/environment_taxonomy_LOCKED.md` |

**Software.** MMseqs2 `18.8cc5c` and Python 3 (pandas 3.0.2, numpy 2.4.4) from the
`/gorilla/home/amujkic/.conda/envs/panaroo` and `.../genesis` environments respectively. 48 cores
available; all steps used 24 threads.

## Analysis set

The standing exclusion of `hab_top ∈ {Simulated-artifact, Lab-artifact}` gives **143,503** plasmids.
Small = `size_bp < 20,000` (**82,261**). Payload-free = zero on all three of `card_n_arg`,
`plasann_n_virulence`, `plasann_n_metal_biocide`, with NaN filled as 0 — **71,414** plasmids, the
same set analysed in notebook cells 16–22.

> **Known definitional leak, carried over from the notebook:** 562 of these 71,414 carry a PlasAnn
> `Antibiotic Resistance` CDS. `card_n_arg` counts only RGI Perfect+Strict hits, and all 863 RGI hits
> on those plasmids were `Loose` (median 70.8% reference coverage). ~171 sit at ≥95% identity and
> ≥90% coverage, i.e. genuine ARGs missed on CARD's curated bitscore cutoff. That is 0.24% of the
> set and does not affect any conclusion below, but it belongs in any methods section using this
> definition.

## Method

### 1. Recovering the dark proteome

PlasAnn's per-gene TSV shards (`data/plasann_run/annot/`) carry gene names and categories but **no
sequence**. The GenBank output (`data/plasann_run/gbk/`) carries `/translation=` for every CDS,
including those categorised `Open reading frame` — this is the only place in the pipeline where the
dark proteome is recoverable.

```bash
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
D=data/dark_orf_run

# dark ORF proteins only        -> headers  >{plasmid_id}|{cds_index}
python scripts/extract_plasann_proteins.py \
    --ids $D/cryptic_small_ids.txt --out $D/faa_dark --mode dark --jobs 24

# every CDS, tagged dark/named  -> headers  >{plasmid_id}|{cds_index}|{D|N}
python scripts/extract_plasann_proteins.py \
    --ids $D/cryptic_small_ids.txt --out $D/faa_all  --mode all  --jobs 24
```

Both complete in ~11 s over the 1,073 tarballs with 24 workers.

- `--mode dark` → **378,552** proteins (after deduplication, see the note below). This matches the notebook's ORF fraction
  (71.5% × 529,341 CDS ≈ 378.5k) to within 0.24%.
- `--mode all` → **456,949** proteins = 378,552 dark + **78,397** named.

> **Worth a methods sentence:** the notebook counted **150,789** *named annotation rows*, but only
> **78,397** of them are protein-coding. The difference is non-coding features — `oriV` (18,817),
> `oriT` (8,443), `RNAI` (6,856), the Col replicon markers and `Tn1681` — which have no
> `/translation`. **The named fraction of the proteome is far smaller than the named fraction of the
> annotation.**

### 2. Clustering

```bash
bash scripts/cluster_dark_orfs.sh     # ~42 s wall, 24 threads
```

which runs, after concatenating the per-shard FASTAs:

```bash
mmseqs easy-cluster $D/dark_orfs.faa $D/dark30 $D/tmp_dark \
    --min-seq-id 0.3 -c 0.8 --cov-mode 0 --threads 24 -v 1

mmseqs easy-cluster $D/all_cds.faa   $D/mix30  $D/tmp_mix  \
    --min-seq-id 0.3 -c 0.8 --cov-mode 0 --threads 24 -v 1
```

30% identity with 80% bidirectional coverage (`--cov-mode 0`) is a conventional protein-family
threshold. **Only this one threshold was run — no sensitivity sweep** (see Limitations).

The script deletes the two `*_all_seqs.fasta` files, which are redundant re-emissions of the input
(~190 MB).

### 3. Summarising

```bash
python scripts/summarize_dark_orf_clusters.py     # ~9 s
```

Lineage spread uses `mob_cluster` (MOB-suite primary cluster) as the lineage proxy, consistent with
`reports/small_cryptic_methodology.md` §3; habitat spread uses `hab_sub`. Within this subset
`mob_cluster` is called for **100%** of plasmids, resolving to **4,174** distinct lineages.

## Results

### R1 — The dark proteome is overwhelmingly recurrent

378,552 dark ORF proteins → **92,752 families**. These are real proteins, not spurious short calls:
median **138 aa**, q25–q75 86–249, only **3.9% under 50 aa**.

| family size | families | dark ORFs | % of dark ORFs |
|---|---:|---:|---:|
| singletons | 57,935 | 57,935 | **15.3%** |
| ≥2 members | 34,817 | 320,617 | **84.7%** |
| ≥5 members | 11,755 | 261,907 | 69.2% |
| ≥10 members | 5,623 | 222,506 | 58.8% |
| ≥50 members | 937 | 131,264 | 34.7% |
| ≥100 members | 394 | 94,220 | 24.9% |
| ≥1000 members | 6 | 9,334 | 2.5% |

**Only 15.3% of the dark plasmidome is idiosyncratic.** The idiosyncratic-junk reading is not
supported.

### R2 — The recurrence is not clonal redundancy

The obvious confound is one lineage sequenced repeatedly — the same trap as the Acinetobacter signal
in notebook cell 10, where a genus-level pattern collapsed to two MOB clusters. It does not apply
here:

| | families | % |
|---|---:|---:|
| families with ≥2 members | 34,817 | — |
| …spanning ≥2 MOB lineages | 15,262 | 43.8% |
| …spanning ≥5 MOB lineages | 2,921 | 8.4% |
| …confined to one lineage | 19,555 | 56.2% |

Because the cross-lineage families are the large ones, they dominate the sequence mass:
**241,727 dark ORFs (63.9% of all of them) sit in families spanning ≥2 independent MOB lineages.**

Restricting to the 5,623 substantial families (≥10 members):

- **79.2% (4,455) span ≥2 MOB lineages**
- median **4 lineages** and **5 habitats** per family

Most lineage-widespread families (full 200 in `darkfam_top200.tsv`):

| representative | proteins | plasmids | lineages | habitats | rep aa |
|---|---:|---:|---:|---:|---:|
| `GenBank_CP048556.1\|10` | 1,089 | 1,079 | **224** | 26 | 199 |
| `IMGPR_plasmid_3300038695_000007\|9` | 1,536 | 1,524 | 184 | 31 | 317 |
| `IMGPR_plasmid_3300022496_000003\|3` | 1,356 | 1,338 | 153 | 31 | 95 |
| `IMGPR_plasmid_3300047678_000089\|1` | 561 | 558 | 151 | 18 | 227 |
| `IMGPR_plasmid_3300007123_000008\|1` | 1,255 | 1,248 | 148 | 31 | 321 |

The top family spans **224 of the 4,174 lineages (5.4%)** and 26 habitats. Conservation of a protein
across hundreds of independent lineages and dozens of habitats is not the signature of neutral junk.

### R3 — The dark families have no named homolog in the dataset

Co-clustering all 456,949 CDS proteins (dark + named together) tests whether the dark ORFs are simply
proteins PlasAnn named elsewhere but missed here. Mostly, they are not:

| | families | dark ORFs absorbed | % of dark |
|---|---:|---:|---:|
| families containing **both** dark and named members | 1,674 | 53,219 | **14.1%** |
| **dark-only** families | 91,159 | 325,333 | **85.9%** |
| dark-only families with ≥10 members | 5,028 | 172,197 | 45.5% |

**85.9% of dark ORFs belong to families with no named homolog anywhere in the annotated set**, and
**5,028 dark-only families have ≥10 members**, holding almost half of all dark protein content.

## Interpretation

The dark plasmidome is not an unstructured bag of unique sequences. It is a **large, bounded,
recurrent protein space** — on the order of 5,000 substantial families — that is conserved across
independent plasmid lineages and habitats, and that current plasmid annotation cannot name.

This reframes the "functional dark matter" thesis. An argument from non-adaptive persistence is not
supported by R2: neutral junk does not stay conserved across 224 lineages. The defensible claim is
positive rather than absential — *the small plasmidome encodes a substantial protein space that
existing annotation databases cannot see* — and it comes with a concrete deliverable, the 5,028
dark-only families ranked by lineage spread in `darkfam_stats.tsv`.

## Limitations

1. **"Dark-only" is scoped to this dataset — and Pfam has since narrowed it.** R3 shows no homolog
   among *proteins PlasAnn named in this corpus*. It does **not** show absence from external
   databases. `reports/pfam_dark_validation.md` ran that test: **28.5% of the dark proteome is named
   outright by Pfam-A 38.2**, so a third of R3's 85.9% is annotation failure rather than novelty.
   What survives is the 57.4% that neither PlasAnn nor Pfam can name, and the 3,343 substantial
   Pfam-negative families in `data/pfam_run/dark_novel_families.tsv`.
2. **One clustering threshold.** 30% identity / 80% coverage only. Family counts are threshold-
   dependent; the qualitative split (recurrent vs. singleton) should be checked at 40% and 50%
   before publication.
3. **`mob_cluster` is a proxy, not a phylogeny.** It is mash-distance based, so "independent
   lineages" means "distinct primary clusters", not demonstrated phylogenetic independence. The
   subset is also lineage-skewed — cluster `AA379` alone holds 19,427 of the 71,414 plasmids.
4. **PlasAnn's ORF calls are inherited, not verified.** Gene calling and the `Open reading frame`
   category come from PlasAnn; no independent ORF calling (e.g. Prodigal) was run as a cross-check.
5. **Per-gene coverage is 92.2%.** 5,569 of the 71,414 plasmids returned no PlasAnn gene rows and
   contribute nothing here.

## Outputs

All under `data/dark_orf_run/` (464 MB total):

| file | size | what |
|---|---:|---|
| `cryptic_small_ids.txt` | 2.0 MB | the 71,414 plasmid IDs |
| `faa_dark/` | 82 MB | per-shard dark ORF proteins (1,073 `.faa`) |
| `faa_all/` | 102 MB | per-shard all-CDS proteins, dark/named tagged |
| `dark_orfs.faa` | 85 MB | concatenated dark proteome, 378,552 seqs |
| `all_cds.faa` | 107 MB | concatenated full proteome, 456,949 seqs |
| `dark30_cluster.tsv` | 24 MB | dark-only clustering, `rep → member` |
| `dark30_rep_seq.fasta` | 19 MB | 92,752 dark family representatives |
| `mix30_cluster.tsv` | 30 MB | dark+named co-clustering, `rep → member` |
| `mix30_rep_seq.fasta` | 21 MB | 97,500 co-clustering representatives |
| `darkfam_stats.tsv` | 5.9 MB | **per-family table**: proteins, plasmids, lineages, habitats, rep length |
| `darkfam_top200.tsv` | 14 KB | 200 most lineage-widespread families |
| `dark_orf_summary.txt` | 1.3 KB | the printed summary |

## Reproducing end to end

```bash
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
D=data/dark_orf_run
PY=/gorilla/home/amujkic/.conda/envs/genesis/bin/python

# cryptic_small_ids.txt must exist (regenerate from the master table if not)
$PY scripts/extract_plasann_proteins.py --ids $D/cryptic_small_ids.txt \
        --out $D/faa_dark --mode dark --jobs 24
$PY scripts/extract_plasann_proteins.py --ids $D/cryptic_small_ids.txt \
        --out $D/faa_all  --mode all  --jobs 24
bash scripts/cluster_dark_orfs.sh
$PY scripts/summarize_dark_orf_clusters.py
```

Total wall time ~1 minute on 24 threads.

## Corrections (2026-08-28)

Two defects were found while running the Pfam validation, and every number above is post-fix. Both
are documented in full in `reports/pfam_dark_validation.md` §6.

1. **Duplicate protein records (0.22%).** 144 plasmids appear in *both* PlasAnn shard series, so
   globbing both — which is mandatory — emitted their CDS twice. Dark ORFs 379,397 → **378,552**;
   recurrence 84.75% → **84.70%**. Fixed by a first-occurrence-wins dedup in
   `scripts/cluster_dark_orfs.sh`. No conclusion changed.
2. **CDS indexing broke the dark↔all-CDS join (10.7%).** `scripts/extract_plasann_proteins.py`
   numbered CDS within the emitted subset, so the same protein had different ids in `dark_orfs.faa`
   and `all_cds.faa`. It affected no number *in this report* (R1–R3 never join the two files), but
   it silently corrupted the first Pfam cross-tab. Fixed by indexing every translated CDS regardless
   of mode; the join is now asserted, not assumed.
