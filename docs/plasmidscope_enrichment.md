# PlasmidScope annotations: what they cover, and how they can enter the pipeline

Produced 2026-09-23 by:

- `tools/plasmidscope_enrichment.py` (log: `data/PlasmidScope/annotation/analysis_set_enrichment.log`)
- `tools/plasmidscope_concordance.py` (log: `data/PlasmidScope/annotation/concordance.log`),
  3,000 random analysis-set plasmids, seed 20260916; cascade comparison against
  `results_test/` (100 test plasmids, 5,504 unique proteins, tiers T1-T4, no nr)

Source: PlasmidScope (Li et al., *Nucleic Acids Res.* 2025, 53:D179), dataset `ALL`,
downloaded from `https://plasmidapi.deepomics.org/api/database/files/ALL/data/ALL.<table>.tsv/`.

## 1. What was downloaded

`data/PlasmidScope/annotation/`, gzipped, covering all 852,600 `ALL` plasmids:

| Table | Content | Size |
|---|---|---|
| `ALL.protein_list` | one row per CDS: coordinates, product, eggNOG-mapper 2.1.12 (COG, GO, KEGG, CAZy, BiGG, Pfam, EC), protein sequence | 5.2 GB |
| `ALL.SP_list` | SignalP 6.0: type, class probabilities, cleavage site | 91 MB |
| `ALL.TMHs_list` | TMHMM 2.0 topology | 371 MB |
| `ALL.ARG_list` | RGI/CARD, including Loose hits | 4.6 GB |
| `ALL.VF_list` | VFDB (DIAMOND) | 102 MB |
| `ALL.CRISPR-Cas_list`, `ALL.SMs_list`, `ALL.trna_list` | CRISPRCasTyper, antiSMASH, ARAGORN | 4 MB |

Not downloaded: ESMFold structures (per plasmid only, no bulk file) and the FASTA/GenBank/GFF
archives (the protein table already carries the sequences).

Reduced to the analysis set: `analysis_set_orfs.tsv.gz` (one row per PlasmidScope ORF) and
`analysis_set_proteins.tsv.gz` (one row per unique protein, `seq_id` computed exactly as
`plasmidann.dereplicate` does, so it joins to every pipeline output).

## 2. Findings

**Coverage.** 143,063 of 143,503 analysis-set plasmids (99.7%); 8,753,749 ORFs;
3,312,766 unique proteins.

**PlasmidScope's ORFs are not one method.** For COMPASS, IMG/PR, PLSDB, mMGE and Kraken2
the ORFs are Prodigal 2.6 calls with Prokka products. For GenBank, RefSeq, DDBJ and EMBL
they are the submitters' deposited CDS (mostly NCBI PGAP). eggNOG-mapper was run on all.
`ps_class` is therefore built from the eggNOG fields only:

These counts use the first class rule (2026-09-23), under which any Pfam family or a COG
outside category S counted as ANNOTATED. The rule in force since 2026-09-24 is in section 3.

| ps_class (unique proteins) | n | % |
|---|---|---|
| ANNOTATED (Pfam, KEGG KO, EC, or COG outside category S) | 2,236,242 | 67.5 |
| UNKNOWN_ORTHOLOG (COG/OG in category S only) | 39,367 | 1.2 |
| NONE (no eggNOG identifier) | 1,037,157 | 31.3 |

**ORF agreement with our S1 caller** (share of our ORFs whose protein is identical to a
PlasmidScope protein on the same plasmid):

| PlasmidScope ORF source | identical | same stop, other start | absent |
|---|---|---|---|
| Prodigal (PlasmidScope's own calls) | 84.4% | 7.6% | 8.0% (median 130 aa) |
| deposited annotation | 74.2% | 11.3% | 14.5% (median 106 aa) |

About one ORF in five that our pipeline calls has no identical counterpart in PlasmidScope.
Those proteins need our own annotation whatever is decided below.

**PlasmidScope class against our cascade** (5,504 test proteins, before nr):

| | our FUNCTIONAL | DOMAIN_ONLY | UNCHARACTERIZED_HOMOLOG | NONE | total |
|---|---|---|---|---|---|
| PS ANNOTATED | 2,284 | 238 | 144 | 182 | 2,848 |
| PS UNKNOWN_ORTHOLOG + NONE | 154 | 46 | 159 | 1,333 | 1,692 |
| not in PlasmidScope | 449 | 44 | 43 | 428 | 964 |

- 9.1% of PlasmidScope-dark proteins are FUNCTIONAL by our cascade: its dark call is not
  final.
- 564 of the 2,617 proteins our cascade left NONE, UNCHARACTERIZED_HOMOLOG or DOMAIN_ONLY
  (21.5%) are ANNOTATED in PlasmidScope. Only NONE and UNCHARACTERIZED_HOMOLOG are
  target-eligible at the quality gate; of those, 326 of 2,289 (14.2%) are ANNOTATED in
  PlasmidScope. eggNOG gives no alignment span, so an eggNOG term says nothing about how
  much of the protein is explained.
- Of the 3,910 test proteins that reach the nr tier (explained fraction below
  `narrow_at` = 0.9 after T4), 39.8% are PS ANNOTATED, 42.0% PS dark and 18.2% not in
  PlasmidScope.

The test set is small (100 plasmids, 67 from IMG/PR), so these rates are indicative.

## 3. Decision (2026-09-23) and first run

Our ORF calls stay. A protein identical to a PlasmidScope protein that PlasmidScope
annotates (`ps_class` ANNOTATED) skips the cascade and is FUNCTIONAL at tier `PS`; every
other protein runs the full cascade, nr included. S4b orthology takes PlasmidScope's
eggNOG terms for those proteins and runs eggNOG-mapper only on the rest. Specification:
`PLASMID_ANALYSIS.md` section 14, Tier 0. Code: `src/plasmidann/plasmidscope.py`, rule
`plasmidscope_import` (`workflow/scripts/plasmidscope_import.py`), and changes to
`prepare_control.py`, `cascade_resolve.py` and `orthology.py`.

Test run on the 100 test plasmids (`results_test_ps/`, SLURM job 6957648, 1 h 45 min,
tiers T1-T4), compared with the run before the change (`results_test/`) by
`tools/compare_test_runs.py`:

| | before | after |
|---|---|---|
| unique proteins found in PlasmidScope | - | 4,691 of 5,304 (88.4%) |
| skipped the cascade (ANNOTATED) | - | 2,956 (55.7%) |
| cascade input (proteins + 100 controls + 100 decoys) | 5,504 | 2,548 |
| positive-control recall / decoy false-positive rate | 1.0 / 0.0 | 1.0 / 0.0 |
| ORFs in `plasmid_annotation.tsv`, all with a class | 5,673 | 5,673 |
| dark set (`dark_ids.txt`) | 2,142 | 1,787 |
| orthology rows from PlasmidScope / eggNOG-mapper / none | - | 2,956 / 418 / 199 |

Of the 355 proteins that left the dark set, 349 are PlasmidScope-annotated. The other six
are not caused by PlasmidScope directly: four now pass the relaxed Pfam tier (E-values
4.7e-6 to 9.3e-6) and two gain a pharokka hit, because `hmmer_z` counts only the sequences
searched and fell from 5,504 to 2,548, which halves every hmmsearch E-value. No protein
entered the dark set.

### Second run (2026-09-24): narrowed class rule, `-Z` over all proteins

Two changes, both made after the first run:

- `ps_class` ANNOTATED now needs a KEGG KO, an EC number, or a Pfam family that passes
  `cascade.is_informative`. A DUF/UPF family alone, or a COG with only a category letter,
  is UNKNOWN_ORTHOLOG and goes through the cascade.
- `hmmer_z` counts every unique protein plus the controls, searched or not (5,504 on the
  test set), so Tier 0 no longer changes any hmmsearch E-value.

Test run `results_test_ps/` (SLURM job 6963002, 1 h 39 min on 16 cores, tiers T1-T4),
compared with `results_test/` by `tools/compare_test_runs.py`:

| | before PlasmidScope | first rule | narrowed rule |
|---|---|---|---|
| skipped the cascade (ANNOTATED) | - | 2,956 | 2,724 (51.4%) |
| cascade input (proteins + 200 controls) | 5,504 | 2,548 | 2,780 |
| positive-control recall / decoy false-positive rate | 1.0 / 0.0 | 1.0 / 0.0 | 1.0 / 0.0 |
| target-eligible proteins | 2,154 | - | 1,959 |
| dark set (`dark_ids.txt`) | 2,142 | 1,787 | 1,947 |
| dark families (broad) | 2,041 | - | 1,855 |
| orthology rows from PlasmidScope / eggNOG-mapper / none | - | 2,956 / 418 / 199 | 2,724 / 489 / 198 |

Of the 195 proteins that left the dark set, 193 are PlasmidScope-annotated under the
narrowed rule. The other two now take a pharokka (T3) label (transcriptional regulator,
E = 5.3e-6; DNA binding protein, E = 6.3e-7). `-Z` does not reach pharokka: it searches the PHROG
profiles as queries against the proteins it is given (`mmseqs search` with the profile
database first, pharokka 1.10.1 `processes.run_mmseqs`), so its E-values scale with the
size of the T3 input, which Tier 0 halved; this is the one remaining way in which Tier 0 can move a dark call. The four
relaxed-Pfam losses of the first run are gone. No protein entered the dark set.

## 4. Other PlasmidScope tables, for later stages

- Stage 12 protein properties: `signal_peptide` from `ALL.SP_list` (SignalP 6.0) and
  `transmembrane_helix_count` from `ALL.TMHs_list` (TMHMM 2.0), joined on
  plasmid_id + Protein_ID, for identical proteins.
- Stage 8 context: VFDB, CRISPR-Cas, antiSMASH and tRNA features as neighbourhood labels.
- Not reused: `ALL.ARG_list` (our RGI run on the full set already exists and PlasmidScope's
  table includes Loose hits); ESMFold structures (no bulk download, and spec §47 uses
  ProstT5/3Di).
- PTUs: PlasmidScope has none. Its `Primary_Cluster_ID` / `Secondary_Cluster_ID` are
  MOB-suite 3.1.8 Mash clusters numbered per source database; plasmid grouping uses our
  own `mob_cluster`.
