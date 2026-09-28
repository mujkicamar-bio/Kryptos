# Rebus

A Snakemake pipeline that takes a collection of complete bacterial plasmids and produces a
complete, evidence-rich annotation of every plasmid ORF and of its dark subset — the
proteins that no database names.

The reference collection is 143,503 complete plasmids (PlasmidScope and IMG/PR, with
simulated and lab-artifact ecosystems excluded), giving 9,317,050 called ORFs and 3,497,616
unique protein sequences.

---

## What this produces, and what it deliberately does not

The deliverable is **complete annotation**. One row per ORF and one row per dark family,
with every piece of evidence side by side.

The pipeline produces no ranking: no `top_1000`, no `candidate_score`, no `novelty_score`
and no `experimental_rank`. Choosing which unknown proteins are worth testing is a separate
downstream decision, made *on* these tables rather than baked into them. The separation is
deliberate — a scoring function embedded in an annotation pipeline is a scientific claim
disguised as an implementation detail, and it cannot be examined by anyone reading the
output.

## The problem the design addresses

Selecting "dark" proteins is easy and nearly worthless on its own, because **absence of
annotation is also what a gene-calling artefact produces**. A shadow ORF on the
reverse-complement strand of a real gene is unannotated by construction — no database
contains it, because it is not a protein — and it is *conserved*, because the real gene
beneath it is conserved. It therefore survives every absence-based filter and looks like an
ideal candidate all the way to the bench.

So the pipeline is built to report something positive rather than something absent. Each
dark ORF accumulates independent, named lines of evidence that it is a real protein:
purifying selection, breadth across independent plasmid lineages, family membership,
genomic context, and a recognisable fold.

Plasmids suit this unusually well. They are small, gene-dense and modular, and one signal
available here is available almost nowhere else: **a dark ORF inside an integron cassette
array is a real gene by construction.** It carries an *attC* recombination site, and it has
been excised, mobilised, re-integrated and then retained under selection — direct evidence
of both existence and function, obtained with no homology at all.

---

## Stages

Output lands in fifteen numbered directories under `outdir`, one per stage.

| rules | directory | what it does |
|---|---|---|
| analysis_set, clonal_registry | `01_analysis_set` | plasmids in scope (plasmids with a eukaryotic host excluded), as ids and sequence; clonal registry over MOB clusters |
| orf_call, orf_index | `02_orf_calling` | Pyrodigal gene calling, with circular-origin repair |
| dereplicate, plasmidscope_import, cascade_selection, check_hmmer_z | `03_dereplication` | exact-identity dereplication, asserted lossless |
| artefact_screen | `04_orf_qc` | AntiFam and low-complexity artefact screen — flags, never discards; an AntiFam-flagged protein skips every annotation tier and is reported as `NOT_SEARCHED`; a low-complexity-flagged one skips the DIAMOND tiers |
| preflight, tier_search, cascade_resolve | `05_annotation_cascade` | the annotation cascade, T1…T5, self-narrowing |
| annotate_plasmids, feature_files | `06_annotation_tables` | the annotated plasmidome, plus GFF3 and GenBank |
| orthology | `07_orthology` | eggNOG-mapper over the named fraction: COG and KEGG terms |
| protein_labels | `08_protein_labels` | every label from every source, as each tool gives it, in one table |
| label_databases | `08_protein_labels` | plasmid label databases: TADB, BacMet, oriTDB, CARD, mobileOG-db, dbAPIS, Anti-CRISPRdb and AMRFinderPlus |
| target_eligibility | `09_target_eligibility` | target eligibility: unnamed, searched, not artefact-flagged |
| dark_set, protein_clustering, protein_families, family_network, plasmid_lineage | `10_clustering` | dark set, then MMseqs2 deep-homology clustering into families |
| recurrence, extract_cds, family_evolution, consensus_recheck | `11_distribution_and_evolution` | CDS recovery, codon alignments, dN/dS, RNAcode, consensus re-check |
| defence_search, defence_gembase, defence_systems, conjugation_systems, integrons, phage_plasmids, is_elements, structure_search, context_features, dark_cooccurrence | `12_context_and_structure` | DefenseFinder, CONJScan, IntegronFinder, ISEScan, geNomad, directons, context terms, dark sequence co-occurrence (`dark_cooccurrence.tsv`: pairs of unique dark protein sequences sharing a plasmid in more lineages than chance predicts), Foldseek + ProstT5 |
| synteny | `13_synteny` | gene-order conservation counted over lineages, at the gene (close) and family (intermediate) level |
| rarity | `14_rarity` | family rarity labels and the saturation curve |
| annotation_report | `15_report` | the deliverable: complete annotation as CSV, per ORF and per dark family (the family table includes each family's co-occurring partners) |


### The cascade

Five tiers, each handed only what the previous one could not explain. A protein stops
being searched once `narrow_at` (0.7, a user decision of 2026-09-25) of its length is
covered. A protein below that reaches every tier, except that T5 does not search proteins
Pfam or Swiss-Prot named (`skip_if_named_by`). It is FUNCTIONAL once `min_coverage` (0.5)
of its length is explained by informative hits. Proteins flagged by AntiFam are not
searched by any tier; proteins flagged for low complexity are not searched by the DIAMOND
tiers.

| tier | method | database | role |
|---|---|---|---|
| T1 | hmmsearch | Pfam-A at `--cut_ga` | Pfam's own assertion of family membership; high precision |
| T2 | hmmsearch | Pfam-A at `-E 1e-5 --domE 1e-5` | divergent homologues below the curatorial bar, kept as a lower-authority claim |
| T3 | pharokka | PHROG, CARD, VFDB | the phage tier |
| T4 | diamond | Swiss-Prot | highest label quality available, and nearly free |
| T5 | diamond | NCBI ClusteredNR (nr clustered at 90% identity), skipping proteins Pfam or Swiss-Prot named | reaches environmental sequence nothing else does |

---

## Quick start

The commands of steps 1 and 2 were run against a fresh clone before being written here.

```bash
git clone git@github.com:mujkicamar-bio/Rebus.git
cd Rebus
```

### 1. What works immediately, with no data and no databases

The workflow engine and the Python libraries install from `pyproject.toml` alone:

```bash
python3 -m venv .venv && .venv/bin/pip install -e ".[dev]"

.venv/bin/snakemake -s workflow/Snakefile --configfile config/test/config.yaml --lint
.venv/bin/python -m pytest -q -m "not slow"
```

The linter reports the workflow is in good condition, and the test suite passes (the slow
tool integration tests are deselected). This verifies
the checkout is complete and internally consistent, which is as far as anyone can get
without the reference data.

A dry run at this point is also informative, and is *meant* to fail:

```bash
.venv/bin/snakemake -s workflow/Snakefile --configfile config/test/config.yaml -n -c 2
```

It builds the DAG, validates every schema, and then stops with `MissingInputException`
naming the input files you have not supplied. That is the pipeline telling you precisely
what step 2 has to provide.

### 2. The environment that can actually run the pipeline

The external executables — DIAMOND, HMMER, MMseqs2, Foldseek, MAFFT and the rest —
are not Python packages, so a working run needs the conda environment rather than the pip
install above:

```bash
conda env create -p envs/plasmidann -f workflow/envs/plasmidann.yaml
export PATH="$PWD/envs/plasmidann/bin:$PATH"
```

One environment for the whole workflow, deliberately: the pre-flight rule cannot check an
environment it is not running in, so per-rule environments would make pre-flight
meaningless.

**The environment that holds the tools also runs Snakemake**, which is why
`workflow/envs/plasmidann.yaml` declares Snakemake itself. Snakemake executes a `script:`
directive with the interpreter running Snakemake, not with the first `python` on `PATH`.
Running Snakemake from somewhere else splits the workflow in two — shell tools resolve one
way, Python scripts another — and the split stays silent until a script imports something
only one of them has.

### 3. What you must supply

Neither the sequence data nor the reference databases are in this repository; together
they are several hundred gigabytes.

| what | where it is configured | approximate size |
|---|---|---|
| plasmid FASTA and master table | `config/config.yaml` → `input` | depends on your collection |
| Pfam-A (with `.dat`) | `config/cascade.yaml` → T1, T2 | 2.1 GB |
| AntiFam | `config/cascade.yaml` → `artefact_screen` | 50 MB |
| pharokka database bundle | `config/cascade.yaml` → T3 | 1.5 GB |
| Swiss-Prot (DIAMOND) | `config/cascade.yaml` → T4 | 264 MB |
| NCBI ClusteredNR (DIAMOND), built by `tools/download_clustered_nr.sh` | `config/cascade.yaml` → T5 | 208 GB |
| eggNOG data | `config/targets.yaml` → `orthology` | 50 GB |
| Foldseek target DB and ProstT5 | `config/config.yaml` → `foldseek_db`, `prostt5_model` | 20 GB |
| MacSyFinder models | `config/config.yaml` → `references.macsyfinder_models` | 100 MB |
| plasmid label databases | `config/config.yaml` → `labels` | 0.9 GB |
| AMRFinderPlus and its database | `config/config.yaml` → `amrfinder` | 0.24 GB, plus 1.0 GB for `envs/amrfinder` |
| CONJScan 2.1.0 and its MacSyFinder 2.1.6 | `config/config.yaml` → `references.conjscan_models`; `config/targets.yaml` → `conjugation` | 17 MB, plus 0.7 GB for `envs/conjscan` |

The last three are installed by two scripts, run once from the repository root with any
Python 3 (both use the standard library only) and `micromamba`, `mamba` or `conda` on `PATH`:

```bash
python tools/download_label_dbs.py    # data/refs/labels/<db>/, with VERSION, SOURCE, MANIFEST.sha256
python tools/install_tool_envs.py     # envs/conjscan, envs/amrfinder, data/refs/conjscan, data/refs/amrfinder
```

pharokka (tier T3) runs from its own environment, `envs/pharokka`, with its database in
`data/refs/pharokka`:

```bash
conda env create -p envs/pharokka -f workflow/envs/pharokka.yaml
envs/pharokka/bin/install_databases.py -o data/refs/pharokka
```

Both are idempotent: a second run downloads and builds nothing. AMRFinderPlus and CONJScan
run from their own environments, named by path in config, because their pins conflict with
the main one - DefenseFinder needs MacSyFinder 2.1.4, CONJScan 2.1.0 needs 2.1.6. Pre-flight
checks both executables, every database directory and its `VERSION`, and that the
MacSyFinder can read CONJScan's model grammar. `labels.required` and `amrfinder.required`
are `true`, so a missing database stops the run in pre-flight rather than turning into an
empty result.

`hmmer_z` in `config/cascade.yaml` is the number of unique proteins of your analysis set,
whether or not the cascade searches them.
hmmsearch reports `E = Z × P(score | null)`, so a stale `-Z` rescales every E-value in the
run. Rule `check_hmmer_z` stops the run after dereplication if the count differs by more
than 2%, and names the value to set.

### 4. Running it

```bash
snakemake -s workflow/Snakefile -n -c 8      # plan
snakemake -s workflow/Snakefile -c 96        # run
sbatch workflow/run_pipeline.sbatch          # on SLURM
```

A core count is **required**: the cascade tiers take their thread count from
`workflow.cores`, so `-c`/`-j` is not optional.

`preflight` runs first and confirms every executable and every configured database is
present, in seconds to a few minutes (a MacSyFinder version check has taken up to 3 min),
before any compute is spent. The submission scripts take the repository root from the
directory they are submitted from, so submit them from the repository root.

### One submission runs everything

```bash
sbatch workflow/run_pipeline.sbatch
```

runs every stage in one job, the structure search (rule structure_search, ProstT5 and
Foldseek) included, on
CPU. ProstT5 is much faster on a GPU, so the split is still available on request:

```bash
STRUCTURE_ON_GPU=1 sbatch workflow/run_pipeline.sbatch   # everything up to structure_search
sbatch workflow/structure_gpu.sbatch                     # structure_search on the GPU
sbatch workflow/run_pipeline.sbatch                      # the rest, resuming
```

A run that asked for a GPU and did not get one fails loudly rather than silently taking ten
times longer.

### The SLURM account

`workflow/run_pipeline.sbatch` and `workflow/structure_gpu.sbatch` carry an `#SBATCH -A`
account and a partition that are specific to the cluster this was developed on. Change
both before submitting anywhere else.

A test-scale configuration is provided in `config/test/`, differing from production only in
what a smaller set forces: the tier list, `hmmer_z` and `outdir`.

---

## Measured performance

From a benchmark run (`config/bench`, `workflow/bench_pipeline.sbatch`) on a 100-plasmid
set (5,504 cascade queries) on one 96-core node, with T5 on full NCBI nr rather than
ClusteredNR and the rule set of that run.

| stage | wall | peak RSS |
|---|---|---|
| T5 nr | 2.91 h | 180.1 GB |
| orthology (eggNOG) | 39.7 min | 7.2 GB |
| T1 Pfam `--cut_ga` | 6.6 min | 0.35 GB |
| T2 Pfam sub-GA | 6.5 min | 0.32 GB |
| T3 pharokka | 2.4 min | 8.2 GB |
| T4 Swiss-Prot | 6.3 s | 2.4 GB |
| all other rules of that run | 8.6 min combined | ≤ 6.8 GB |

Two properties matter when sizing a larger run.

**Memory carries; wall time does not.** A search tier's peak memory is set by the database,
not the query count. Two independent measurements of the nr tier at query counts differing
by 25× agreed within 2%, at 177–180 GB. So `mem_mb` measured on a small set is a real
number for a large one.

**The nr tier is mostly fixed cost.** Across three measurements — 152, 3,808 and 3,910
queries — the pass costs about `2,400 s + 1.9 s per query`. At small query counts nearly
all of that is the single pass over the 357 GB nr database; at large ones the slope dominates.
This is why the cascade runs one job per tier and is not sharded: sharding into *N* pieces
would pay the fixed cost *N* times.

---

## Tests

```bash
python -m pytest -q -m "not slow"    # unit and script-harness tests
python -m pytest -q                  # adds the tool integration tests
```

`tests/conftest.py` supplies the `snakemake` global that Snakemake injects, so a workflow
script can run against a small fixture in milliseconds. `tests/test_workflow_wiring.py`
resolves the rule graph of the test configuration, which needs Snakemake but no data.

---

## Sources

- Rodríguez del Río Á. *et al.* Functional and evolutionary significance of unknown genes from uncultivated taxa. *Nature* **626**, 377–384 (2024)
- Pavlopoulos G.A. *et al.* Unraveling the functional dark matter through global metagenomics. *Nature* **622**, 594–602 (2023)
- Mistry J. *et al.* Pfam: the protein families database in 2021. *Nucleic Acids Res.* **49**, D412–D419 (2021)
- Eberhardt R.Y. *et al.* AntiFam: a tool to help identify spurious ORFs. *Database* **2012**, bas003
- Tesson F. *et al.* Systematic and quantitative view of the antiviral arsenal of prokaryotes. *Nat. Commun.* **13**, 2561 (2022)
- Eddy S.R. Accelerated profile HMM searches. *PLoS Comput. Biol.* **7**, e1002195 (2011)
- Buchfink B., Reuter K. & Drost H.-G. Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nat. Methods* **18**, 366–368 (2021)
- Islam H., Sharma A., Blair J. & Lopatkin A.J. PlasAnn: a curated plasmid-specific database and annotation pipeline for standardized gene and function analysis. *Nucleic Acids Res.* **54**, gkaf1507 (2026) - its identity and coverage tiers only
- Alcock B.P. *et al.* CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Res.* **51**, D690–D699 (2023)
- Feldgarden M. *et al.* AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Sci. Rep.* **11**, 12728 (2021)
- Guan J. *et al.* TADB 3.0: an updated database of bacterial toxin–antitoxin loci and associated mobile genetic elements. *Nucleic Acids Res.* **52**, D784–D790 (2024)
- Pal C. *et al.* BacMet: antibacterial biocide and metal resistance genes database. *Nucleic Acids Res.* **42**, D737–D743 (2014)
- Liu G. *et al.* oriTDB: a database of the origin-of-transfer regions of bacterial mobile genetic elements. *Nucleic Acids Res.* **53**, D163–D168 (2025)
- Brown C.L. *et al.* mobileOG-db: a manually curated database of protein families mediating the life cycle of bacterial mobile genetic elements. *Appl. Environ. Microbiol.* **88**, e00991-22 (2022)
- Yan Y., Zheng J., Zhang X. & Yin Y. dbAPIS: a database of anti-prokaryotic immune system genes. *Nucleic Acids Res.* **52**, D419–D425 (2024)
- Dong C. *et al.* Anti-CRISPRdb v2.2: an online repository of anti-CRISPR proteins including information on inhibitory mechanisms, activities and neighbors of curated anti-CRISPR proteins. *Database* **2022**, baac010 (2022)
- Cury J. *et al.* Identifying conjugative plasmids and integrative conjugative elements with CONJscan. *Methods Mol. Biol.* **2075**, 265–283 (2020)
- Coluzzi C., Garcillán-Barcia M.P., de la Cruz F. & Rocha E.P.C. Evolution of plasmid mobility: origin and fate of conjugative and nonconjugative plasmids. *Mol. Biol. Evol.* **39**, msac115 (2022)
- Benjamini Y. & Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* **57**, 289–300 (1995)
- Li Y. *et al.* PlasmidScope: a comprehensive plasmid database with rich annotations and online analytical tools. *Nucleic Acids Res.* **53**, D179 (2025) - the primary data source
- Schoch C.L. *et al.* NCBI Taxonomy: a comprehensive update on curation, resources and tools. *Database* **2020**, baaa062 (2020) - taxdump of 2026-04-05, eukaryotic hosts
- Robertson J. & Nash J.H.E. MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies. *Microb. Genom.* **4**, e000206 (2018)
- Hyatt D. *et al.* Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics* **11**, 119 (2010)
- Larralde M. Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. *J. Open Source Softw.* **7**, 4296 (2022) - pyrodigal 3.7.1
- Nayfach S. *et al.* CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nat. Biotechnol.* **39**, 578–585 (2021) - the terminal-repeat criterion
- Smillie C. *et al.* Mobility of plasmids. *Microbiol. Mol. Biol. Rev.* **74**, 434–452 (2010) - the size cut-off
- Frith M.C. A new repeat-masking method enables specific detection of homologous sequences. *Nucleic Acids Res.* **39**, e23 (2011) - tantan 51
- Steinegger M. & Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nat. Biotechnol.* **35**, 1026–1028 (2017) - MMseqs2 18.8cc5c
- Bouras G. *et al.* Pharokka: a fast scalable bacteriophage annotation tool. *Bioinformatics* **39**, btac776 (2023) - pharokka 1.10.1
- Cantalapiedra C.P. *et al.* eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Mol. Biol. Evol.* **38**, 5825–5829 (2021) - eggNOG-mapper 2.1.12
- Suzek B.E. *et al.* UniRef clusters: a comprehensive and scalable alternative for improving sequence similarity searches. *Bioinformatics* **31**, 926–932 (2015)
- Durairaj J. *et al.* Uncovering new families and folds in the natural protein universe. *Nature* **622**, 646–653 (2023) - the family network
- Ondov B.D. *et al.* Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol.* **17**, 132 (2016) - Mash 2.3
- Katoh K. & Standley D.M. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Mol. Biol. Evol.* **30**, 772–780 (2013) - MAFFT 7.526
- Yang Z. & Nielsen R. Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. *Mol. Biol. Evol.* **17**, 32–43 (2000) - yn00, PAML 4.10.7
- Washietl S. *et al.* RNAcode: robust discrimination of coding and noncoding regions in comparative sequence data. *RNA* **17**, 578–594 (2011) - RNAcode 0.3.1
- Néron B. *et al.* MacSyFinder v2: improved modelling and search engine to identify molecular systems in genomes. *Peer Community J.* **3**, e28 (2023) - MacSyFinder 2.1.4 and 2.1.6
- Néron B. *et al.* IntegronFinder 2.0: identification and analysis of integrons across bacteria, with a focus on antibiotic resistance in Klebsiella. *Microorganisms* **10**, 700 (2022) - IntegronFinder 2.0.6
- Xie Z. & Tang H. ISEScan: automated identification of insertion sequence elements in prokaryotic genomes. *Bioinformatics* **33**, 3340–3347 (2017) - ISEScan 1.7.3
- Camargo A.P. *et al.* Identification of mobile genetic elements with geNomad. *Nat. Biotechnol.* **42**, 1303–1312 (2024) - geNomad 1.12.0, database 1.9
- van Kempen M. *et al.* Fast and accurate protein structure search with Foldseek. *Nat. Biotechnol.* **42**, 243–246 (2024) - Foldseek 10.941cd33
- Heinzinger M. *et al.* Bilingual language model for protein sequence and structure. *NAR Genom. Bioinform.* **6**, lqae150 (2024) - ProstT5

## Licence

MIT. See [`LICENSE`](LICENSE). The MIT licence covers this code only, not the reference
data it downloads.

### Reference data whose terms restrict redistribution or use

None of these files is in the repository; each is downloaded by the installers above, and
its terms are recorded in `data/refs/labels/<db>/SOURCE` or the package metadata. Check
them before redistributing the databases, or outputs that reproduce their content, and
before any commercial use.

| data | terms |
|---|---|
| BacMet 2.0 | website footer "Copyright 2013-2018 All rights reserved"; no data licence stated, although the paper describes the database as freely available |
| CARD | free for non-commercial research or academic use by academic, government or non-profit institutions; commercial use needs a licence from McMaster University (card.mcmaster.ca/about, Terms of Use, sections 4 and 5) |
| CONJScan models | CC BY-NC-SA 4.0 (`data/refs/conjscan/CONJScan/metadata.yml`): non-commercial use, and derivatives under the same licence |
| dbAPIS | no data licence stated; website footer "Copyright 2023 YIN LAB, UNL. All rights reserved"; the article is CC BY 4.0 |
| Anti-CRISPRdb v2.2 | no data licence stated; website footer "Copyright CEFG 2021 All rights reserved"; the article is CC BY 4.0. The original hosts no longer serve the data, and the core dataset was taken from an Internet Archive capture |
| TADB 3.0, oriTDB 2.0 | no data licence stated; the papers describe the databases as freely available |
| mobileOG-db beatrix-1.6 | CC BY 4.0 (Zenodo archive) |
| AMRFinderPlus | public domain (NCBI) |
