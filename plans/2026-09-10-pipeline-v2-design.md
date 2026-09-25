Dark ORF Discovery Pipeline

Agent Implementation Specification

Purpose: Build a discovery-oriented annotation pipeline for plasmid-encoded ORFs from a large plasmid collection, with the primary goal of producing a complete, auditable evidence table from which approximately 1,000 putative "dark ORFs" can later be prioritized for experimental functional assays.

Core principle: The annotation pipeline must answer "What do we know about this ORF?". It must not decide which 1,000 ORFs are experimentally best. Experimental prioritization is a downstream workflow.

1. Project Goals

1.1 Primary goal

Given a collection of plasmid sequences and associated metadata:

predict ORFs;

reduce redundant computation while retaining all plasmid-level occurrences;

identify obvious spurious ORFs with AntiFam;

comprehensively annotate proteins using multiple evidence sources;

distinguish annotated proteins from dark (uncharacterized) proteins;

construct protein families;

characterize genomic context;

quantify evolutionary conservation and distribution;

characterize structural relationships;

record physical/sequence properties;

normalize database annotations and metadata;

produce complete, auditable output tables.

1.2 Non-goals

The first version must not:

produce a final biological ranking score;

remove proteins because they lack a functional hypothesis;

use genomic context as proof of function;

use structural similarity as proof of function;

treat every database hit as experimentally validated;

propagate a functional annotation blindly through a protein family;

equate "no structural hit" with "novel fold";

perform experimental candidate selection inside the annotation pipeline.

1.3 Discovery philosophy

The pipeline deliberately favors sensitivity and evidence preservation.

Every meaningful observation should be retained so that reporting thresholds and experimental-selection rules can be changed later without rerunning the entire analysis.

2. Terminology

2.1 ORF

A protein-coding ORF predicted by Prodigal/Pyrodigal.

2.2 ORF occurrence

A specific occurrence of an ORF on one plasmid at one set of genomic coordinates.

2.3 Unique protein sequence

A dereplicated protein sequence used to avoid repeated execution of expensive annotation tools.

One unique protein sequence can correspond to many ORF occurrences.

2.4 Annotation

An evidence-supported assignment of a protein to a known function, family, domain, orthologous group, or other interpretable biological category.

2.5 Dark / Uncharacterized

A protein is dark when its biological function is not sufficiently characterized by the complete annotation workflow. Uncharacterized, hypothetical, unknown-function, and DUF-like proteins are therefore included within the dark category, provided they pass the ORF/artifact criteria.

Darkness is an annotation state, not a statement that the protein is novel.

A dark protein may have:

many homologues;

a known structural fold;

a conserved genomic context;

strong evolutionary conservation;

repeated database deposition.

The pipeline may retain sublabels describing why a protein is dark, for example:

DARK_HYPOTHETICAL
DARK_DARK_DUF
DARK_NO_SEQUENCE_FUNCTION

These are subclasses of DARK, not separate annotation states.

2.7 Structural darkness

Two structural states are used:

DARK_NO_STRUCTURE
DARK_FOLD_KNOWN

DARK_FOLD_KNOWN means the protein remains functionally unresolved but has a convincing structural relationship to a known fold/structural family.

DARK_NO_STRUCTURE means no convincing known structural relationship was detected by the configured structural workflow.

Do not rename DARK_NO_STRUCTURE to "novel fold" automatically.

3. Input Data

3.1 Primary input

A collection of approximately 143,504 plasmid sequences retrieved from PlasmidScope.

3.2 Metadata

Metadata may be supplemented from:

PlasmidScope

PLSDB

NCBI

IMG/PR

MOB-suite and related plasmid resources

other explicitly configured sources

Typical metadata:

plasmid_id
sequence_source
accession
version
length
topology
host
host_taxonomy
habitat
MOB_class
MOB_cluster
plasmid_type
collection/source metadata
database provenance

3.3 Input data rules

Every source field must be retained.

Do not overwrite raw metadata.

Use:

raw_<field>
normalized_<field>

where normalization is required.

Example:

raw_host = "Escherichia coli K-12"
normalized_host_species = "Escherichia coli"

4. Data Model

The pipeline should use a normalized relational-style internal data model even if final deliverables are TSV/CSV.

4.1 Required logical tables

plasmids
orf_occurrences
unique_proteins
orf_annotations
orf_families
orf_context
orf_evolution
orf_structure
orf_properties
annotation_normalization
evidence_records

4.2 Required identifiers

Every entity must have a stable identifier.

Recommended:

plasmid_id
orf_occurrence_id
protein_id
family_id
annotation_id
context_id
structure_id

Do not use array indices as biological identifiers.

IDs must remain stable between pipeline reruns when inputs are unchanged.

5. Global Pipeline Rules

5.1 Raw data preservation

Never overwrite raw tool output.

Store:

raw/
intermediate/
normalized/
final/
logs/

Tool outputs should be copied or parsed into normalized tables, but the original outputs should remain available for auditing where practical.

5.2 Reproducibility

Every major output must contain or be associated with:

software version;

database version;

database date;

command/configuration;

pipeline version;

input dataset hash when appropriate.

5.3 Missing values

Use one explicit missing-value representation in final tables:

NA

Do not mix:

.
-
unknown
none
blank
NULL
NaN

NA means "not available/not observed/not applicable according to the field semantics."

The normalization layer must translate source-specific missing values to NA.

5.4 Evidence versus labels

A label summarizes an evidence state.

The underlying measurements must still be retained.

Example:

defence_associated = TRUE

must coexist with:

defence_system = CBASS
defence_system_source = DefenseFinder
distance_to_dark_orf = 2

6. Stage 1 — ORF Prediction

6.1 Tool

Use:

Prodigal / Pyrodigal

No second gene caller is required in version 1.

6.2 Input

Plasmid nucleotide FASTA.

6.3 Output

One record per predicted ORF occurrence.

Required fields:

plasmid_id
orf_occurrence_id
protein_id_pre_derep
start
end
strand
nucleotide_length
protein_length
protein_sequence
partial
start_type

Additional fields may be retained when available from Prodigal.

6.4 Important distinction

orf_occurrence_id is occurrence-specific.

protein_id should identify the dereplicated protein sequence after deduplication.

7. Stage 2 — ORF QC / Artifact Screen

7.1 Scope

Version 1 uses only:

Prodigal/Pyrodigal output properties;

AntiFam.

No additional artifact detector should be introduced without a deliberate pipeline revision.

7.2 AntiFam

Run AntiFam against unique protein sequences.

Record:

antifam_hit
antifam_model
antifam_score
antifam_evalue
antifam_description

7.3 Hard exclusion rules for dark-candidate discovery

The following are excluded from the dark-candidate discovery universe:

AntiFam-positive ORFs;

partial coding sequences;

proteins shorter than 20 amino acids.

These records must remain in the master annotation dataset with their exclusion reason.

7.4 Why <20 aa is excluded

The 20-aa cutoff is an operational discovery threshold, not a claim that proteins below 20 aa cannot be real.

Below approximately this size:

random ORFs become increasingly frequent;

sequence similarity is difficult to interpret;

profile/domain evidence becomes sparse;

structure prediction becomes unreliable or uninformative;

family clustering becomes more vulnerable to chance matches.

However, short microproteins can be biologically real and important.

Therefore:

protein_length < 20

must be labeled:

discovery_excluded_short = TRUE

rather than deleted from the master dataset.

Future versions may create a dedicated microprotein branch.

7.5 Other properties are labels only

Do not exclude automatically because of:

low_complexity
cysteine_rich
transmembrane
signal_peptide
overlap
unusual_amino_acid_composition

These can be biologically meaningful.

Record them later as properties/liabilities.

8. Stage 3 — Dereplication

8.1 Purpose

Avoid sending identical protein sequences through the expensive annotation workflow multiple times.

8.2 Critical rule

Dereplication changes computation, not biological occurrence counts.

For every unique protein:

protein_id
protein_sequence

retain a mapping to all source ORFs:

protein_id
orf_occurrence_id
plasmid_id

8.3 Required mapping

protein_id -> many orf_occurrence_id
orf_occurrence_id -> exactly one plasmid_id

Do not lose occurrence information.

9. Stage 4 — Annotation Normalization Layer

This stage is required before annotation adjudication.

9.1 Purpose

Database annotation strings are heterogeneous and often misleading.

Examples:

hypothetical protein
hypothetical protein
uncharacterized protein
putative protein
predicted protein
unknown protein
conserved hypothetical protein
protein of unknown function
DUF1234

must be normalized into controlled classes.

9.2 Required normalized annotation classes

Initial vocabulary:

FUNCTIONALLY_CHARACTERIZED
CURATED_UNKNOWN_FUNCTION
HYPOTHETICAL_PROTEIN
UNCHARACTERIZED_PROTEIN
PUTATIVE_PROTEIN
PREDICTED_PROTEIN
UNKNOWN_PROTEIN
DUF_OR_UNKNOWN_FAMILY
PSEUDOGENE
OTHER

9.3 Raw strings must remain

Store:

raw_annotation
normalized_annotation_class
normalization_rule
normalization_version

Never destroy the original annotation.

9.4 Functional characterization policy

An annotation is not considered functionally characterized merely because:

it appears in Swiss-Prot;

it has a name;

it has high sequence identity;

another database copied the description.

Evidence provenance must be inspected.

10. Stage 5 — Functional Annotation Hub

The goal is broad evidence collection.

Do not collapse the results immediately into one field.

Run and retain results from multiple complementary evidence sources.

10.1 T1 — Pfam-A high-confidence domain detection

Tool:

hmmscan

Database:

Pfam-A

Threshold:

--cut_ga

Purpose:

HIGH_CONFIDENCE_DOMAIN_ANNOTATION

Required fields:

pfam_id
pfam_acc
pfam_description
pfam_ga_score
pfam_evalue
pfam_domain_score
pfam_domain_evalue
pfam_alignment_start
pfam_alignment_end
pfam_query_coverage

Do not interpret a Pfam domain as complete functional characterization automatically.

10.2 T2 — Pfam relaxed search

Run a more permissive Pfam-A search.

Purpose:

CANDIDATE_REMOTE_OR_PARTIAL_HOMOLOGY

Required fields analogous to T1.

The relaxed result should be labeled:

pfam_relaxed_hit = TRUE

and should not automatically become the primary annotation.

10.3 T3 — InterPro

Use InterPro-compatible annotations/signatures where available.

Purpose:

INTEGRATED_DOMAIN / SIGNATURE EVIDENCE

Retain:

interpro_id
interpro_type
interpro_description
source_database
score/evalue where available
coverage

InterPro output should not be treated as independent biological evidence if it is derived from the same underlying family/domain model as another hit.

10.4 T4 — eggNOG-mapper

Run eggNOG-mapper on the ORF protein set.

Record:

eggnog_ortholog_group
eggnog_description
eggnog_best_og
eggnog_function
eggnog_evalue
eggnog_score
eggnog_propagated_annotation

Classification must preserve the distinction between:

orthology detected

and

experimentally supported function

10.5 T5 — DIAMOND Swiss-Prot

Purpose:

CURATED SEQUENCE HOMOLOGY

Required fields:

diamond_swissprot_target
diamond_swissprot_description
diamond_swissprot_bitscore
diamond_swissprot_evalue
diamond_swissprot_identity
diamond_swissprot_query_coverage
diamond_swissprot_subject_coverage
diamond_swissprot_alignment_length

A high-scoring hit to a hypothetical protein is not a validated functional annotation.

10.6 T6 — Remote HMM-HMM label

Use an HMM-HMM method such as HH-suite/HHsearch if configured.

Purpose:

REMOTE_HOMOLOGY_DETECTION

Important rule:

HMM-HMM is a label/evidence subtype, not a separate independent evidence category.

Do not count:

Pfam
HHsearch
DIAMOND
MMseqs2

as four independent biological observations when they represent the same underlying homology.

Required fields:

hmm_hmm_match
hmm_hmm_target
hmm_hmm_score
hmm_hmm_probability
hmm_hmm_evalue
hmm_hmm_coverage

10.7 T7 — DIAMOND broader database

Run a broad protein database search such as NCBI nr or another configured comprehensive protein database.

Purpose:

BROAD_HOMOLOGY_DISCOVERY

Do not use broad-database labels as automatic functional truth.

Required fields:

diamond_nr_target
diamond_nr_description
diamond_nr_bitscore
diamond_nr_evalue
diamond_nr_identity
diamond_nr_query_coverage
diamond_nr_subject_coverage

11. Search and Reporting Thresholds

The pipeline must separate:

Search threshold

Determines whether a tool is computationally allowed to report a candidate hit.

Evidence threshold

Determines whether a hit is biologically credible.

Reporting threshold

Determines whether a hit can be used as the primary annotation.

These must not be represented as one universal score.

11.1 No universal "90%" rule

Do not use a universal 90% explained threshold across:

Pfam;

HMM-HMM;

DIAMOND;

eggNOG;

InterPro.

Scores are not directly comparable.

11.2 Recommended implementation

Run broad/permissive searches where computationally practical.

Store all relevant raw results.

Allow adjudication after all evidence is collected.

This permits thresholds to be modified without rerunning the full pipeline when possible.

12. Annotation Adjudication

This stage converts raw hits into interpretable status labels.

12.1 Evidence hierarchy

A starting hierarchy:

1. experimentally characterized functional evidence
2. strong curated family/domain evidence
3. high-confidence orthology with informative function
4. strong curated sequence homology
5. remote sequence/profile homology
6. broad homology to known proteins
7. contextual/structural hypothesis
8. unresolved

This is an evidence hierarchy, not a scalar score.

12.2 Primary annotation

Every protein may have:

primary_annotation
primary_annotation_source
primary_annotation_confidence

while preserving all secondary observations.

12.3 Example

If:

Pfam = weak domain
Swiss-Prot = strong experimentally characterized homolog

then:

both results retained
Swiss-Prot = primary annotation

12.4 Transitive annotation safeguard

Do not automatically propagate:

A = known function
B = similar to A
C = similar to B

into:

C = known function

unless the evidence meets the configured annotation-transfer criteria.

13. Stage 6 — Dark/Uncharacterized Classification

This stage defines the central discovery population.

13.1 Top-level states

ANNOTATED
DARK

13.2 ANNOTATED

Use when the complete evidence set supports a sufficiently specific functional/family assignment according to the configured adjudication rules.

13.3 DARK

Use when no sufficiently supported functional assignment is available after the complete annotation workflow. This includes uncharacterized, hypothetical, unknown-function, conserved hypothetical, DUF-like and similarly unresolved proteins.

Dark does not imply novelty.

13.5 Recommended additional field

dark_evidence_level

Suggested semantic scale:

0 = not dark
1 = weakly supported unresolved sequence
2 = repeated sequence/homology but unresolved
3 = independently recurrent hypothetical/uncharacterized protein
4 = conserved dark family
5 = strongly conserved/recurrent dark protein with independent biological context

The exact numerical boundaries should be treated as configuration and revised after inspecting the dataset.

Do not use this level as the final experimental score.

14. Hypothetical / Unknown Annotation Evidence

A critical rule:

A protein deposited repeatedly as:

hypothetical protein
uncharacterized protein
unknown protein

is not equivalent to an isolated predicted ORF with no external observations.

The pipeline must estimate annotation evidence strength from:

number of independent database records;

number of independent source genomes/plasmids;

sequence recurrence;

sequence family size;

taxonomic breadth;

presence in independent database sources;

presence across independent research/database submissions.

Required fields may include:

hypothetical_deposit_count
independent_source_count
independent_database_count
independent_plasmid_count
hypothetical_recurrence_level

These increase confidence that the protein is real while leaving its function unresolved.

15. Stage 7 — MMseqs2 Protein Family Construction

Run MMseqs2 after the full annotation phase or at the configured point where all protein sequences are available.

15.1 Purpose

Identify protein families among all ORFs.

15.2 Family philosophy

Families are descriptive sequence relationships.

A family is not automatically:

a functional family;

an evolutionary clade;

an experimentally validated unit.

15.3 Required family columns

For each protein:

mmseq_family

For each family:

family_size
dark_member_count
percentage_dark_in_family
annotated_member_count
plasmid_count
host_count
genus_count
MOB_count
lineage_count
habitat_count
max_annotation_strength

15.4 Dark families

Families with:

percentage_dark_in_family = 100

are especially interesting.

But:

family darkness is a label, not an exclusion or ranking criterion.

15.5 Nested family option

The implementation should allow multiple clustering resolutions.

Recommended conceptual fields:

mmseq_family_close
mmseq_family_intermediate
mmseq_family_broad

Exact thresholds should be configurable rather than hard-coded into downstream logic.

16. Stage 8 — Genomic Context

Genomic context is independent from direct sequence annotation.

16.1 Input

For every dark ORF:

plasmid;

ORF coordinates;

neighboring ORFs.

16.2 Neighborhood extraction

Internally extract at least:

±10 ORFs

when available.

Generate a standard summary for:

±3 ORFs

The larger window is retained so that future analyses do not require re-parsing the plasmids.

16.3 Required neighbor fields

neighbor_orf
neighbor_distance_genes
neighbor_distance_bp
neighbor_strand
neighbor_relative_orientation
neighbor_annotation
neighbor_function_class

16.4 Context categories

Examples:

relaxase
mobilization
replication
partition
toxin_antitoxin
defense
restriction
integron
conjugation
plasmid_backbone
transcriptional_regulation
membrane
metabolic
hypothetical
other

16.5 DefenseFinder

Run DefenseFinder/system detection where configured.

Record:

defence_system
defence_system_type
defence_system_coordinates
defence_distance_to_dark_orf
defence_source

16.6 IntegronFinder

Run IntegronFinder where configured.

Record:

integron_id
integron_type
integron_coordinates
attC_present
distance_to_dark_orf

16.7 Context labels

Create labels:

defence_associated
integron_associated
mobilization_associated
toxin_antitoxin_associated
replication_associated
partition_associated
backbone_associated

These are hypothesis labels.

They are not definitive functional annotations.

17. Genomic Context Conservation / Synteny

The same dark protein may occur in multiple plasmids.

Compare neighborhoods across occurrences and homologous family members.

17.1 Example

A - B - DARK - C - D
A - B - DARK - C - D
A - B - DARK - C
A - B - DARK - C - D

This is evidence of conserved genomic context.

17.2 Required measurements

left_neighbor_conservation
right_neighbor_conservation
neighborhood_conservation
operon_like_conservation
context_recurrence
synteny_conservation

17.3 Important rule

Context recurrence is evidence that a locus may be biologically maintained.

It is not proof of the function of the dark protein.

18. Stage 9 — Evolutionary Analysis

Evolutionary analyses are primarily performed at the protein-family level.

18.1 Distribution variables

For every dark protein/family:

plasmid_count
host_count
species_count
genus_count
MOB_count
lineage_count
habitat_count

18.2 Sequence conservation

Measure where appropriate:

mean_pairwise_identity
identity_distribution
alignment_coverage
conserved_positions
conserved_motifs
amino_acid_conservation
nucleotide_conservation

18.3 RNAcode

Where sufficient homologous nucleotide sequence is available:

obtain codon-aware alignments;

run RNAcode;

retain result as comparative coding evidence.

Fields:

rnacode_support
rnacode_score
rnacode_pvalue
rnacode_coordinates

RNAcode supports evidence of conserved coding potential, not functional identity.

18.4 dN/dS

Use dN/dS only when:

the family has enough independent homologous sequences;

alignment quality is adequate;

codon alignment is defensible;

the evolutionary comparison is appropriate.

Do not force dN/dS on every dark ORF.

Fields:

dnds_available
dnds_value
dn_value
ds_value
selection_model
selection_significance

dN/dS is evidence, not a universal ranking variable.

19. Stage 10 — Protein Structural Analysis

All dark proteins should be considered for structural analysis, not only family representatives.

19.1 Structure prediction / representation

Use the configured structure workflow, including ProstT5 or equivalent sequence-to-3Di representation.

Record:

structure_model_id
structure_model_source
structure_confidence
model_length

19.2 Foldseek

Search predicted structures against an appropriate structural database.

Record:

foldseek_hit
foldseek_target
foldseek_description
foldseek_score
foldseek_tm_score_if_available
foldseek_alignment_coverage
foldseek_evalue_if_available

19.3 Structural states

At minimum:

NO_STRUCTURAL_MODEL
STRUCTURAL_MODEL_AVAILABLE
DARK_FOLD_KNOWN
DARK_NO_STRUCTURE

19.4 Structural interpretation

If a convincing structural match exists:

DARK_FOLD_KNOWN

If no convincing structural match exists:

DARK_NO_STRUCTURE

Do not state:

novel fold

unless an explicit downstream structural-novelty analysis supports that claim.

19.5 Structural annotation

A structural match generates a:

structural_hypothesis

not necessarily a:

functional_annotation

20. Stage 11 — Sequence / Physical Properties

Record descriptors for all dark proteins.

Required/desired fields:

protein_length
molecular_weight
isoelectric_point
hydrophobicity_summary
transmembrane_helix_count
signal_peptide
low_complexity_fraction
disorder_fraction
cysteine_fraction
coiled_coil_prediction
secondary_structure_summary

20.1 Liabilities

Generate labels such as:

SHORT
MEMBRANE
SIGNAL_PEPTIDE
DISORDERED
LOW_COMPLEXITY
CYS_RICH
UNCERTAIN_ORF
OTHER

These are not exclusion criteria unless explicitly requested by a future analysis.

21. Stage 12 — Distribution and Background Normalization

This is required for later comparative analyses.

21.1 Problem

Raw counts are confounded by:

plasmid length;

gene count;

taxonomic composition;

source database;

sampling bias;

habitat;

sequencing/project origin.

21.2 Required denominator/background fields

At plasmid level:

plasmid_gene_count
host_taxonomy
MOB_cluster
habitat
source_database

21.3 Raw and normalized quantities

For relevant features, retain both:

raw_count
raw_prevalence
normalized_prevalence

21.4 Examples of future comparisons

dark ORF enrichment in MOB group
dark family enrichment in host taxon
dark family enrichment in habitat
context enrichment by plasmid type

21.5 Large-plasmid confounding

A larger plasmid contains more genes and therefore has a higher chance of containing any feature.

Do not call a feature "enriched" without accounting for the relevant denominator.

The normalization layer should support methods such as:

per-gene normalization;

per-plasmid normalization;

size-matched backgrounds;

stratified comparison;

appropriate statistical models.

Exact statistical tests belong to a later analysis module, not the core annotation engine.

22. Stage 13 — Rarity versus Conservation

Record distribution rather than assigning "good" or "bad."

Every protein/family should be describable as:

rare
lineage_specific
plasmid_family_specific
widely_conserved
cross_MOB
cross_host
cross_taxon

These labels must be derived from retained measurements.

22.1 Important principle

A protein occurring in many plasmids is not automatically more interesting than one occurring in seven plasmids.

Examples:

highly conserved dark proteins may indicate ancient/maintained plasmid biology;

rare recurrent proteins may represent specialized functions;

lineage-specific proteins may represent niche adaptation;

cross-MOB dark proteins may indicate broadly mobile biology.

23. Stage 14 — Evidence Integration

The pipeline must not collapse all evidence into one numerical score.

23.1 Evidence categories

Use named evidence categories:

sequence_homology
domain_homology
orthology
genomic_context
evolutionary_conservation
distribution
structural_evidence
protein_properties
orf_qc

23.2 Nested evidence rule

Do not count redundant observations multiple times.

For example:

Pfam
InterPro/Pfam-derived domain
HHsearch
DIAMOND
MMseqs

may all derive from one underlying homology relationship.

They should remain individually visible, but the evidence framework must identify them as potentially non-independent.

23.3 Evidence summary

Create:

evidence_count_total
evidence_categories_present
independent_evidence_categories

Where:

independent_evidence_categories

counts categories rather than individual database hits.

24. Stage 15 — Functional Hypothesis Layer

The pipeline may generate broad hypotheses when evidence supports them.

Examples:

defence_associated
mobilization_associated
toxin_antitoxin_associated
replication_associated
partition_associated
regulatory_associated
membrane_associated
metabolic_associated
structural_associated
unknown

Each hypothesis must have traceable evidence.

Example:

functional_hypothesis = defence_associated

support:
- repeated adjacency to DefenseFinder system
- median distance = 2 genes
- observed in 17 independent plasmids
- conserved neighborhood in MMseqs family

Do not write:

annotation = defence protein

unless direct sequence/domain/experimental evidence justifies that statement.

25. Stage 16 — Final ORF Annotation Table

Primary deliverable:

final/annotation_complete.tsv

Recommended one-row-per-ORF-occurrence structure.

25.1 Core identity

plasmid_id
orf_occurrence_id
protein_id
family_id

25.2 Coordinates

start
end
strand
protein_length
partial

25.3 QC

antifam_hit
antifam_model
discovery_excluded_short
discovery_excluded_partial

25.4 Annotation state

annotation_status
primary_annotation
primary_annotation_source
primary_annotation_confidence
dark_status
dark_evidence_level

25.5 Pfam

pfam_id
pfam_ga_score
pfam_evalue
pfam_domain_score
pfam_domain_evalue
pfam_query_coverage
pfam_relaxed_hit

25.6 HMM-HMM

hmm_hmm_match
hmm_hmm_target
hmm_hmm_score
hmm_hmm_probability
hmm_hmm_coverage

25.7 DIAMOND Swiss-Prot

diamond_swissprot_target
diamond_swissprot_bitscore
diamond_swissprot_evalue
diamond_swissprot_identity
diamond_swissprot_query_coverage
diamond_swissprot_subject_coverage

25.8 DIAMOND broad DB

diamond_nr_target
diamond_nr_bitscore
diamond_nr_evalue
diamond_nr_identity
diamond_nr_query_coverage
diamond_nr_subject_coverage

25.9 eggNOG / InterPro

eggnog_ortholog_group
eggnog_function
eggnog_evalue
interpro_id
interpro_description

25.10 Family

mmseq_family
family_size
dark_member_count
percentage_dark_in_family
family_plasmid_count
family_host_count
family_genus_count
family_MOB_count
family_lineage_count

25.11 Context

defence_associated
integron_associated
mobilization_associated
toxin_antitoxin_associated
replication_associated
partition_associated
backbone_associated
neighborhood_conservation
synteny_conservation

25.12 Evolution

plasmid_count
host_count
species_count
genus_count
MOB_count
lineage_count
mean_pairwise_identity
alignment_coverage
conserved_motif_count
rnacode_support
dnds_available
dnds_value

25.13 Structure

structural_status
structure_model_id
structure_confidence
foldseek_hit
foldseek_target
foldseek_score
foldseek_alignment_coverage

25.14 Properties/liabilities

transmembrane_helix_count
signal_peptide
low_complexity_fraction
disorder_fraction
cysteine_fraction
liabilities

25.15 Hypothesis

functional_hypothesis
hypothesis_support_summary

26. Final Dark Family Table

Deliverable:

final/dark_families_complete.tsv

One row per dark family.

Required:

family_id
family_size
dark_member_count
percentage_dark_in_family
annotated_member_count
plasmid_count
host_count
species_count
genus_count
MOB_count
lineage_count
habitat_count
max_annotation_strength
family_structural_status
dominant_context
context_conservation
dark_evidence_level_summary

Also retain a normalized member table:

final/dark_family_members.tsv

with:

family_id
protein_id
orf_occurrence_id
plasmid_id
dark_status
annotation_status

27. Plasmid-Level Summary

Deliverable:

final/plasmid_annotation_summary.tsv

Recommended fields:

plasmid_id
plasmid_length
gene_count
annotated_gene_count
uncharacterized_gene_count
dark_gene_count
dark_family_count
defence_system_count
integron_count
mobilization_features
MOB_cluster
host
habitat
source_database

This provides the bridge from ORF-level discovery to plasmid-level biology.

28. Recommended Repository Structure

dark_orf_pipeline/
├── README.md
├── PIPELINE.md
├── pyproject.toml
├── config/
│   ├── config.yaml
│   ├── databases.yaml
│   └── normalization.yaml
├── workflow/
│   └── Snakefile
├── src/
│   └── dark_orf/
│       ├── io/
│       ├── orf/
│       ├── antifam/
│       ├── annotation/
│       ├── normalization/
│       ├── families/
│       ├── context/
│       ├── evolution/
│       ├── structure/
│       ├── properties/
│       ├── evidence/
│       └── utils/
├── agents/
│   ├── AGENT_GUIDE.md
│   ├── TASKS.md
│   └── HANDOFF.md
├── scripts/
├── tests/
├── data/
│   ├── raw/
│   ├── intermediate/
│   └── normalized/
├── results/
│   ├── annotation/
│   ├── families/
│   ├── context/
│   ├── evolution/
│   ├── structure/
│   └── final/
├── logs/
└── docs/

29. Agent Architecture

Coding should be divided into agents/modules with clear contracts.

Agent 1 — Input / Metadata Agent

Responsibilities:

validate plasmid FASTA;

ingest metadata;

normalize identifiers;

generate plasmids.tsv;

preserve source provenance;

implement missing-value normalization.

Deliverables:

plasmids.tsv
metadata validation report

Must not modify ORF annotations.

Agent 2 — ORF Agent

Responsibilities:

run Prodigal/Pyrodigal;

parse output;

generate orf_occurrences.tsv;

translate proteins;

assign stable identifiers.

Must retain coordinates and plasmid mapping.

Agent 3 — ORF QC Agent

Responsibilities:

run AntiFam;

mark partial ORFs;

apply <20 aa discovery exclusion;

generate QC flags.

Must not delete records from the master dataset.

Agent 4 — Annotation Agent

Responsibilities:

Pfam GA;

Pfam relaxed;

InterPro;

eggNOG;

Swiss-Prot;

broad database;

optional HMM-HMM.

Must save raw and normalized results independently.

Agent 5 — Annotation Normalization / Adjudication Agent

Responsibilities:

normalize source labels;

distinguish hypothetical vs functionally characterized;

assign annotation status;

assign primary annotation;

assign dark status;

assign dark evidence level;

prevent unsafe transitive annotation.

Agent 6 — Family Agent

Responsibilities:

run MMseqs2;

assign family IDs;

calculate family composition;

calculate dark-member percentage;

calculate family distribution.

Agent 7 — Context Agent

Responsibilities:

extract neighborhoods;

map neighbor annotations;

run DefenseFinder;

run IntegronFinder;

calculate context recurrence;

calculate synteny conservation;

assign contextual hypothesis labels.

Agent 8 — Evolution Agent

Responsibilities:

construct family alignments;

calculate conservation;

calculate lineage breadth;

run RNAcode where appropriate;

run dN/dS where appropriate;

write evolutionary evidence tables.

Agent 9 — Structure Agent

Responsibilities:

generate structure representations/models;

run Foldseek;

record structure status;

distinguish known-fold from no-known-match.

Agent 10 — Properties Agent

Responsibilities:

calculate sequence properties;

transmembrane predictions;

signal peptide;

disorder;

low complexity;

cysteine-richness;

experimental liability labels.

Agent 11 — Evidence Integration Agent

Responsibilities:

merge evidence tables;

calculate named evidence categories;

prevent duplicate counting;

generate hypothesis support summaries;

create annotation_complete.tsv.

Must not rank the top 1,000.

Agent 12 — Validation Agent

Responsibilities:

validate row counts;

validate identifier mappings;

check missingness;

check impossible values;

check annotation consistency;

generate QC report.

30. Agent Rules

Every coding agent must obey:

Rule 1

Do not invent biological thresholds that are not defined in this document or configuration.

Rule 2

Prefer configuration values over hard-coded thresholds.

Rule 3

Do not silently remove records.

Any filtering must produce:

filter_reason
filter_stage

Rule 4

Never overwrite raw tool output.

Rule 5

Never replace a raw database annotation with a normalized label.

Rule 6

Never convert a contextual hypothesis into a direct functional annotation without explicit evidence.

Rule 7

Never collapse evidence sources prematurely.

Rule 8

Never create a global candidate score in the annotation pipeline.

Rule 9

Every output table must have a schema and validation test.

Rule 10

Every tool execution must record:

software
version
database
database_version/date
parameters
input
output

31. Configuration

A central YAML configuration should contain all thresholds.

Example:

orf:
  min_dark_length_aa: 30
  exclude_partial: true

antifam:
  enabled: true

pfam:
  enabled: true
  ga: true
  relaxed: true
  relaxed_evalue: configurable

annotation:
  retain_all_hits: true

mmseqs:
  enabled: true
  family_thresholds:
    close: configurable
    intermediate: configurable
    broad: configurable

context:
  neighborhood_genes: 10
  summary_neighborhood_genes: 3

structure:
  enabled: true
  foldseek: true

evolution:
  rnacode: true
  dnds: conditional

output:
  missing_value: "NA"

Thresholds should be documented separately from biological interpretations.

32. Validation Requirements

Before producing final outputs, the pipeline must run validation checks.

32.1 Identity checks

Verify:

every ORF occurrence maps to one plasmid
every unique protein maps to >=1 ORF occurrence
every family member maps to one family

32.2 Sequence checks

Verify:

protein sequences contain valid amino-acid symbols;

translated ORFs have expected lengths;

no duplicate IDs with different sequences;

no orphaned proteins.

32.3 Annotation checks

Verify:

all annotation sources retain raw identifiers;

primary_annotation exists only when evidence criteria are met;

dark proteins have no qualifying primary functional annotation;

hypothetical labels are not interpreted as functional assignments.

32.4 QC checks

Verify:

AntiFam-positive ORFs are flagged;

partial ORFs are flagged;

<20 aa proteins are flagged;

excluded records remain recoverable.

32.5 Family checks

Verify:

0 <= percentage_dark_in_family <= 100
dark_member_count <= family_size
annotated_member_count <= family_size

32.6 Distribution checks

Verify:

host_count <= plasmid_count
genus_count <= host_count where definitions require it
MOB_count <= plasmid_count

Exact inequalities depend on the unit definitions and must be encoded explicitly.

33. Provenance Requirements

Each final annotation should be traceable to:

plasmid sequence
ORF coordinates
protein sequence
annotation tool
database
database version
database date
search parameters
raw hit
normalization rule
adjudication rule

The user must be able to take one row from:

annotation_complete.tsv

and trace it back to its source plasmid and original tool output.

34. Interpretation Rules

These rules are mandatory.

34.1 "No Pfam"

Means:

no qualifying Pfam result

not:

novel protein

34.2 "No Swiss-Prot hit"

Means:

no sufficiently significant Swiss-Prot similarity detected

not:

unknown biology

34.3 "Hypothetical protein"

Means:

sequence exists in a database without an established functional description

and can constitute high-confidence evidence of protein existence when independently recurrent.

34.4 "Known fold"

Means:

structural relationship detected

not:

function proven

34.5 "No structural hit"

Means:

no convincing structural relationship detected under the configured workflow

not necessarily:

novel fold

34.6 "Dark family"

Means:

family members are predominantly or completely unresolved

not:

family is biologically important

35. Annotation Output Versus Experimental Prioritization

This boundary must remain explicit.

Annotation pipeline output

Answers:

What do we know about this ORF?

Produces:

complete evidence
labels
scores
annotations
families
context
evolution
structure
properties
provenance

Downstream prioritization pipeline

Answers:

Which unresolved ORFs are most worth testing?

This later workflow may use:

conservation;

cross-MOB distribution;

family size;

family darkness;

genomic context;

defence association;

structural novelty;

expression/evidence if added;

experimental feasibility;

technical assay suitability.

But those choices belong downstream.

36. Future Candidate Prioritization — Not Part of Version 1

The eventual top-1,000 workflow should consume:

final/annotation_complete.tsv
final/dark_families_complete.tsv

and generate a candidate portfolio.

It should not modify the underlying annotation records.

Possible future candidate classes:

highly conserved dark proteins
dark proteins with conserved genomic context
defence-associated dark proteins
mobilization-associated dark proteins
dark proteins with known fold
dark proteins with no known structural match
rare recurrent dark proteins
lineage-specific dark proteins
small-protein candidates
membrane-associated candidates

The exact composition should be decided after inspecting the complete dataset.

37. Minimal End-to-End Workflow

143,504 plasmids
        │
        ▼
metadata normalization
        │
        ▼
Prodigal/Pyrodigal
        │
        ▼
ORF occurrences
        │
        ▼
protein dereplication
        │
        ▼
AntiFam
        │
        ├── AntiFam positive ──> flagged / discovery excluded
        │
        ▼
partial / <20 aa filtering for dark discovery
        │
        ▼
annotation hub
        ├── Pfam GA
        ├── Pfam relaxed
        ├── InterPro
        ├── eggNOG
        ├── Swiss-Prot
        ├── HMM-HMM
        └── broad DB
        │
        ▼
annotation normalization
        │
        ▼
annotation adjudication
        │
        ├── ANNOTATED
        └── DARK
                 │
                 ├── MMseqs2 family
                 ├── genomic context
                 ├── DefenseFinder
                 ├── IntegronFinder
                 ├── synteny/conservation
                 ├── evolutionary analysis
                 ├── structure/Foldseek
                 └── sequence properties
                         │
                         ▼
                 evidence integration
                         │
                         ▼
              complete annotation tables

38. Version 1 Deliverables

At minimum:

final/plasmids_complete.tsv
final/orf_occurrences_complete.tsv
final/protein_annotations_complete.tsv
final/annotation_complete.tsv
final/dark_families_complete.tsv
final/dark_family_members.tsv
final/plasmid_annotation_summary.tsv
final/normalization_dictionary.tsv
final/pipeline_qc_report.tsv

Optional raw/intermediate deliverables:

raw tool outputs
normalized tool outputs
family tables
context tables
evolution tables
structure tables

39. Definition of Success for Version 1

The pipeline is successful if:

every input plasmid is either fully processed or has an explicit failure record;

every predicted ORF has a stable identity;

every unique protein is linked back to all original plasmid occurrences;

AntiFam status is available;

partial and <20 aa records are flagged and recoverable;

annotation evidence from all enabled tools is preserved;

raw database labels remain available;

normalized annotation classes are available;

annotated/dark states are reproducible;

hypothetical proteins are not incorrectly treated as functionally characterized;

MMseqs2 family membership is available;

genomic context is available for dark ORFs;

evolutionary measurements are available where valid;

structural status is available where attempted;

sequence/physical properties are available;

distribution values are normalized or raw values are sufficient to perform normalization later;

evidence is named rather than collapsed into one final score;

complete tables can be regenerated from the workflow;

no downstream candidate ranking has contaminated the annotation layer;

every final row is auditable back to the original sequence and tool evidence.

40. Design Principle to Preserve

The most important property of the entire system is:

The pipeline should maximize the amount of interpretable information retained about every ORF while minimizing irreversible decisions.

The annotation pipeline should therefore end with a complete evidence universe, not a list of "interesting" proteins.

The top-1,000 experimental candidates are a separate downstream decision.