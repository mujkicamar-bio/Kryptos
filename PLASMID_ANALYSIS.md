Dark ORF Discovery Pipeline

Agent Implementation Specification — Final Revised Design

Status: Pre-implementation specification
Primary purpose: Comprehensive discovery-oriented annotation of plasmid-encoded ORFs
Experimental objective: Generate a complete, auditable population of putative dark proteins from which a later analysis can select approximately 1,000 candidates for functional assays in Escherichia coli.

1. Scientific Objective

The project analyzes a large plasmid sequence collection with the goal of identifying plasmid-encoded proteins for which no established biological function is currently available.

The central scientific question is:

Which plasmid-encoded proteins remain functionally unresolved after a broad, sensitive, evidence-preserving annotation workflow, and what biological information can be extracted about them without prematurely filtering away novel biology?

The pipeline has three conceptually separate layers.

Layer A — Annotation

What do we know about this ORF?

Layer B — Dark-protein characterization

What evidence exists that this unresolved protein is real, conserved, biologically associated, or structurally interpretable?

Layer C — Experimental prioritization

Which unresolved proteins should ultimately be tested experimentally?

Layer C is not part of the core annotation pipeline.

The core pipeline produces a complete evidence universe from which Layer C operates later.

This separation is a fundamental project requirement.

2. Design Principles

2.1 Preserve information

The pipeline should maximize information retained about every ORF while minimizing irreversible decisions.

Unusual proteins must not be removed merely because they are:

short;

membrane-associated;

low-complexity;

cysteine-rich;

highly conserved;

restricted to one lineage;

structurally unusual.

These properties should generally be represented as labels or measurements.

2.2 Sensitivity before prioritization

The annotation pipeline is a discovery pipeline.

It should favor broad detection and evidence preservation.

The eventual experimental ranking belongs to a separate downstream workflow.

2.3 No global biological score

Do not collapse all evidence into a single:

dark_score
novelty_score
candidate_score
biological_importance_score

Evidence must remain named and interpretable.

2.4 Search thresholds and interpretation thresholds are different

A computational search threshold controls whether additional analysis is performed.

An interpretation threshold controls whether an observed result is sufficiently strong to affect annotation.

A reporting threshold controls whether a result becomes the primary annotation.

These are separate concepts.

2.5 Cost-aware cascade

Do not run every annotation method against every protein merely because the method exists.

Use the validated cascade to reduce the search population for deeper/expensive tiers.

The current project configuration uses:

narrow_at = 0.9
min_explained = 0.5

unless later benchmarked evidence leads to a formal revision.

2.6 Recurrence is not independence

A protein appearing many times in sequence databases does not mean that the observations are biologically independent.

The pipeline must distinguish:

database recurrence
unique plasmid recurrence
independent plasmid sequence clusters
host breadth
taxonomic breadth
MOB breadth

2.7 Context is evidence, not direct annotation

A dark ORF repeatedly located next to a defense system may be:

defence_associated

It must not automatically be called:

defence_protein

2.8 Structure is an evidence dimension

Structural similarity can provide a useful biological hypothesis even when sequence homology does not produce a clear functional interpretation.

Structural and sequence evidence are related and must not be treated as statistically independent observations merely because different tools were used.

2.9 Missing values retain meaning

A test that:

was not run;

was not applicable;

lacked enough sequences;

saturated;

failed;

must not be represented as the same state.

Numeric values are nullable; companion status fields retain meaning.

3. Dataset

3.1 Primary plasmid collection

The analysis begins from the PlasmidScope plasmid collection.

Current approximate scale:

~143,504 plasmids
~9.3 million ORF occurrences
~3.5 million unique protein sequences

The final analysis-set denominator must be generated from the actual configured input set and recorded in the run manifest.

The study's focus is the plasmids below 20 kb (82,261 of the 143,503; section 13.3). All
plasmids are processed; the small ones decide what the cascade annotates.

Do not hard-code the total plasmid count in downstream code.

3.2 Metadata

Metadata may be obtained from:

PlasmidScope

PLSDB

NCBI

IMG/PR

MOB-suite

other explicitly configured sources

Typical fields:

plasmid_id
accession
sequence_source
source_database
length
topology
host
host_taxonomy
habitat
MOB_class
MOB_cluster
plasmid_type
BioProject
submitter/source provenance

Raw values must remain available.

4. Data Architecture

The implementation must use a normalized relational logical model.

The giant integrated occurrence table is a view/export, not the authoritative storage model.

4.1 Authoritative tables

plasmids
orf_occurrences
proteins
annotations
annotation_hits
families
family_members
plasmid_lineages
context_features
context_occurrences
evolution
structure
properties
normalization_dictionary
evidence
controls
run_manifest

4.2 Storage

Primary implementation:

Flat TSV per stage, CSV for the two integrated deliverables.

A columnar store was specified here once and never wired into any stage. Measured on the
validation run it compressed no better than gzip, so the declaration was removed rather
than implemented, and the storage model is now what the code does.

Do not require loading the complete occurrence dataset into pandas.

This requirement stands and is not satisfied by the file format alone. Projected at the
full corpus, annotation_complete is ~2.6 GB over ~8.1M rows, which pandas would expand
several times over in memory. Readers query it out of core instead, directly over the
delimited text, plain or gzipped, with no conversion step.

4.3 Information levels

Plasmid-level

plasmid_length
host
habitat
MOB_cluster
plasmid_lineage_cluster

Occurrence-level

orf_occurrence_id
plasmid_id
start
end
strand
partial
neighbors

Protein-level

protein_id
protein_sequence
CDS_nucleotide_sequence
annotation_status
dark_status
annotation evidence
structure
properties

Family-level

family_id
family_size
dark_member_count
percentage_dark_in_family
family_plasmid_count
family_context_conservation

Do not silently duplicate family-level facts into every occurrence row.

5. Identifier System

Identifiers must remain stable when the biological input is unchanged.

5.1 Plasmid ID

Use a source accession where available.

Otherwise use a deterministic content-derived identifier.

5.2 ORF occurrence ID

Do not use an ordinal number.

Use:

<plasmid_id>:<start>-<end>:<strand>

The implementation must correctly represent circular origin-spanning ORFs.

5.3 Protein ID

Use a deterministic sequence hash:

protein_id = SHA256(normalized_protein_sequence)[:32]

This must be mandatory.

5.4 Family ID

Family IDs must not depend on result ordering.

Recommended:

<family_resolution>:<representative_protein_id>

or a deterministic hash of the sorted member set.

The selected convention must be fixed in configuration.

5.5 Observation IDs

Annotation/context/structure observations should use content-derived IDs.

Example:

SHA256(
    protein_id,
    tool,
    database,
    database_version,
    target_accession,
    query_start,
    query_end
)

Tool version remains provenance but does not automatically change the observation ID.

6. Run Manifest

Every production run must create:

run_manifest.yaml

Required:

pipeline_version
git_commit
input_dataset_id
input_dataset_hash
plasmid_count
orf_count
unique_protein_count
configuration_file
database_versions
database_paths
software_versions
SLURM configuration
thread counts
date

The run manifest is an input dependency for all resumable stages.

7. Missing-Value Semantics

7.1 Numeric fields

Numeric columns are nullable.

Example:

dnds_value = NULL
dnds_status = SATURATED

7.2 Status fields

Use explicit status values such as:

NOT_RUN
NO_HIT
TOO_FEW_MEMBERS
NO_DIVERGENCE
SATURATED
NO_OUTPUT
FAILED
NOT_APPLICABLE
SUCCESS

7.3 Text fields

Use:

NA

for genuinely unavailable text values.

Do not use NA as a substitute for a meaningful status.

8. Stage 1 — ORF Prediction

8.1 Tool

Use:

Prodigal / Pyrodigal

Version 1 does not introduce a second gene caller.

8.2 Inputs

Plasmid nucleotide FASTA.

8.3 Required fields

plasmid_id
orf_occurrence_id
protein_id_pre_derep
start
end
strand
nucleotide_length
protein_length
protein_sequence
CDS_nucleotide_sequence
partial
start_type

8.3b OPEN ISSUE (2026-09-24) — translation table

Meta mode chooses a gene model per plasmid, and some models use translation table 4
(TGA read as Trp), which is correct only for Mycoplasma/Spiroplasma. On 2,000 random
analysis-set plasmids, 22 (1.1%) were called with table 4, four of them at 41-46% GC; on
the test plasmid GenBank_CP124069.1, 9 of 64 genes contain an internal TGA read as Trp.
Such genes can read through a real stop, match nothing, and enter the dark set.

Until this is decided, the table is recorded on every ORF (translation_table: 11 or 4,
carried to plasmid_annotation.tsv and annotation_complete.csv) and orf_call reports the
count of table-4 genes. Options under consideration: force table 11; allow table 4 only
for Mollicute hosts; keep and flag. Downstream stages (extract_cds, dN/dS, RNAcode)
assume table 11.

8.4 Circular plasmids

Origin-spanning ORFs must be explicitly supported.

The CDS nucleotide sequence must be reconstructed across the origin.

Do not assume:

end > start

for every ORF.

DECIDED 2026-09-24: terminal repeats. A circular record whose first and last >= 20 bp are
identical (orf.min_terminal_repeat_bp; CheckV's direct-terminal-repeat criterion, Nayfach
et al. 2021) still carries the assembler's overlap: 400 of 400 sampled 'direct terminal
repeat' records do, and in 76% the repeat is not a multiple of 3 long, so joining the ends
through both copies shifts the frame of every gene across the junction. S0 writes such a
record with the last copy removed and lists it in 01_analysis_set/terminal_repeats.tsv
(record_bp, repeat_bp, molecule_bp). Every later stage, gene calling included, reads the
molecule once; coordinates of origin-spanning genes wrap at molecule_bp, not size_bp.

9. Stage 2 — ORF QC and Artifact Screening

Version 1 intentionally keeps the explicit artifact tool set limited to:

Prodigal/Pyrodigal-derived information
AntiFam

No large collection of additional artifact detectors is introduced in this version.

However, inexpensive ORF-QC measurements are retained because they are valuable evidence.

9.1 AntiFam

Run AntiFam against dereplicated protein sequences.

Use the established project threshold:

--cut_ga

Do not substitute an arbitrary universal E-value.

9.2 AntiFam fields

antifam_hit
antifam_model
antifam_description
antifam_score
antifam_evalue

9.3 Overlap QC

Calculate:

opposite_strand_overlap_bp
opposite_strand_overlap_fraction
opposite_strand_overlapping_orf_id
same_strand_overlap_fraction

These are QC features, not automatic exclusion criteria except where explicitly defined below.

Opposite-strand overlap is especially important because a shadow ORF can be recurrent and conserved while representing the translation of an unrelated real gene.

9.4 Origin-spanning QC

Record:

spans_origin

and retain the correctly reconstructed CDS.

10. Discovery Eligibility

The primary dark-protein discovery population excludes:

AntiFam-positive ORFs
partial-only proteins
proteins <20 aa

These records remain in the master dataset.

They are flagged rather than deleted.

11. Minimum Protein Length

The agreed discovery cutoff is:

<20 amino acids

11.1 Configuration

orf:
    min_call_length_aa: 20

discovery:
    min_dark_candidate_length_aa: 20

11.2 Rationale

The cutoff is an operational discovery boundary.

Proteins below 20 aa are increasingly difficult to distinguish from accidental ORFs using conventional:

homology;

family clustering;

profile evidence;

structure prediction.

However, genuine bacterial microproteins can occur below this length.

Therefore:

discovery_excluded_short = TRUE

rather than deleting the sequence from the dataset.

A future dedicated microprotein analysis can revisit this population.

12. Partial CDS Handling

Partial status is occurrence-level information.

Example:

Protein X
    Plasmid A = partial
    Plasmid B = complete
    Plasmid C = complete

The protein remains eligible if it has at least one complete eligible occurrence.

Required fields:

partial_occurrence_count
complete_occurrence_count
has_complete_occurrence

A protein must not become globally partial merely because one occurrence is incomplete.

13. Stage 3 — Dereplication

Dereplication is only a computational optimization.

It must not change biological occurrence counts.

13.1 Unique protein table

protein_id
protein_sequence
protein_length

13.2 Occurrence mapping

protein_id
orf_occurrence_id
plasmid_id

Every occurrence must remain recoverable.

13.3 ADDED 2026-09-24. Small-plasmid focus: families first, then what the cascade searches

The study is about small plasmids: size_bp below input.max_plasmid_size_bp (20,000 bp, the
antimode of this collection's size distribution; docs/PARAMETER_PROVENANCE.md). Every
plasmid is still processed - gene calling, Tier 0, families, lineages, defence, integrons,
IS elements - because the large plasmids are where a small-plasmid family's relatives live.
S0 writes 01_analysis_set/small_plasmids.txt; "on a small plasmid" means on at least one of
those, everywhere downstream.

Three steps decide what the cascade (Stage 4) spends its time on:

S2f  protein_clustering   every unique protein, three resolutions (close 90/80,
                          intermediate 50/80, broad 30/50; cov-mode 1, cluster-mode 2).
                          A family is a sequence cluster, so it is made before any
                          annotation and no member is lost to an earlier filter.
S2s  cascade_selection    SELECT the proteins Tier 0 does not annotate, in families
                          holding a small-plasmid protein Tier 0 does not annotate. Then
                          SEARCH only the representatives of a 90% clustering of those
                          (90% identity, 80% coverage of both; UniRef90).
S4r  cascade_resolve      each 90% member takes its representative's result
                          (annot_source = representative, annot_representative); a protein
                          outside every selected family is NOT_SEARCHED, which is neither
                          dark nor annotated.

The family level only chooses; no annotation is copied between proteins below 90% identity.
Because every unexplained member of a selected family is annotated by the same
cascade, "known" and "unknown" mean the same on the small and the large side of every such
family.

Measured on PlasmidScope's own proteins (2026-09-24): 1,239,766 unique proteins without a
Tier 0 annotation, 148,284 of them on small plasmids; 356,959 selected; 216,546 searched.
nr, the dominant cost, is then about 93 hours on one 96-core node.

DECIDED 2026-09-25: a protein family is the INTERMEDIATE clustering (50% identity, 80%
coverage; clustering.primary). Selection, dark families, synteny and the report's rarity
labels all use it; close and broad are computed and reported, never used as the family.
The counts above were measured with broad families; intermediate families are smaller, so
fewer proteins are selected and the nr estimate is an upper bound until re-measured.

Family table (Stage 5) additions: n_small_members, n_large_members, n_not_searched,
family_small_plasmid_count, scope (small_only_known | small_only_unknown | mixed_known |
mixed_unknown) and known_from (small | large | both). A family is written only when it holds
a small-plasmid protein; a dark family is one with a dark small-plasmid member, and carries
all its dark members (members) and the small-plasmid ones (small_members).

Synteny (Stage 9) and evolution (S7b) report every measurement twice: over all occurrences
or members, and with the prefix small_ over those on small plasmids. Synteny compares
neighbours by family, not by label: every gene has a family, dark genes included, and
it means the same on small and large plasmids.

The 1,000 proteins for experimental follow-up are chosen by hand from these tables; the
pipeline does not select them.

14. Stage 4 — Annotation Cascade

The annotation architecture is a cost-aware cascade.

It is not a full all-against-all annotation hub.

The production cascade is:

PlasmidScope eggNOG (precomputed, identical sequence only)
    ↓
Pfam GA
    ↓
Pfam relaxed
    ↓
pharokka
    ↓
Swiss-Prot
    ↓
eggNOG
    ↓
nr / broad database

The exact ordering may be adjusted after benchmarking, but every change must be documented.

No HMM-HMM or InterPro stage is included in the production pipeline.

ADDED 2026-09-23. Tier 0 — PlasmidScope eggNOG.

PlasmidScope (Li et al., Nucleic Acids Res. 2025, 53:D179) published eggNOG-mapper
2.1.12 results for every protein of its ALL set. A protein whose sequence is identical to a
PlasmidScope protein (same seq_id) takes that result; the ORFs remain our own Stage 1 calls.

A protein is resolved at Tier 0 when PlasmidScope gives it a KEGG KO, an EC number, or a
Pfam family whose name passes the cascade's own is_informative test. A DUF or UPF family
alone, or a COG/OG with only a category letter, names no function and does not resolve a
protein. It is then:

FUNCTIONAL, annot_tier PS, annot_completeness NOT_MEASURED

and it is not searched by any later tier. Every other protein - PlasmidScope-dark, or
absent from PlasmidScope - is annotated by the full cascade below, nr included, when its
family holds an unexplained small-plasmid protein (13.3).

Tier 0 does not change the HMMER search space: hmmer_z (-Z) counts every unique protein
plus the controls, searched or not. E-values therefore mean what they would if the whole
collection had been searched, do not depend on PlasmidScope's coverage, and stay on the
conservative side; a -Z counting only the searched proteins halved every E-value on the
test set and let weak hits remove proteins from the dark set.

The PlasmidScope product name is not used: its ORFs are Prodigal calls for some source
databases and deposited PGAP/submitter CDS for others, so only the eggNOG fields mean the
same thing on every row.

Known cost: eggNOG reports no alignment span, so a Tier 0 protein has no explained
fraction. A protein that the cascade alone would have left target-eligible (NONE or
UNCHARACTERIZED_HOMOLOG) but that PlasmidScope annotates leaves the dark set; on the
100-plasmid test set, with the rule above, the dark set fell from 2,142 to 1,947 proteins,
193 of them for this reason, measured before nr (docs/plasmidscope_enrichment.md). The
other two gained a pharokka label: pharokka searches its profiles against the proteins it
is given, so its E-values depend on the T3 input size, which -Z does not control.

15. Cascade Control

Two concepts must remain separate:

narrow_at
min_explained

15.1 narrow_at

Controls computational stopping.

A protein that is sufficiently covered by a current tier may avoid deeper searches.

15.2 min_explained

Controls scientific classification.

The current project value:

min_explained = 0.5

is explicitly the threshold on merged informative annotation coverage.

The exact calculation of informative coverage must be implemented as a named function and tested.

15.3 Current cascade value

narrow_at = 0.9

These parameters must be stamped into the run manifest and relevant output rows.

16. Tier 1 — Pfam GA

Run:

hmmsearch
Pfam-A
--cut_ga

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

Pfam evidence is strong domain/family evidence but does not automatically prove a precise molecular function.

17. Tier 2 — Pfam Relaxed

Run on proteins not sufficiently explained by Tier 1.

Purpose:

remote_or_partial_domain_candidate

Required fields:

pfam_relaxed_hit
pfam_relaxed_evalue
pfam_relaxed_score
pfam_relaxed_query_coverage
pfam_relaxed_alignment_start
pfam_relaxed_alignment_end

The relaxed threshold must be explicitly configured.

18. Tier 3 — pharokka (phage protein families)

A phage tier is included because phage/plasmid-associated proteins may be biologically relevant to the discovery problem.

REVISED 2026-09-18. The tier is executed by pharokka in protein mode, which annotates
against the prokaryotic virus protein families as a MMseqs2 sequence database and as
HMMER3 profiles searched with pyhmmer, and carries a functional category for every family.
The families' own distribution is HH-suite HHM profiles, which cannot be searched without
the HMM-HMM tooling section 79 excludes; an MMseqs2 conversion of those profiles and an
HH-suite build were both implemented and abandoned, and nothing of either remains in the
pipeline. pharokka is the tool; the families are what it searches.

Run pharokka in protein mode against its configured database.

Required:

pharokka_hit
pharokka_family_id
pharokka_annotation
pharokka_category
pharokka_evalue

No query coverage or alignment span is required, because pharokka does not report one: it
deletes its raw alignment tables on exit. A hit is therefore a FAMILY-LEVEL assignment -
the families are whole-protein clusters - and is treated as such: FUNCTIONAL when the
family is named, with annot_completeness recorded as NOT_MEASURED rather than derived from
a span that was never observed.

pharokka belongs to the broader sequence/domain homology evidence dimension.

It is not independent evidence from Pfam merely because it uses a different database.

19. Tier 4 — DIAMOND Swiss-Prot

Use the configured Swiss-Prot database.

Required:

diamond_swissprot_target
diamond_swissprot_description
diamond_swissprot_bitscore
diamond_swissprot_evalue
diamond_swissprot_identity
diamond_swissprot_query_coverage
diamond_swissprot_subject_coverage
diamond_swissprot_alignment_length

Where available, also retain target-level functional provenance such as:

protein_existence
annotation_evidence
curation/provenance

A high-scoring hit to an uncharacterized target is not itself functional characterization.

20. Tier 5 — eggNOG-mapper

Run eggNOG-mapper according to the configured cascade strategy.

Record:

eggnog_ortholog_group
eggnog_description
eggnog_best_og
eggnog_function
eggnog_evalue
eggnog_score

Automated orthology transfer must not automatically outrank direct experimentally supported/curated functional evidence.

21. Tier 6 — DIAMOND Broad Database

Use the configured broad protein database, e.g. nr.

DECIDED 2026-09-24: the broad database is NCBI ClusteredNR (release 2026-08-30; nr clustered
at 90% identity and 90% length, 545.6M representatives, 173.6G residues), replacing the
2025-03-03 nr (707M sequences, 273G residues). It is newer, so fewer proteins are called
dark only because their relatives were deposited after the snapshot, and it is about 0.64x
the residues to search. A hit is to a cluster representative and carries its title. Built
by tools/download_clustered_nr.sh; version in references.nr_version.

Purpose:

BROAD_HOMOLOGY_DISCOVERY

Required:

diamond_nr_target
diamond_nr_description
diamond_nr_bitscore
diamond_nr_evalue
diamond_nr_identity
diamond_nr_query_coverage
diamond_nr_subject_coverage
diamond_nr_alignment_length

Broad-database descriptions remain provenance information.

DECIDED 2026-09-24: a protein that Pfam (Tiers 1-2) or Swiss-Prot (Tier 4) named - any
informative hit - is not searched against the broad database (skip_if_named_by in
config/cascade.yaml). The broad database is spent only on what the curated databases could
not name, and cannot replace a curated label with free text. The unexplained part of a
DOMAIN_ONLY protein that Pfam named is therefore not searched there. Every DIAMOND hit
records identity, alignment length, bit score and target length (hits.tsv), as required
above and in section 19.

Do not automatically turn:

hypothetical protein

into functional annotation.

22. Removed Annotation Methods

The following are deliberately absent from the production pipeline:

HMM-HMM / HH-suite
InterPro / InterProScan

Rationale:

computational cost at the current dataset scale;

substantial overlap with existing domain/homology evidence;

they are not required for the initial definition of the dark protein population;

they can be reconsidered in a future revision after a benchmark.

No agent may reintroduce either method without a formal specification revision.

23. Annotation Evidence Retention

The pipeline should retain:

tool
database
database_version
target
target_description
score
evalue
identity
query_coverage
subject_coverage
alignment_length
alignment_coordinates
hit_rank

Raw tool output retention must follow the configured storage policy.

The implementation must not silently discard raw outputs merely because only the best hit is used in final interpretation.

24. Annotation Normalization

Normalization occurs after annotation results exist.

24.1 Raw values

Always preserve:

raw_annotation
raw_database
raw_description

24.2 Controlled annotation classes

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

The normalization dictionary must be versioned.

25. Hypothetical and Uncharacterized Proteins

The following are considered dark when no stronger functional evidence exists:

hypothetical protein
uncharacterized protein
conserved hypothetical protein
protein of unknown function
DUF / unknown-function family

This is a core project decision.

A protein may be:

well established as real
+
highly recurrent
+
highly conserved
+
functionally unresolved
=
DARK

This distinction is central to the project.

26. Top-Level Annotation Status

There are only two top-level states:

ANNOTATED
DARK

26.1 ANNOTATED

The complete configured annotation workflow provides sufficiently strong evidence for a biologically interpretable assignment under the configured min_explained rule.

26.2 DARK

The protein remains functionally unresolved after the configured annotation cascade.

This includes uncharacterized/hypothetical/DUF proteins.

There is no separate terminal UNCHARACTERIZED state.

27. Dark Subclasses

Dark proteins may receive descriptive labels:

DARK_HYPOTHETICAL
DARK_DUF
DARK_NO_SEQUENCE_FUNCTION
DARK_FOLD_KNOWN
DARK_NO_STRUCTURE

These labels do not form an ordinal ranking.

28. Operational Definition of DARK

A protein is:

DARK

when:

it is in the eligible discovery population;

it is not AntiFam-positive;

at least one complete occurrence exists;

protein length is at least 20 aa;

final merged informative annotation coverage is below:

min_explained = 0.5

after the configured annotation cascade has completed.

This includes proteins whose only homologues are themselves unnamed
(UNCHARACTERIZED_HOMOLOG). FESNov (Rodriguez del Rio et al. 2024) calls a family unknown
only when it has no homologue at all; the quality gate reports both counts
(target_eligible and dark_no_homologue).

The exact informative-coverage calculation must be defined in the adjudication module and covered by regression tests.

29. Primary Annotation

Every protein may have:

primary_annotation
primary_annotation_source
primary_annotation_confidence

only when the configured criteria are satisfied.

Suggested authority order:

experimentally supported/curated functional assignment
    ↓
strong curated family/domain assignment
    ↓
informative orthology assignment
    ↓
strong curated sequence homology
    ↓
weaker sequence evidence
    ↓
unknown/uncharacterized

This is an authority hierarchy, not a numerical score.

IMPLEMENTED 2026-09-24: annot_label comes from the most authoritative tier that named the
protein, in cascade order (Pfam GA, Pfam relaxed, pharokka, Swiss-Prot, broad database),
and within a tier from the hit with the best E-value. The label with the best E-value
across all tiers is kept beside it (best_evalue_label, best_evalue_tier). A protein whose
only informative names are domain-level - "X domain-containing protein", "X family
protein", which NCBI PGAP assigns from domain or family models (Li W. et al. 2021) - is
DOMAIN_ONLY whatever their coverage (named_by_domain_only = 1).

30. Annotation Transfer Safeguards

Do not perform uncontrolled transitive annotation.

Example:

A = known function
B = strong homolog of A
C = weak homolog of B

does not automatically imply:

C = function of A

Configured transfer criteria must include:

identity
coverage
E-value
alignment quality
target evidence quality
domain consistency

31. Stage 5 — MMseqs2 Protein Families

Construct protein families after the main annotation layer.

Family membership is a descriptive sequence relationship.

It does not automatically establish function.

31.1 Required fields

family_id
family_resolution

31.2 Family statistics

family_size
dark_member_count
annotated_member_count
percentage_dark_in_family

31.3 Distribution

family_plasmid_count
family_host_count
family_species_count
family_genus_count
family_MOB_count
family_plasmid_lineage_count
family_habitat_count

31.4 Multiple resolutions

The implementation should support configured resolutions such as:

close
intermediate
broad

Thresholds must be explicit.

31b. Stage 5b — Family Network (ADDED 2026-09-24)

A map of the protein families, built as Durairaj et al. (Nature 2023, 622:646) built
theirs over UniRef50, so the numbers can be compared with theirs.

Nodes: the intermediate clusters (50% identity, 80% coverage), the UniRef50 analogue.
Edges: MMseqs2 all-against-all over the cluster representatives; an alignment covering
at least 50% of either protein at E < 1e-4; at most four outbound edges per node; edge
weight from the E-value. MMseqs2 search sensitivity is left at its default, which Durairaj
et al. do not state.

Node attributes: members, ORFs, plasmids, brightness (the annotation coverage the
best-annotated member reaches; a FUNCTIONAL protein with an unmeasured span counts as 1),
dark (brightness <= 0.05), dark fraction, family, most frequent informative label,
community (asynchronous label propagation, seeded), degree.

Summary: dark nodes connected, dark nodes connected to a bright node, and broad-resolution
dark singletons with an edge - the measurement that decides whether the family thresholds
leave related proteins apart.

The network is description: it changes no family and no dark call.

32. Dark Families

A dark-only family is defined as:

percentage_dark_in_family = 100%

Dark-family status is descriptive.

It is not itself proof of biological importance.

33. Stage 6 — Plasmid Lineage / Independence

MOB classification and plasmid sequence similarity are separate concepts.

Retain:

MOB_cluster
plasmid_lineage_cluster

33.1 Plasmid sequence clustering

Use a configured sequence-based approach such as Mash/sketch or an appropriate ANI-like clustering method.

Generate:

plasmid_lineage_cluster

33.2 Interpretation

Example:

plasmid_occurrences = 2,143
MOB_clusters = 8
independent_plasmid_clusters = 47

These values represent different biological properties.

MOB breadth does not automatically equal evolutionary independence.

34. Stage 7 — Distribution and Recurrence

For each dark protein/family calculate:

plasmid_occurrence_count
unique_plasmid_count
independent_plasmid_cluster_count
host_count
species_count
genus_count
MOB_count
habitat_count

The host of each plasmid comes from three sources in order of preference: PLSDB species,
PlasmidScope's per-record host (which carries IMG/PR's), and the GenBank/RefSeq source
organism (plasmidann.hosts). Together they name a host for 61.1% of the analysis set (98.6%
of isolate plasmids, 24.0% of metagenomic ones) against 31.9% for PLSDB alone. Every
family also reports n_plasmids_with_host, and host_count_status is NOT_MEASURED rather
than a count of 0 when none of its plasmids has a recorded host.

DECIDED 2026-09-25: MOB-suite's predicted host range (mob_host_range) is reported beside the
host, as a separate measurement over every plasmid, hosted or not
(n_plasmids_with_predicted_range, predicted_host_range_count, predicted_host_ranges; per ORF
predicted_host_range, beside host_species and host_genus). It is a prediction at any rank - genus
to several phyla - and never enters host_count, genus_count or the CROSS_HOST and
CROSS_TAXON labels. It adds a range for 18,170 of the 55,785 plasmids without a host.

34.1 Database recurrence

Retain, where possible:

database_record_count
hypothetical_record_count
database_source_count

34.2 Provenance-aware recurrence

Where possible distinguish:

unique accession
unique upstream source
unique BioProject
unique plasmid lineage cluster

Database record counts must never be treated as independent biological observations.

35. Stage 8 — Genomic Context

The pipeline builds a context representation for all ORFs so that dark/background comparisons are possible.

35.1 Compact all-ORF context

Maintain a compact neighborhood representation covering approximately:

±10 ORFs

for all ORFs.

This provides the background topology without generating an unnecessarily large detailed table.

35.2 Detailed dark context

For dark ORFs, create detailed per-neighbor records for:

±3 ORFs

where:

one row = one dark occurrence × one neighbor

This balances information retention with dataset size.

36. Neighbor Fields

For detailed context:

neighbor_orf_id
distance_in_genes
distance_in_bp
strand
relative_orientation
neighbor_annotation
normalized_annotation_class
functional_role

Do not store the entire neighborhood as a single free-text field.

37. Context Functional Roles

Use a controlled vocabulary:

RELAXASE
MOBILIZATION
REPLICATION
PARTITION
TOXIN_ANTITOXIN
DEFENSE
INTEGRON
TRANSPOSITION
CONJUGATION
PLASMID_BACKBONE
TRANSCRIPTIONAL_REGULATION
MEMBRANE
METABOLIC
HYPOTHETICAL
OTHER

The source of each role must be retained.

38. DefenseFinder

Run DefenseFinder as configured.

Record:

defence_system
defence_system_type
system_coordinates
dark_orf_distance
occurrence_id
source

Generate:

defence_associated

as a contextual hypothesis.

Do not automatically assign the defense function to the dark ORF.

39. IntegronFinder

Run IntegronFinder as configured.

Record:

integron_id
integron_type
integron_coordinates
attC_present
dark_orf_distance

Generate:

integron_associated

as a contextual label.

39b. ISEScan (ADDED 2026-09-24)

Run ISEScan (Xie & Tang, Bioinformatics 2017, 33:3340) on the analysis set with its
published defaults, so partial elements are reported as well as complete ones.

Record:

is_id
is_family
is_cluster
is_coordinates
is_complete
tir

ISEScan calls its own genes; elements are joined to our ORFs by coordinates, as integron
arrays are.

Generate:

is_associated

as a contextual label (the IS element is an island in Stage 8 context), and, per ORF, the
IS families of the elements it lies inside.

A dark ORF inside an IS element is often a degenerate transposase fragment or an IS-borne
accessory gene; one beside it may be a passenger. Neither is assigned as the ORF's
function, and neither removes it from the dark set: flag, never discard.

40. Other Context Labels

Possible labels:

mobilization_associated
toxin_antitoxin_associated
replication_associated
partition_associated
backbone_associated

Every positive label must link to underlying evidence.

41. Context State

Every context feature supports:

TRUE
FALSE
NOT_ASSESSED

Missing output must never be silently interpreted as FALSE.

42. Stage 9 — Synteny and Context Conservation

For dark families with multiple occurrences compare neighborhoods.

Example:

A - B - DARK - C - D
A - B - DARK - C - D
A - B - DARK - C
A - B - DARK - C - D

Calculate:

left_neighbor_conservation
right_neighbor_conservation
neighborhood_conservation
operon_like_conservation
synteny_conservation
context_recurrence

These measurements remain separate from annotation.

Left and right are upstream and downstream on the dark gene's own strand, so one
arrangement written in the two orientations counts once. On a circular plasmid the
neighbour window wraps across the origin. operon_like means sharing a directon with a
FUNCTIONAL partner, the definition Stage 8 uses.

43. Stage 10 — Evolutionary Analysis

Evolutionary analyses are conditional.

Not every dark protein/family will have sufficient information for every test.

43.1 Distribution

plasmid_count
independent_plasmid_cluster_count
host_count
species_count
genus_count
MOB_count

43.2 Sequence conservation

Where appropriate:

mean_pairwise_identity
identity_distribution
alignment_coverage
conserved_positions
conserved_motifs
amino_acid_conservation
nucleotide_conservation

44. RNAcode

RNAcode is used where suitable homologous nucleotide sequences exist.

The pipeline must retain nucleotide CDS sequences from Stage 1.

44.1 Both-strand analysis

Run/record RNAcode for both relevant strands where the comparison is applicable.

Retain:

rnacode_support
rnacode_score
rnacode_pvalue
rnacode_coordinates
rnacode_strand
rnacode_status

The strand-specific analysis is evidence/QC.

It is not an automatic ORF deletion rule.

44.2 Multiple testing

Any large-scale significance analysis must apply an appropriate multiple-testing correction.

45. dN/dS

dN/dS is supporting evidence, not a universal requirement.

Run only where explicit configured criteria are satisfied.

Potential criteria:

minimum_sequence_count
minimum_alignment_coverage
alignment_quality
sequence_divergence
codon_alignment_quality

Required status:

NOT_RUN
TOO_FEW_MEMBERS
NO_DIVERGENCE
SATURATED
NO_OUTPUT
SUCCESS
FAILED

Do not force dN/dS onto families where it is not interpretable.

46. Stage 11 — Structural Analysis

Structural analysis is a separate evidence dimension.

46.1 Structural model state

NO_MODEL
MODEL_AVAILABLE

46.2 Structural relationship

NO_STRUCTURAL_MATCH
KNOWN_STRUCTURAL_RELATIONSHIP
NOT_ASSESSED

These are separate fields.

47. ProstT5 / 3Di

Use the configured ProstT5 workflow for structural representation.

Do not create a pLDDT field unless a separate structure-prediction system actually supplies such a metric.

Required:

structure_model_id
structure_model_source
structure_status

48. Foldseek

Search structural representations against the configured structural database.

Record:

foldseek_hit
foldseek_target
foldseek_description
foldseek_score
foldseek_evalue_if_available
foldseek_alignment_coverage
foldseek_tm_score_if_available

A convincing relationship results in:

KNOWN_STRUCTURAL_RELATIONSHIP

No convincing relationship results in:

NO_STRUCTURAL_MATCH

Do not automatically call this a novel fold.

49. Structural Compute Strategy

Because structure generation is expensive:

Discovery scale

Analyze dark-family representatives where practical.

Candidate scale

Perform detailed individual structural analysis for candidates entering downstream prioritization.

Targeted exceptions

Unusual individual dark proteins may be escalated earlier.

50. Stage 12 — Protein Properties

Calculate:

protein_length
molecular_weight
isoelectric_point
hydrophobicity
transmembrane_helix_count
signal_peptide
low_complexity_fraction
disorder_fraction
cysteine_fraction
coiled_coil_prediction
secondary_structure_summary

These are descriptors, not automatic filters.

51. Liability Labels

Possible:

SHORT
MEMBRANE
SIGNAL_PEPTIDE
LOW_COMPLEXITY
DISORDERED
CYS_RICH
UNCERTAIN_ORF
OTHER

Liabilities are retained for downstream experimental planning.

52. Stage 13 — Background Normalization

Distribution and context analyses require explicit backgrounds.

At minimum consider:

plasmid_length
gene_count
host_taxonomy
MOB_group
plasmid_lineage_cluster
habitat
source_database
clonal redundancy

52.1 Raw values

Always retain:

raw_count
raw_prevalence

52.2 Normalized values

When normalized values are generated, record:

normalization_method
normalization_version
background_definition

No normalized prevalence field is valid without its normalization definition.

DECIDED 2026-09-25: the pipeline produces no normalized values. The stratified background
(banded plasmid_length and gene_count) was built for the context enrichment of section 53
and was removed with it.

53. Context Enrichment

For example:

dark ORFs associated with defense systems

must be compared against an appropriate background of non-dark ORFs or another explicitly justified reference population.

Potential approaches:

size-matched background
stratified background
permutation background
gene-level normalization
plasmid-level normalization

The exact statistical model must be configured and documented.

DECIDED 2026-09-25: no context enrichment is computed. The Fisher test, its stratified
background and the label-category layer were removed. S8c (context_features) reports only
descriptive per-family rates - the fraction of a family's plasmids on which a member lies
in a defence system, an integron or an IS element, has an annotated neighbour, shares a
directon with an annotated gene, or is in a two-gene directon (cons_* columns). They are
not corrected for plasmid size. Neighbour identity per family is described by synteny
(Stage 9).

54. Stage 14 — Rarity and Conservation

Rarity and conservation remain separate descriptors.

Potential labels:

RARE
LINEAGE_SPECIFIC
WIDELY_CONSERVED
SINGLE_MOB
CROSS_MOB
SINGLE_HOST
CROSS_HOST
CROSS_TAXON

Thresholds must be explicit.

These labels are not experimental rankings.

DECIDED 2026-09-25: lineage breadth - the rarity labels and the multi_lineage reality test
(section 56) - is counted in Stage 6 Mash lineages only. MOB-suite assigns every plasmid the
cluster of its nearest reference however distant (53% of plasmids lie beyond its 0.06
threshold; cluster AA379 holds 21,817 unrelated small plasmids of the analysis set), so it
is reported only as SINGLE_MOB (one cluster) or CROSS_MOB (two or more);
PLASMID_FAMILY_SPECIFIC was removed. Host labels use observed hosts only: CROSS_HOST at two
or more species, CROSS_TAXON at two or more genera, SINGLE_HOST only when every plasmid of
the family is named to one species. MOB-suite's predicted host range is a separate
measurement, reported over all of a family's plasmids and per ORF, and never counted as a
host.

55. Rarefaction

Approximately 1,000 candidates is an experimental-budget objective, not a biological assumption.

The annotation pipeline must calculate a rarefaction curve of dark-family discovery versus plasmids sampled.

The plasmids sampled are all small plasmids, those with no dark family included; a family
is discovered on a plasmid through its dark small-plasmid members.

Deliver:

final/dark_family_rarefaction.tsv
final/dark_family_rarefaction.pdf

This estimates whether additional plasmids continue to reveal new dark families.

56. Stage 15 — Evidence Integration

Evidence integration combines results without producing an experimental ranking.

56.1 Evidence dimensions

Use named dimensions:

ORF_QC
SEQUENCE_HOMOLOGY
ORTHOLOGY
GENOMIC_CONTEXT
EVOLUTIONARY_CONSERVATION
DISTRIBUTION
STRUCTURAL_RELATIONSHIP
PROTEIN_PROPERTIES

56.2 Dependencies

The pipeline must recognize that these are not statistically independent.

For example:

Pfam
pharokka
Swiss-Prot
eggNOG
nr
MMseqs2
Foldseek

share evolutionary information to varying degrees.

The objective is not to manufacture "independent evidence" but to preserve distinct measurements and avoid double-counting them.

56.3 Descriptive evidence summary

Possible:

evidence_dimensions_present
supporting_observations_count

These are descriptive, not scores.

57. Functional Hypothesis Layer

Possible hypotheses:

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

Each hypothesis must have traceable support.

Example:

functional_hypothesis = defence_associated

DefenseFinder_system = CBASS
distance = 2 genes
recurrent_occurrences = 17
synteny_conserved = TRUE

The protein remains:

DARK

unless direct sequence/domain evidence supports a specific function.

DECIDED 2026-09-25: no functional_hypothesis is produced. It was derived from the top
enriched context category (section 53), which no longer exists. GENOMIC_CONTEXT remains an
evidence dimension (section 56.1), present when the S8c rates or synteny were measured.

58. Negative Controls

The pipeline must include both positive and negative controls.

58.1 Positive controls

Known functionally characterized proteins should be inserted as:

record_class = positive_control

and must be recovered as functionally annotated.

58.2 Negative controls

Include appropriately constructed decoy sequences designed to test whether the dark-classification workflow behaves as expected for non-biological sequences.

Possible decoy categories:

shuffled controls
reverse-complement-derived controls
other sequence-derived decoys

The exact construction must be fixed and documented.

Expected behavior:

record_class = negative_control
no functional annotation
DARK and/or explicit decoy/QC status

The negative control must not be treated as an observed biological protein.

59. Quality Gates

The pipeline must stop if positive-control annotation recall falls below the configured validated threshold.

Positive controls test:

Can the pipeline still recover known biology?

Negative controls test:

Does the pipeline avoid treating decoy sequences as convincing biological annotations?

Neither control establishes the expected size of the dark population.

60. Manual Review

A manually reviewed subset of:

DARK
ANNOTATED

proteins should be used as a validation set.

The review protocol must specify:

sample size
sampling method
reviewer/role
review criteria
pass criterion
date/version

Manual review is validation, not candidate ranking.

61. Deterministic Sweep Cohort

Threshold sensitivity analyses must use deterministic cohort membership based on the protein identifier.

The same protein must belong to the same sweep cohort across reruns when the input protein universe is unchanged.

Do not use uncontrolled random sampling for reproducibility-critical threshold sweeps.

62. Data Products

62.1 Plasmids

final/plasmids.tsv

62.2 ORF occurrences

final/orf_occurrences.tsv

62.3 Proteins

final/proteins.tsv

62.4 Annotations

final/annotations.tsv
final/annotation_hits.tsv

62.5 Families

final/families.tsv
final/family_members.tsv

62.6 Plasmid lineages

final/plasmid_lineages.tsv

62.7 Context

final/context_features.tsv
final/context_occurrences.tsv

62.8 Evolution

final/evolution.tsv

62.9 Structure

final/structure.tsv

62.10 Properties

final/properties.tsv

62.11 Evidence

final/evidence.tsv

62.12 Normalization dictionary

final/normalization_dictionary.tsv

62.13 Controls

final/control_results.tsv

63. Main Integrated Annotation View

Generate:

final/annotation_complete.csv

This is the authoritative integrated view, and it is the primary storage model rather than
an export of one. One row per ORF, every piece of evidence side by side, nothing filtered
and nothing ranked.

64. Dark Family Output

Generate:

final/dark_families_complete.csv

One row per dark family.

Fields include:

family_id
family_resolution
family_size
dark_member_count
percentage_dark_in_family
annotated_member_count
plasmid_count
independent_plasmid_cluster_count
host_count
species_count
genus_count
MOB_count
lineage_count
habitat_count
family_annotation_summary
family_context_summary
family_structural_summary

65. Dark Family Membership

Generate:

final/dark_family_members.tsv

One row per:

family_id × protein_id

Required:

family_id
protein_id
dark_status
annotation_status

66. Plasmid Annotation Summary

Generate:

final/plasmid_annotation_summary.tsv
final/plasmid_annotation_summary.tsv

Possible fields:

plasmid_id
plasmid_length
gene_count
annotated_gene_count
dark_gene_count
dark_family_count
defence_system_count
integron_count
mobilization_features
MOB_cluster
plasmid_lineage_cluster
host
habitat
source_database

67. Occurrence-Level Output

The occurrence-level authoritative table should contain occurrence-specific fields.

plasmid_id
orf_occurrence_id
protein_id
start
end
strand
partial
protein_length
discovery_excluded_short
antifam_hit
spans_origin

Protein/family-level measurements should remain in their own tables.

68. Agent Architecture

Agent 1 — Input / Metadata

Responsibilities:

ingest plasmid sequences;

validate FASTA;

normalize accession IDs;

ingest metadata;

preserve provenance;

resolve analysis-set denominator.

Agent 2 — ORF Prediction

Responsibilities:

run Prodigal/Pyrodigal;

generate coordinates;

reconstruct CDS sequences;

handle origin-spanning ORFs;

generate occurrence IDs.

Agent 3 — ORF QC / AntiFam

Responsibilities:

run AntiFam;

identify partials;

flag <20 aa;

calculate overlap QC;

generate QC statuses.

Agent 4 — Dereplication

Responsibilities:

generate protein sequence hashes;

construct unique protein table;

maintain occurrence mapping;

verify losslessness.

Agent 5 — Annotation Cascade

Responsibilities:

Pfam GA;

Pfam relaxed;

pharokka;

Swiss-Prot;

eggNOG;

nr/broad database.

No HMM-HMM or InterPro stages are included in the production specification.

Agent 6 — Annotation Normalization

Responsibilities:

preserve raw annotations;

map descriptions to controlled classes;

identify hypothetical/uncharacterized states.

Agent 7 — Annotation Adjudication

Responsibilities:

calculate merged informative coverage;

apply min_explained;

assign ANNOTATED/DARK;

determine primary annotation;

enforce transfer safeguards.

Agent 8 — Family Construction

Responsibilities:

run MMseqs2;

produce deterministic family IDs;

calculate family statistics;

identify dark-only families.

Agent 9 — Plasmid Lineage

Responsibilities:

perform plasmid sequence clustering;

create lineage cluster identifiers;

distinguish lineage clusters from MOB groups.

Agent 10 — Context

Responsibilities:

create compact all-ORF context;

create detailed dark ±3 context;

annotate neighbors;

run DefenseFinder;

run IntegronFinder.

Agent 11 — Synteny

Responsibilities:

compare neighborhoods;

calculate context and synteny conservation.

Agent 12 — Evolution

Responsibilities:

construct family alignments;

conservation;

RNAcode both strands where applicable;

dN/dS where applicable;

maintain explicit statuses.

Agent 13 — Structure

Responsibilities:

generate structure representations;

run Foldseek;

assign structural relationship state.

Agent 14 — Properties

Responsibilities:

sequence descriptors;

transmembrane prediction;

signal peptide;

low complexity;

disorder;

cysteine richness;

liability labels.

Agent 15 — Background Normalization

Responsibilities:

construct appropriate backgrounds;

account for plasmid size;

account for gene count;

account for taxonomy;

account for clonal/lineage structure;

calculate normalized prevalence/enrichment.

Agent 16 — Evidence Integration

Responsibilities:

merge evidence tables;

produce integrated views;

preserve evidence dependencies;

generate hypothesis summaries.

Must not rank experimental candidates.

Agent 17 — Validation

Responsibilities:

identifier validation;

row-count validation;

sequence validation;

annotation validation;

family validation;

context validation;

control validation;

provenance validation.

69. Agent Implementation Rule

Agents are implementation modules, not necessarily 17 independent programs.

Prefer implementing them as:

Snakemake rules
Python modules/scripts
shared schemas
shared utility functions

rather than creating 17 disconnected applications.

Each stage must have an explicit input/output contract.

70. Compute Strategy

The production strategy is:

cheap broad annotation
        ↓
cascade narrowing
        ↓
deeper annotation of unresolved proteins
        ↓
family analysis
        ↓
context/evolution
        ↓
representative-scale structure
        ↓
complete evidence integration

Avoid full-population expensive stages unless benchmarked.

71. Benchmarking

Before the production run, benchmark at least:

Pfam GA
Pfam relaxed
pharokka
Swiss-Prot
eggNOG
nr
MMseqs2
DefenseFinder
IntegronFinder
RNAcode
dN/dS
ProstT5
Foldseek

For each benchmark record:

input_count
runtime
CPU_hours
peak_RAM
I/O
output_size

The nr benchmark is mandatory before finalizing the scheduling configuration.

72. One job per stage

No stage is sharded. The pipeline takes one input FASTA; S0 applies the analysis-set
exclusion to it once, and every later stage is a single job over one file, given every core
the run has.

Decided 2026-09-22. A search against a streamed database has a fixed cost per invocation -
DIAMOND reads and indexes the whole of nr each time it is run - and the measurement in
workflow/bench_nr.sbatch (2,000 queries, 16 threads, more than 12 hours) showed that cost
dominates the run. Sharding the query set multiplies it by the shard count. Resumability is
per stage: --rerun-incomplete keeps every stage that finished and rebuilds the one that did
not.

73. Resume Semantics

Every stage must write a completion record containing:

expected_row_count
actual_row_count
input_hash
config_hash

Incomplete files must never be mistaken for successful stage outputs.

74. Provenance

Every important observation must be traceable through:

source plasmid
ORF coordinates
protein sequence
tool
database
database version
parameters
raw result
normalization rule
adjudication rule

The final integrated dataset must support end-to-end provenance tracing.

75. Scientific Interpretation Rules

No Pfam

Does not mean:

novel protein

No Swiss-Prot

Does not mean:

unknown biology

Hypothetical protein

Means:

sequence observed
+
function unresolved

DARK

Means:

no sufficiently supported functional annotation after the configured cascade

Known fold

Means:

structural relationship detected

not:

function proven

No structural match

Means:

no convincing relationship detected under the configured structural workflow

not:

novel fold proven

Many database records

Means:

repeated sequence representation

not:

independent observations = record count

MOB breadth

Means:

distribution across mobility classifications

not:

independent evolutionary lineages

Dark family

Means:

family members remain functionally unresolved

not:

family is biologically important

76. Experimental Prioritization Boundary

The core pipeline ends with:

A complete evidence-rich representation of the plasmid protein universe and its dark subset.

It does not produce:

top_1,000
candidate_score
novelty_score
experimental_rank

The downstream prioritization workflow consumes the complete evidence tables.

77. Future Experimental Prioritization

The downstream workflow may use:

dark status
conservation
independent recurrence
plasmid lineage breadth
MOB breadth
genomic context
defense association
mobilization association
structural relationship
structural state
protein properties
family diversity
experimental feasibility

Potential candidate classes include:

highly conserved dark proteins
dark proteins with conserved synteny
defence-associated dark proteins
mobilization-associated dark proteins
dark proteins with known fold but unknown function
dark proteins with no known structural relationship
rare recurrent dark proteins
lineage-specific dark proteins
small proteins >=20 aa
membrane-associated dark proteins

These are downstream selection categories, not annotation classifications.

78. Final Production Workflow

PLASMIDS
    │
    ▼
metadata normalization
    │
    ▼
Prodigal/Pyrodigal
    │
    ├── ORF coordinates
    ├── protein sequence
    └── CDS nucleotide sequence
    │
    ▼
protein dereplication
    │
    ▼
ORF QC
    ├── AntiFam
    ├── partial
    ├── <20 aa
    ├── overlap
    └── origin flags
    │
    ▼
ANNOTATION CASCADE
    ├── Pfam GA
    ├── Pfam relaxed
    ├── pharokka
    ├── Swiss-Prot
    ├── eggNOG
    └── nr
    │
    ▼
annotation normalization
    │
    ▼
annotation adjudication
    │
    ├────────────────┐
    ▼                ▼
ANNOTATED           DARK
                       │
                       ▼
                 MMseqs2 families
                       │
             ┌─────────┼──────────┐
             ▼         ▼          ▼
       plasmid       context    evolution
        lineage         │          │
             │      DefenseFinder  │
             │      IntegronFinder │
             │      ±3 neighbors  │
             │      synteny        │
             │                   RNAcode
             │                   dN/dS
             │
             └─────────┬──────────┘
                       │
                       ▼
                   structure
                 representatives
                       │
                       ▼
                   properties
                       │
                       ▼
              background normalization
                       │
                       ▼
                evidence integration
                       │
                       ▼
              COMPLETE DATASET
                       │
                       ▼
           separate prioritization
                       │
                       ▼
          experimental candidates

79. Definition of Success

Version 1 is successful when:

every analysis-set plasmid is processed or has an explicit failure record;

every ORF has a stable occurrence identifier;

every ORF is linked to its source plasmid;

every unique protein has a deterministic sequence-derived ID;

CDS nucleotide sequences are retained;

AntiFam status is available;

overlap/origin/partial/length QC information is retained;

proteins <20 aa are flagged and retained;

all eligible proteins enter the configured annotation cascade;

narrow_at and min_explained are explicit and recorded;

ANNOTATED/DARK classification is reproducible;

hypothetical/uncharacterized/DUF proteins are DARK when unresolved;

raw annotation provenance is retained;

pharokka is included;

HMM-HMM and InterPro are absent unless the pipeline specification is formally revised;

MMseqs2 family membership is reproducible;

family IDs are stable;

plasmid sequence lineage clusters are available;

MOB classification is retained separately;

context is available for all ORFs at the compact level;

detailed dark context is available for ±3 neighbors;

context background comparisons are possible;

RNAcode can distinguish strands where applicable;

evolutionary analyses have explicit applicability statuses;

structural state and structural relationship are separate;

structure is performed at representative scale in the production run;

protein properties/liabilities are retained;

recurrence is provenance-aware;

missing values retain semantic status;

positive controls pass;

negative controls behave according to the configured expectation;

manual validation is documented;

rarefaction is generated;

final outputs are queryable without loading the entire dataset into pandas;

every final result is traceable to its source evidence;

the core pipeline produces no experimental ranking.

80. Central Project Principle

The annotation pipeline should not decide which unknown proteins are interesting. It should build the most complete and defensible representation possible of what is known, unknown, recurrent, conserved, contextualized, and structurally interpretable about every plasmid ORF.

The project deliberately preserves three different questions:

WHAT IS KNOWN?

WHAT IS DARK?

WHAT IS WORTH TESTING?

The first two are outputs of this pipeline.

The third is a separate downstream scientific decision.