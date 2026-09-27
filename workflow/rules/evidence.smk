# =====================================================================================
# S5-S8 and the report: the evidence recorded for every dark family. Nothing here ranks or
# selects proteins; the 1,000 for experimental follow-up are chosen by hand from these
# tables.
#
# S5  target eligibility  unnamed, searched, not artefact-flagged
# S6  dark set, families  MMseqs2 deep-homology clustering, family network, lineages
# S7  evolutionary        recurrence, CDS recovery, codon alignments, dN/dS, RNAcode
# S8  context, structure  DefenseFinder, CONJScan, IntegronFinder, ISEScan, directons,
#                         context terms, Foldseek
# S9+ synteny, rarity, report
# =====================================================================================

rule clonal_registry:
    """S0b: per plasmid, its MOB-suite cluster, topology, observed host and predicted host
    range. Independent occurrences are counted over Stage 6 lineages (plasmid_lineage),
    not over these clusters."""
    input:
        master=config["input"]["master_table"],
        ids=f"{OUT}/01_analysis_set/analysis_set.txt",
        # Host names beyond PLSDB's (plasmidann.hosts).
        ps_hosts=config["input"]["host_provenance"],
        working_set=config["input"]["working_set"],
    output:
        f"{OUT}/01_analysis_set/clonal_registry.tsv",
    resources:
        mem_mb=4000,
        runtime=60,
    benchmark:
        f"{OUT}/benchmarks/clonal_registry.tsv"
    log:
        f"{OUT}/logs/01_analysis_set/clonal_registry.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/clonal_registry.py"


rule target_eligibility:
    """S5: which proteins are screening candidates: unnamed, searched, not artefacts."""
    input:
        prot=f"{OUT}/05_annotation_cascade/protein_annotation.tsv",
        artefact=f"{OUT}/04_orf_qc/artefact_flags.tsv",
    output:
        flags=f"{OUT}/09_target_eligibility/target_eligibility.tsv",
        report=f"{OUT}/09_target_eligibility/target_eligibility.txt",
    resources:
        # ~1.4 KB per protein row measured, ~4.9 GB at 3.5 M proteins.
        mem_mb=16000,
        runtime=60,
    benchmark:
        f"{OUT}/benchmarks/target_eligibility.tsv"
    log:
        f"{OUT}/logs/09_target_eligibility/target_eligibility.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/target_eligibility.py"


rule dark_set:
    """S6a: the proteins that are screening candidates at all."""
    input:
        flags=f"{OUT}/09_target_eligibility/target_eligibility.tsv",
        faa=f"{OUT}/03_dereplication/unique_proteins.faa",
        map=f"{OUT}/03_dereplication/protein_map.tsv",
        index=f"{OUT}/02_orf_calling/orf_index.tsv",
    output:
        faa=f"{OUT}/10_clustering/dark_proteins.faa",
        ids=f"{OUT}/10_clustering/dark_ids.txt",
    resources:
        mem_mb=16000,
        runtime=120,
    benchmark:
        f"{OUT}/benchmarks/dark_set.tsv"
    log:
        f"{OUT}/logs/10_clustering/dark_set.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/dark_set.py"


rule plasmid_lineage:
    """Stage 6: cluster plasmids by sequence similarity into independent lineages.

    Separate from MOB class: MOB typing describes the relaxase a plasmid carries and says
    nothing about whether two records are the same molecule sequenced twice.
    """
    input:
        fasta=f"{OUT}/01_analysis_set/analysis_set.fna",
    output:
        tsv=f"{OUT}/10_clustering/plasmid_lineage.tsv",
    params:
        lineage=targets["lineage"],
    threads: 16
    resources:
        mem_mb=32000,
        runtime=720,
    benchmark:
        f"{OUT}/benchmarks/plasmid_lineage.tsv"
    log:
        f"{OUT}/logs/10_clustering/plasmid_lineage.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/plasmid_lineage.py"


rule protein_families:
    """Stage 5: the family table for every unique protein, not only the dark set.

    The clusters are made before the cascade (S2f); this adds the annotation and the
    small/large scope. The distribution counts are in recurrence.tsv (Stage 7). dark_member_count, annotated_member_count,
    percentage_dark_in_family and the dark-only family (100% dark) need the annotated
    members present.

    Two outputs: the complete table across every configured resolution, and the derived
    dark-family subset at the primary resolution that the dark stages read.
    """
    input:
        # Made before the cascade (S2f, protein_clustering); annotation is added here.
        clusters=expand(f"{OUT}/10_clustering/families_{{res}}_cluster.tsv",
                        res=targets["clustering"]["resolutions"]),
        prot=f"{OUT}/05_annotation_cascade/protein_annotation.tsv",
        small_ids=f"{OUT}/01_analysis_set/small_plasmids.txt",
        dark_ids=f"{OUT}/10_clustering/dark_ids.txt",
        map=f"{OUT}/03_dereplication/protein_map.tsv",
        registry=f"{OUT}/01_analysis_set/clonal_registry.tsv",
    output:
        families=f"{OUT}/10_clustering/protein_families.tsv",
        dark_families=f"{OUT}/10_clustering/dark_families.tsv",
    params:
        clustering=targets["clustering"],
    resources:
        mem_mb=64000,
        runtime=1440,
    benchmark:
        f"{OUT}/benchmarks/protein_families.tsv"
    log:
        f"{OUT}/logs/10_clustering/protein_families.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/protein_families.py"


rule family_network:
    """Stage 5b: 50%-identity clusters linked by sequence similarity (Durairaj et al.
    2023), with communities and an annotation state per node - a map of where the dark
    plasmidome sits relative to the known. Reads the intermediate clustering Stage 5
    already made; changes no family and no dark call.
    """
    input:
        families=f"{OUT}/10_clustering/protein_families.tsv",
        reps=f"{OUT}/10_clustering/families_{targets['network']['node_resolution']}_rep_seq.fasta",
        clusters=f"{OUT}/10_clustering/families_{targets['network']['node_resolution']}_cluster.tsv",
        prot=f"{OUT}/05_annotation_cascade/protein_annotation.tsv",
        dark_ids=f"{OUT}/10_clustering/dark_ids.txt",
        map=f"{OUT}/03_dereplication/protein_map.tsv",
    output:
        nodes=f"{OUT}/10_clustering/network_nodes.tsv",
        edges=f"{OUT}/10_clustering/network_edges.tsv",
        summary=f"{OUT}/10_clustering/network_summary.tsv",
    params:
        network=targets["network"],
        primary=targets["clustering"]["primary"],
        node_resolution=targets["network"]["node_resolution"],
        seed=config["seed"],
    threads: workflow.cores
    resources:
        mem_mb=64000,
        runtime=1440,
    benchmark:
        f"{OUT}/benchmarks/family_network.tsv"
    log:
        f"{OUT}/logs/10_clustering/family_network.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/family_network.py"


rule rarity:
    """Stage 14: rarity labels per family, and the dark-family rarefaction curve.

    The curve answers whether the collection has saturated - whether more plasmids would
    keep revealing new dark families - which is what says if the dark count is a lower
    bound.
    """
    input:
        recurrence=f"{OUT}/11_distribution_and_evolution/recurrence.tsv",
        dark_families=f"{OUT}/10_clustering/dark_families.tsv",
        map=f"{OUT}/03_dereplication/protein_map.tsv",
        # The rarefaction axis: every small plasmid, with a dark family or without.
        small_ids=f"{OUT}/01_analysis_set/small_plasmids.txt",
    output:
        rarity=f"{OUT}/14_rarity/family_rarity.tsv",
        rarefaction=f"{OUT}/15_report/dark_family_rarefaction.tsv",
    params:
        rarity=targets["rarity"],
        seed=config["seed"],
    resources:
        mem_mb=16000,
        runtime=240,
    benchmark:
        f"{OUT}/benchmarks/rarity.tsv"
    log:
        f"{OUT}/logs/14_rarity/rarity.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/rarity.py"


rule synteny:
    """Stage 9: does a dark family's gene order RECUR across its occurrences?

    Distinct from Stage 8, which asks what one ORF sits next to once. Conserved gene order
    survives because the arrangement matters, so it is a much stronger claim than
    adjacency. Six measurements, kept separate, counted over Stage 6
    lineages, at every level of synteny.levels (close and intermediate).
    """
    input:
        annotation=f"{OUT}/06_annotation_tables/plasmid_annotation.tsv",
        # Only family_id is read: the primary-level set must equal it.
        families=f"{OUT}/10_clustering/dark_families.tsv",
        map=f"{OUT}/03_dereplication/protein_map.tsv",
        # One cluster file per level, in the order of synteny.levels: a gene and its
        # neighbours are named by their cluster at that level.
        clusters=expand(f"{OUT}/10_clustering/families_{{res}}_cluster.tsv",
                        res=targets["synteny"]["levels"]),
        # The counting unit: one vote per independent lineage (Stage 6).
        lineage=f"{OUT}/10_clustering/plasmid_lineage.tsv",
        dark_ids=f"{OUT}/10_clustering/dark_ids.txt",
        small_ids=f"{OUT}/01_analysis_set/small_plasmids.txt",
        # Topology: the neighbour window wraps across the origin of a circular plasmid.
        registry=f"{OUT}/01_analysis_set/clonal_registry.tsv",
        # Directons merge across the origin of a circular plasmid, which needs its length.
        lengths=f"{OUT}/01_analysis_set/plasmid_lengths.tsv",
    output:
        tsv=f"{OUT}/13_synteny/synteny.tsv",
    params:
        context=targets["context"],
        primary=targets["clustering"]["primary"],
        synteny=targets["synteny"],
    resources:
        mem_mb=32000,
        runtime=480,
    benchmark:
        f"{OUT}/benchmarks/synteny.tsv"
    log:
        f"{OUT}/logs/13_synteny/synteny.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/synteny.py"


rule dark_cooccurrence:
    """S8g: do two dark families travel together more often than chance predicts?

    Together = a member ORF of each on the same plasmid; counted once per Stage 6 lineage;
    tested by the hypergeometric upper tail over lineages, with Benjamini-Hochberg across
    the pairs together in at least cooccurrence.min_lineages_together lineages. Only those
    tested pairs are written. See plasmidann.cooccurrence.
    """
    input:
        # The primary-resolution dark families and their dark members.
        families=f"{OUT}/10_clustering/dark_families.tsv",
        map=f"{OUT}/03_dereplication/protein_map.tsv",
        lineage=f"{OUT}/10_clustering/plasmid_lineage.tsv",
    output:
        f"{OUT}/12_context_and_structure/dark_cooccurrence.tsv",
    params:
        cooccurrence=targets["cooccurrence"],
    resources:
        # Memory is ~100 B per distinct family pair sharing a plasmid (measured, a Counter
        # of tuple keys) plus ~1 KB per tested pair. results_test (job 6985208, 100
        # plasmids): 533 pair occurrences, 8.8 per small plasmid and 0.15 per large one,
        # 25 MB peak, under 1 s. Scaled to 82,261 small and 61,242 large plasmids that is
        # ~0.73 M pair occurrences (< 1 GB). The ceiling, if large plasmids carried as many
        # dark families as when every dark protein made a family (results_bench, 2,631
        # pairs per large plasmid), is ~162 M (~16 GB for the counts). 64 GB is a
        # scheduling figure above that ceiling, not a measurement at full scale.
        mem_mb=64000,
        # Each tested pair's tail takes 2-25 us (measured); ten million pairs take minutes.
        runtime=120,
    benchmark:
        f"{OUT}/benchmarks/dark_cooccurrence.tsv"
    log:
        f"{OUT}/logs/12_context_and_structure/dark_cooccurrence.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/dark_cooccurrence.py"


rule recurrence:
    """Stage 7: distribution and recurrence, counted over independent units.

    Seven counts per family, never collapsed, because database record counts are not
    independent biological observations.
    """
    input:
        families=f"{OUT}/10_clustering/protein_families.tsv",
        map=f"{OUT}/03_dereplication/protein_map.tsv",
        registry=f"{OUT}/01_analysis_set/clonal_registry.tsv",
        lineage=f"{OUT}/10_clustering/plasmid_lineage.tsv",
        master=config["input"]["master_table"],
    output:
        tsv=f"{OUT}/11_distribution_and_evolution/recurrence.tsv",
    resources:
        mem_mb=16000,
        runtime=240,
    benchmark:
        f"{OUT}/benchmarks/recurrence.tsv"
    log:
        f"{OUT}/logs/11_distribution_and_evolution/recurrence.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/recurrence.py"


rule extract_cds:
    """S7a: recover nucleotide CDS - dN/dS needs codons, and we store protein only."""
    input:
        ids=f"{OUT}/10_clustering/dark_ids.txt",
        map=f"{OUT}/03_dereplication/protein_map.tsv",
        index=f"{OUT}/02_orf_calling/orf_index.tsv",
        fasta=f"{OUT}/01_analysis_set/analysis_set.fna",
    output:
        f"{OUT}/11_distribution_and_evolution/dark_cds.fna",
    resources:
        mem_mb=16000,
        runtime=240,
    benchmark:
        f"{OUT}/benchmarks/extract_cds.tsv"
    log:
        f"{OUT}/logs/11_distribution_and_evolution/extract_cds.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/extract_cds.py"


rule family_evolution:
    """S7b: dN/dS per family - the strongest evidence a dark ORF is a real protein."""
    input:
        families=f"{OUT}/10_clustering/dark_families.tsv",
        faa=f"{OUT}/10_clustering/dark_proteins.faa",
        cds=f"{OUT}/11_distribution_and_evolution/dark_cds.fna",
    output:
        tsv=f"{OUT}/11_distribution_and_evolution/family_evolution.tsv",
        # The family consensus, built from the protein alignment this rule already makes.
        # S7c re-searches it: a family can be collectively recognisable while every member
        # individually misses the cut, and Pavlopoulos removed 6.5% of clusters that way.
        consensus=f"{OUT}/11_distribution_and_evolution/family_consensus.faa",
    params:
        evolution=targets["evolution"],
    threads: 16
    resources:
        mem_mb=16000,
        runtime=2880,
    benchmark:
        f"{OUT}/benchmarks/family_evolution.tsv"
    log:
        f"{OUT}/logs/11_distribution_and_evolution/family_evolution.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/family_evolution.py"


rule consensus_recheck:
    """S7c: is the family collectively novel, or only individually unmatched?

    A label, never a filter. A family whose consensus hits Pfam keeps its row and gains
    `collectively_novel = 0` plus the name of what it matched.
    """
    input:
        consensus=f"{OUT}/11_distribution_and_evolution/family_consensus.faa",
        # Every dark family gets a row; those without a consensus say NOT_RUN.
        families=f"{OUT}/10_clustering/dark_families.tsv",
    output:
        f"{OUT}/11_distribution_and_evolution/consensus_recheck.tsv",
    params:
        db=TIER_BY_ID[TIER_IDS[0]]["db"],
        args=TIER_BY_ID[TIER_IDS[0]]["args"],
        hmmer_z=cascade["hmmer_z"],
    threads: 8
    resources:
        mem_mb=16000,
        runtime=240,
    benchmark:
        f"{OUT}/benchmarks/consensus_recheck.tsv"
    log:
        f"{OUT}/logs/11_distribution_and_evolution/consensus_recheck.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/consensus_recheck.py"


rule defence_search:
    """S8a phase 1: which proteins look like defence components.

    Runs on the DEREPLICATED set with --db-type unordered, MacSyFinder's
    "components only, no system calling" mode. Safe to dereplicate here and only here:
    phase 1 asks a per-protein question that depends on identity alone.
    """
    input:
        faa=f"{OUT}/03_dereplication/unique_proteins.faa",
    output:
        tsv=f"{OUT}/12_context_and_structure/defence_components.tsv",
    params:
        models_dir=config["references"]["macsyfinder_models"],
        # When false and the models are absent, the stage records NOT_RUN rather than
        # halting: a missing optional database must not be fatal to a deliverable that
        # does not depend on it.
        required=targets["defence"]["required"],
    threads: 16
    resources:
        mem_mb=16000,
        runtime=720,
    benchmark:
        f"{OUT}/benchmarks/defence_search.tsv"
    log:
        f"{OUT}/logs/12_context_and_structure/defence_search.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/defence_search.py"


rule defence_gembase:
    """S8a phase 1.5: propagate component labels, prune, write genomic order.

    A component hit on a unique protein applies to every ORF sharing that sequence, so one
    search covers all copies. Plasmids carrying no component are pruned - they cannot meet
    any model's quorum. Measured on the test run this kept 41 of 100 plasmids but 83% of
    their ORFs, so phase 2 reads ~1.2x less.
    """
    input:
        components=f"{OUT}/12_context_and_structure/defence_components.tsv",
        map=f"{OUT}/03_dereplication/protein_map.tsv",
        index=f"{OUT}/02_orf_calling/orf_index.tsv",
    output:
        faa=f"{OUT}/12_context_and_structure/defence_candidates.faa",
        map=f"{OUT}/12_context_and_structure/defence_gembase_map.tsv",
    resources:
        mem_mb=24000,
        runtime=180,
    benchmark:
        f"{OUT}/benchmarks/defence_gembase.tsv"
    log:
        f"{OUT}/logs/12_context_and_structure/defence_gembase.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/defence_gembase.py"


rule defence_systems:
    """S8a phase 2: call systems from gene adjacency.

    MacSyFinder is driven directly rather than through `defense-finder run`, because the
    wrapper does not pass --replicon-topology through and 94% of these plasmids are
    circular. Under linear topology a system spanning the origin is invisible.
    """
    input:
        faa=f"{OUT}/12_context_and_structure/defence_candidates.faa",
        map=f"{OUT}/12_context_and_structure/defence_gembase_map.tsv",
    output:
        tsv=f"{OUT}/12_context_and_structure/defence_systems.tsv",
    params:
        models_dir=config["references"]["macsyfinder_models"],
        required=targets["defence"]["required"],
    # ONE MacSyFinder process over every candidate replicon, the cores to --worker: chunks
    # made the calls depend on -c, because HMMER's i-evalue scales with the database size
    # (defence_systems.py).
    threads: workflow.cores
    resources:
        mem_mb=48000,
        runtime=1440,
    benchmark:
        f"{OUT}/benchmarks/defence_systems.tsv"
    log:
        f"{OUT}/logs/12_context_and_structure/defence_systems.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/defence_systems.py"


rule conjugation_systems:
    """S8f: conjugation and mobilisation systems (CONJScan 2.1.0, Plasmids models) on
    every plasmid, and each plasmid's mobility class (pCONJ, pdCONJ, pMOB, pMOBless).

    Every ORF of every plasmid in genomic order, as ONE MacSyFinder database: HMMER's
    i-evalue scales with the database size, so per-core chunks made the calls depend on -c.
    MacSyFinder 2.1.6 runs from envs/conjscan, named by path (conjugation.exe).
    """
    input:
        index=f"{OUT}/02_orf_calling/orf_index.tsv",
        # Pre-flight has checked the executable against the models' grammar.
        preflight=f"{OUT}/05_annotation_cascade/preflight.tsv",
    output:
        systems=f"{OUT}/12_context_and_structure/conjugation_systems.tsv",
        classes=f"{OUT}/12_context_and_structure/conjugation_plasmid_class.tsv",
    params:
        models_dir=config["references"]["conjscan_models"],
        exe=targets["conjugation"]["exe"],
        required=targets["conjugation"]["required"],
        version=targets["conjugation"]["version"],
    threads: workflow.cores
    resources:
        # Extrapolated from the test set, not a production measurement: ~7 GB
        # and ~25 min on 16 workers at 9.3 M ORFs (~0.77 ms CPU per ORF).
        mem_mb=16000,
        runtime=240,
    benchmark:
        f"{OUT}/benchmarks/conjugation_systems.tsv"
    log:
        f"{OUT}/logs/12_context_and_structure/conjugation_systems.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/conjugation_systems.py"


rule integrons:
    """S8b: integron cassette arrays - the strongest plasmid-specific signal available.

    IntegronFinder walks the replicons one at a time and threads only its HMM searches,
    so the analysis set is split into one chunk per core, each run on one thread.
    """
    input:
        fasta=f"{OUT}/01_analysis_set/analysis_set.fna",
    output:
        f"{OUT}/12_context_and_structure/integrons.tsv",
    threads: workflow.cores
    resources:
        # ~0.3 GB per IntegronFinder process (test run), one per core.
        mem_mb=48000,
        runtime=4320,
    benchmark:
        f"{OUT}/benchmarks/integrons.tsv"
    log:
        f"{OUT}/logs/12_context_and_structure/integrons.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/integrons.py"


rule is_elements:
    """S8e: insertion sequence elements - boundaries, IS family, complete or partial.

    ISEScan on the whole analysis set in one job, split into one chunk per core with each
    chunk on one thread: its own threading kept 3.8 of 16 cores busy on the test run
    (1,693 CPU-s in 441 s, 5.4 Mbp). The analysis set is 8.7 Gbp, 1,607x that, so ~760
    core-hours, ~8 h on 96 cores if the chunks balance.
    """
    input:
        fasta=f"{OUT}/01_analysis_set/analysis_set.fna",
    output:
        f"{OUT}/12_context_and_structure/is_elements.tsv",
    threads: workflow.cores
    resources:
        # ~0.8 GB per ISEScan process (test run), one per core.
        mem_mb=96000,
        runtime=4320,
    benchmark:
        f"{OUT}/benchmarks/is_elements.tsv"
    log:
        f"{OUT}/logs/12_context_and_structure/is_elements.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/is_elements.py"


rule structure_search:
    """S8d: structural homology for the dark set, via Foldseek + ProstT5.

    Scope is family representatives by default: ProstT5 is a transformer, and the query
    count is the cost of this stage.
    """
    input:
        faa=f"{OUT}/10_clustering/dark_proteins.faa",
        families=f"{OUT}/10_clustering/dark_families.tsv",
    output:
        f"{OUT}/12_context_and_structure/structure_hits.tsv",
    params:
        structure=targets["structure"],
        target_db=config["foldseek_db"],
        prostt5=config["prostt5_model"],
    # preflight already failed the run if structure.required and the databases are absent,
    # so reaching this rule means either the databases exist or structure is optional.
    threads: 16
    resources:
        mem_mb=32000,
        runtime=1440,
    benchmark:
        f"{OUT}/benchmarks/structure_search.tsv"
    log:
        f"{OUT}/logs/12_context_and_structure/structure.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/structure_search.py"


rule context_features:
    """S8c: genomic context per ORF, as one row of descriptive rates per family (defence,
    conjugation, integron and IS element membership, annotated neighbours, operons), and
    the context terms per family counted over lineages (family_context_terms.tsv): the
    plasmid label databases and the defence and conjugation systems, never KEGG. Rows for
    the dark families and for every known family at the primary resolution, which is the
    benchmark tools/calibrate_context.py reads after the run. No enrichment test."""
    input:
        annotation=f"{OUT}/06_annotation_tables/plasmid_annotation.tsv",
        families=f"{OUT}/10_clustering/dark_families.tsv",
        map=f"{OUT}/03_dereplication/protein_map.tsv",
        defence=f"{OUT}/12_context_and_structure/defence_systems.tsv",
        conjugation=f"{OUT}/12_context_and_structure/conjugation_systems.tsv",
        # The terms: each protein's labels (kind, label, sub_label).
        labels=f"{OUT}/08_protein_labels/protein_labels.tsv",
        # Every family at the primary resolution, for the known-family benchmark rows.
        all_families=f"{OUT}/10_clustering/protein_families.tsv",
        lineage=f"{OUT}/10_clustering/plasmid_lineage.tsv",
        integrons=f"{OUT}/12_context_and_structure/integrons.tsv",
        is_elements=f"{OUT}/12_context_and_structure/is_elements.tsv",
        master=config["input"]["master_table"],
        # Directons merge across the origin of a circular plasmid, which needs its length.
        lengths=f"{OUT}/01_analysis_set/plasmid_lengths.tsv",
    output:
        families=f"{OUT}/12_context_and_structure/family_context.tsv",
        terms=f"{OUT}/12_context_and_structure/family_context_terms.tsv",
    params:
        context=targets["context"],
        primary=targets["clustering"]["primary"],
    resources:
        # ~2.3 KB per ORF measured, ~20 GB at 8.3 M ORFs.
        mem_mb=48000,
        runtime=480,
    benchmark:
        f"{OUT}/benchmarks/context_features.tsv"
    log:
        f"{OUT}/logs/12_context_and_structure/context.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/context_features.py"


rule annotation_report:
    """The deliverable: every annotation the run produced, as CSV.

    Not a shortlist. One row per ORF and one row per dark family, with every piece of
    evidence side by side, so that choosing candidates is a decision made ON this table
    rather than one baked into a rule. Nothing is filtered and nothing is ranked.
    """
    input:
        annotation=f"{OUT}/06_annotation_tables/plasmid_annotation.tsv",
        map=f"{OUT}/03_dereplication/protein_map.tsv",
        orthology=f"{OUT}/07_orthology/orthology.tsv",
        families=f"{OUT}/10_clustering/dark_families.tsv",
        evolution=f"{OUT}/11_distribution_and_evolution/family_evolution.tsv",
        recheck=f"{OUT}/11_distribution_and_evolution/consensus_recheck.tsv",
        context=f"{OUT}/12_context_and_structure/family_context.tsv",
        structure=f"{OUT}/12_context_and_structure/structure_hits.tsv",
        # Stages 7, 9 and 14: distribution counts, synteny and rarity labels.
        recurrence=f"{OUT}/11_distribution_and_evolution/recurrence.tsv",
        synteny=f"{OUT}/13_synteny/synteny.tsv",
        rarity=f"{OUT}/14_rarity/family_rarity.tsv",
        is_elements=f"{OUT}/12_context_and_structure/is_elements.tsv",
        registry=f"{OUT}/01_analysis_set/clonal_registry.tsv",
        # Per ORF: the close-level synteny row of the ORF's close cluster.
        clusters_close=f"{OUT}/10_clustering/families_close_cluster.tsv",
        # Per ORF: the plasmid label databases (S4d) and the CONJScan calls (S8f).
        labels_plasmid=f"{OUT}/08_protein_labels/protein_labels_plasmid.tsv",
        conjugation=f"{OUT}/12_context_and_structure/conjugation_systems.tsv",
        conjugation_class=f"{OUT}/12_context_and_structure/conjugation_plasmid_class.tsv",
        # Per family: partners it travels with (S8g).
        cooccurrence=f"{OUT}/12_context_and_structure/dark_cooccurrence.tsv",
    output:
        annotation=f"{OUT}/15_report/annotation_complete.csv",
        families=f"{OUT}/15_report/dark_families_complete.csv",
    params:
        evolution=targets["evolution"],
        cooccurrence=targets["cooccurrence"],
    resources:
        mem_mb=32000,
        runtime=240,
    benchmark:
        f"{OUT}/benchmarks/annotation_report.tsv"
    log:
        f"{OUT}/logs/15_report/annotation_report.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/annotation_report.py"
