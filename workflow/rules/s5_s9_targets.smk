# =====================================================================================
# S5-S9: from the annotated plasmidome to 1,000 screening candidates.
#
# S5  quality gate        positive control (run-halting) + backbone stop-list
# S6  dark set, families  MMseqs2 deep-homology clustering
# S7  evolutionary        CDS recovery, codon alignments, dN/dS
# S8  context, structure  DefenseFinder, IntegronFinder, directons, Foldseek
# S9  prioritisation      composite score, weight sweep, stratified portfolio
# S9b library design      codon-optimised synthesis order
#
# Rationale: plans/2026-09-10-pipeline-v2-design.md
# =====================================================================================

rule clonal_registry:
    """S0b: which plasmids count as independent observations.

    Every later count of independent occurrences is over MOB clusters, not raw plasmids.
    Without this, a family on forty plasmids may be one clone sequenced forty times.
    """
    input:
        master=config["input"]["master_table"],
        ids=f"{OUT}/s0/analysis_set.txt",
    output:
        f"{OUT}/s0/clonal_registry.tsv",
    resources:
        mem_mb=4000,
        runtime=60,
    log:
        f"{OUT}/logs/s0/clonal_registry.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/clonal_registry.py"


rule quality_gate:
    """S5: halts the run if known plasmid biology comes out dark (SC2)."""
    input:
        prot=f"{OUT}/s3/protein_annotation.tsv",
        artefact=f"{OUT}/s2b/artefact_flags.tsv",
    output:
        flags=f"{OUT}/s5/target_eligibility.tsv",
        report=f"{OUT}/s5/quality_gate.txt",
    params:
        gate=targets["quality_gate"],
    resources:
        mem_mb=8000,
        runtime=60,
    log:
        f"{OUT}/logs/s5/quality_gate.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/quality_gate.py"


rule dark_set:
    """S6a: the proteins that are screening candidates at all."""
    input:
        flags=f"{OUT}/s5/target_eligibility.tsv",
        faa=f"{OUT}/s2/unique_proteins.faa",
        map=f"{OUT}/s2/protein_map.tsv",
        index=f"{OUT}/s1/orf_index.tsv",
    output:
        faa=f"{OUT}/s6/dark_proteins.faa",
        ids=f"{OUT}/s6/dark_ids.txt",
    resources:
        mem_mb=16000,
        runtime=120,
    log:
        f"{OUT}/logs/s6/dark_set.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/dark_set.py"


rule protein_families:
    """Stage 5: cluster EVERY unique protein into families (spec section 31).

    Not only the dark set. Section 31.2 requires dark_member_count,
    annotated_member_count and percentage_dark_in_family, and section 32 derives a
    dark-only family as 100% dark - all four need the annotated members present.

    Two outputs: the complete table across every configured resolution, and the derived
    dark-family subset at the primary resolution that the dark stages read.
    """
    input:
        faa=f"{OUT}/s2/unique_proteins.faa",
        dark_ids=f"{OUT}/s6/dark_ids.txt",
        map=f"{OUT}/s2/protein_map.tsv",
        registry=f"{OUT}/s0/clonal_registry.tsv",
    output:
        families=f"{OUT}/s6/protein_families.tsv",
        dark_families=f"{OUT}/s6/dark_families.tsv",
    params:
        clustering=targets["clustering"],
    threads: 16
    resources:
        mem_mb=64000,
        runtime=1440,
    log:
        f"{OUT}/logs/s6/protein_families.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/protein_families.py"


rule extract_cds:
    """S7a: recover nucleotide CDS - dN/dS needs codons, and we store protein only."""
    input:
        ids=f"{OUT}/s6/dark_ids.txt",
        map=f"{OUT}/s2/protein_map.tsv",
        index=f"{OUT}/s1/orf_index.tsv",
        fasta=config["input"]["working_set_fasta"],
    output:
        f"{OUT}/s7/dark_cds.fna",
    resources:
        mem_mb=16000,
        runtime=240,
    log:
        f"{OUT}/logs/s7/extract_cds.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/extract_cds.py"


rule family_evolution:
    """S7b: dN/dS per family - the strongest evidence a dark ORF is a real protein."""
    input:
        families=f"{OUT}/s6/dark_families.tsv",
        faa=f"{OUT}/s6/dark_proteins.faa",
        cds=f"{OUT}/s7/dark_cds.fna",
    output:
        tsv=f"{OUT}/s7/family_evolution.tsv",
        # The family consensus, built from the protein alignment this rule already makes.
        # S7c re-searches it: a family can be collectively recognisable while every member
        # individually misses the cut, and Pavlopoulos removed 6.5% of clusters that way.
        consensus=f"{OUT}/s7/family_consensus.faa",
    params:
        evolution=targets["evolution"],
    threads: 16
    resources:
        mem_mb=16000,
        runtime=2880,
    log:
        f"{OUT}/logs/s7/family_evolution.log",
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
        consensus=f"{OUT}/s7/family_consensus.faa",
    output:
        f"{OUT}/s7/consensus_recheck.tsv",
    params:
        db=TIER_BY_ID[TIER_IDS[0]]["db"],
        args=TIER_BY_ID[TIER_IDS[0]]["args"],
        hmmer_z=cascade["hmmer_z"],
    threads: 8
    resources:
        mem_mb=16000,
        runtime=240,
    log:
        f"{OUT}/logs/s7/consensus_recheck.log",
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
        faa=f"{OUT}/s2/unique_proteins.faa",
    output:
        tsv=f"{OUT}/s8/defence_components.tsv",
    params:
        models_dir=config["references"]["macsyfinder_models"],
        # When false and the models are absent, the stage records NOT_RUN rather than
        # halting: a missing optional database must not be fatal to a deliverable that
        # does not depend on it (spec section 7.2).
        required=targets["defence"]["required"],
    threads: 16
    resources:
        mem_mb=16000,
        runtime=720,
    log:
        f"{OUT}/logs/s8/defence_search.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/defence_search.py"


rule defence_gembase:
    """S8a phase 1.5: propagate component labels, prune, write genomic order.

    A component hit on a unique protein applies to every ORF sharing that sequence, so one
    search covers all copies. Plasmids carrying no component are pruned - they cannot meet
    any model's quorum - which is roughly a 4-5x reduction in what phase 2 must read.
    """
    input:
        components=f"{OUT}/s8/defence_components.tsv",
        map=f"{OUT}/s2/protein_map.tsv",
        index=f"{OUT}/s1/orf_index.tsv",
    output:
        faa=f"{OUT}/s8/defence_candidates.faa",
        map=f"{OUT}/s8/defence_gembase_map.tsv",
    resources:
        mem_mb=24000,
        runtime=180,
    log:
        f"{OUT}/logs/s8/defence_gembase.log",
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
        faa=f"{OUT}/s8/defence_candidates.faa",
        map=f"{OUT}/s8/defence_gembase_map.tsv",
    output:
        tsv=f"{OUT}/s8/defence_systems.tsv",
    params:
        models_dir=config["references"]["macsyfinder_models"],
        required=targets["defence"]["required"],
    threads: 16
    resources:
        mem_mb=16000,
        runtime=1440,
    log:
        f"{OUT}/logs/s8/defence_systems.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/defence_systems.py"


rule integrons:
    """S8b: integron cassette arrays - the strongest plasmid-specific signal available."""
    input:
        fasta=lambda wc: SHARD_PATHS[wc.shard],
    output:
        f"{OUT}/s8/integrons/{{shard}}.tsv",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=240,
    log:
        f"{OUT}/logs/s8/integrons/{{shard}}.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/integrons.py"


rule structure_search:
    """S8d: structural homology for the dark set, via Foldseek + ProstT5.

    Scope is representatives by default. Spec section 49 sets that as the discovery-scale
    strategy and section 79 makes it a success criterion; ProstT5 is a transformer and the
    query count is the cost of this stage.
    """
    input:
        faa=f"{OUT}/s6/dark_proteins.faa",
        families=f"{OUT}/s6/dark_families.tsv",
    output:
        f"{OUT}/s8/structure_hits.tsv",
    params:
        structure=targets["structure"],
        target_db=config.get("foldseek_db", "data/refs/foldseek/pdb"),
        prostt5=config.get("prostt5_model", "data/refs/foldseek/prostt5"),
    # preflight already failed the run if structure.required and the databases are absent,
    # so reaching this rule means either the databases exist or structure is optional.
    threads: 16
    resources:
        mem_mb=32000,
        runtime=1440,
    log:
        f"{OUT}/logs/s8/structure.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/structure_search.py"


rule context_features:
    """S8c: genomic context per ORF, aggregated to families against a background."""
    input:
        annotation=f"{OUT}/s4/plasmid_annotation.tsv",
        families=f"{OUT}/s6/dark_families.tsv",
        map=f"{OUT}/s2/protein_map.tsv",
        defence=f"{OUT}/s8/defence_systems.tsv",
        integrons=expand(f"{OUT}/s8/integrons/{{shard}}.tsv", shard=SHARDS),
        labels=f"{OUT}/s4c/protein_labels.tsv",
    output:
        families=f"{OUT}/s8/family_context.tsv",
        background=f"{OUT}/s8/context_background.tsv",
    params:
        context=targets["context"],
        # null means no grouping: every (kind, label) is its own category. The biological
        # grouping is derived from results/s4c/protein_labels.tsv after a full annotation
        # run and enabled by pointing this at config/label_categories.yaml.
        categories=config["references"]["label_categories"],
    resources:
        mem_mb=32000,
        runtime=480,
    log:
        f"{OUT}/logs/s8/context.log",
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
        annotation=f"{OUT}/s4/plasmid_annotation.tsv",
        map=f"{OUT}/s2/protein_map.tsv",
        orthology=f"{OUT}/s4b/orthology.tsv",
        families=f"{OUT}/s6/dark_families.tsv",
        evolution=f"{OUT}/s7/family_evolution.tsv",
        recheck=f"{OUT}/s7/consensus_recheck.tsv",
        context=f"{OUT}/s8/family_context.tsv",
        structure=f"{OUT}/s8/structure_hits.tsv",
    output:
        annotation=f"{OUT}/final/annotation_complete.csv",
        families=f"{OUT}/final/dark_families_complete.csv",
    params:
        evolution=targets["evolution"],
    resources:
        mem_mb=32000,
        runtime=240,
    log:
        f"{OUT}/logs/final/annotation_report.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/annotation_report.py"
