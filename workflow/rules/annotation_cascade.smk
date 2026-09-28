# =====================================================================================
# The annotation cascade
#
# Each tier searches whatever the previous tier could not explain, and hands on the
# residue. Tier order is authority order.
#
# A protein stops being searched at the tier where its explained fraction reaches
# `narrow_at`.
# =====================================================================================

rule preflight:
    """Fail in seconds if a tool or database is missing, rather than after 45 hours.

    A dependency of every tier and of the artefact screen, so no search can start until
    every configured tool and database has been confirmed present.
    """
    output:
        f"{OUT}/05_annotation_cascade/preflight.tsv",
    params:
        tiers=cascade["tiers"],
        artefact=cascade["artefact_screen"],
        structure=targets["structure"],
        orthology=targets["orthology"],
        foldseek_db=config["foldseek_db"],
        prostt5=config["prostt5_model"],
        # label_databases and conjugation_systems: tools named by path in their own
        # environments, and their databases.
        labels=config["labels"],
        amrfinder=config["amrfinder"],
        conjugation=targets["conjugation"],
        conjscan_models=config["references"]["conjscan_models"],
        defence=targets["defence"],
        macsyfinder_models=config["references"]["macsyfinder_models"],
        genomad=config["genomad"],
    conda:
        "../envs/plasmidann.yaml"
    resources:
        mem_mb=2000,
        # Minutes. `macsyfinder --version` from envs/conjscan took between 1.4 s and
        # 2 min 53 s on /gorilla (slow Python imports).
        runtime=90,
    benchmark:
        f"{OUT}/benchmarks/preflight.tsv"
    log:
        f"{OUT}/logs/05_annotation_cascade/preflight.log",
    script:
        "../scripts/preflight.py"


rule check_hmmer_z:
    """Refuse a declared -Z that does not describe the protein set (unique_proteins.faa),
    before the artefact screen or any tier searches."""
    input:
        unique=f"{OUT}/03_dereplication/unique_proteins.faa",
    output:
        f"{OUT}/03_dereplication/hmmer_z_checked.tsv",
    params:
        hmmer_z=cascade["hmmer_z"],
    conda:
        "../envs/plasmidann.yaml"
    resources:
        mem_mb=2000,
        runtime=30,
    log:
        f"{OUT}/logs/03_dereplication/check_hmmer_z.log",
    script:
        "../scripts/check_hmmer_z.py"


rule tier_search:
    """Search one tier, then hand the next tier only what stayed unexplained.

    One job per tier, with every core the run has. A search against a streamed database
    - DIAMOND reads the whole of nr per invocation - has a per-pass cost that sharding the
    query set multiplies by the shard count, which is why the cascade is not sharded.
    """
    input:
        faa=tier_query,
        spans=tier_spans,
        # Earlier tiers' hits, where this tier skips what they named (skip_if_named_by).
        named=tier_named,
        # The DIAMOND tiers skip artefact-flagged proteins.
        artefact=f"{OUT}/04_orf_qc/artefact_flags.tsv",
        preflight=f"{OUT}/05_annotation_cascade/preflight.tsv",
    output:
        hits=f"{OUT}/05_annotation_cascade/{{tier}}/hits.tsv",
        unresolved=f"{OUT}/05_annotation_cascade/{{tier}}/unresolved.faa",
        spans=f"{OUT}/05_annotation_cascade/{{tier}}/spans.tsv",
    params:
        spec=lambda wc: TIER_BY_ID[wc.tier],
        narrow_at=cascade["narrow_at"],
        hmmer_z=cascade["hmmer_z"],
        max_target_seqs=cascade["max_target_seqs"],
    threads: workflow.cores
    conda:
        "../envs/plasmidann.yaml"
    resources:
        # PER TIER, and MEASURED rather than declared, on full nr (T5 searches ClusteredNR,
        # which has not been measured).
        #
        #   nr     136 GB at 152 queries, 185 GB at 3,808 (bench_nr, job 6928733) and 180 GB
        #          at 3,910 in the pipeline. DIAMOND holds one database block (-b 16), but the
        #          peak also grows with the query count, so the production peak at ~200,000
        #          queries is NOT known. 220 GB is therefore a scheduling figure, not a bound:
        #          nothing else runs beside the tier (threads: workflow.cores), so it has the
        #          allocation's 550 GB in practice.
        #   others   8 GB measured (pharokka, the largest); the HMMER tiers peak at 0.35 GB
        #            because hmmsearch memory-maps the pressed library rather than loading it.
        mem_mb=lambda wc: 220000 if TIER_BY_ID[wc.tier]["source"] == "nr" else 24000,
        # nr on full nr: ~45 min per pass and ~1.7 s per query on 96 threads
        # standalone (bench_nr, -k 5); ~40 min and ~1.9 s inside the pipeline
        # (results_bench). Production uses -k 25, which was not benchmarked. Local
        # execution does not enforce runtime; the figure documents the expectation.
        runtime=lambda wc: 5760 if TIER_BY_ID[wc.tier]["source"] == "nr" else 1440,
    benchmark:
        f"{OUT}/benchmarks/tier_search.{{tier}}.tsv"
    log:
        f"{OUT}/logs/05_annotation_cascade/{{tier}}.log",
    script:
        "../scripts/tier_search.py"


rule cascade_resolve:
    """Assign functional_class per protein from every tier's hits."""
    input:
        hits=expand(f"{OUT}/05_annotation_cascade/{{tier}}/hits.tsv", tier=TIER_IDS),
        # The last tier's spans carry the cumulative explained fraction for every protein
        # the cascade ever saw, because each tier writes forward everything it inherited.
        spans=f"{OUT}/05_annotation_cascade/{TIER_IDS[-1]}/spans.tsv",
        faa=f"{OUT}/03_dereplication/search_representatives.faa",
        # The proteins that skipped the cascade, so the table covers every protein.
        ps=f"{OUT}/03_dereplication/plasmidscope_proteins.tsv",
        # Search-cluster members, and the proteins no selected family holds.
        selection=f"{OUT}/03_dereplication/selection.tsv",
    output:
        f"{OUT}/05_annotation_cascade/protein_annotation.tsv",
    params:
        thresholds=cascade,
        tier_order=TIER_IDS,
    resources:
        mem_mb=32000,
        runtime=480,
    benchmark:
        f"{OUT}/benchmarks/cascade_resolve.tsv"
    log:
        f"{OUT}/logs/05_annotation_cascade/resolve.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/cascade_resolve.py"


rule annotate_plasmids:
    """S4: the primary deliverable - every ORF with its annotation and provenance."""
    input:
        prot=f"{OUT}/05_annotation_cascade/protein_annotation.tsv",
        index=f"{OUT}/02_orf_calling/orf_index.tsv",
        map=f"{OUT}/03_dereplication/protein_map.tsv",
        artefact=f"{OUT}/04_orf_qc/artefact_flags.tsv",
    output:
        f"{OUT}/06_annotation_tables/plasmid_annotation.tsv",
    resources:
        mem_mb=24000,
        runtime=240,
    benchmark:
        f"{OUT}/benchmarks/annotate_plasmids.tsv"
    log:
        f"{OUT}/logs/06_annotation_tables/annotate.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/annotate_plasmids.py"


rule orthology:
    """S4b: COG and KEGG terms (eggNOG-mapper) for the proteins the cascade named, so
    that a dark ORF's neighbours can be described in aggregatable terms."""
    input:
        prot=f"{OUT}/05_annotation_cascade/protein_annotation.tsv",
        faa=f"{OUT}/03_dereplication/unique_proteins.faa",
        # PlasmidScope's own eggNOG result for the proteins it annotated.
        ps=f"{OUT}/03_dereplication/plasmidscope_proteins.tsv",
    output:
        f"{OUT}/07_orthology/orthology.tsv",
    params:
        orthology=targets["orthology"],
    threads: 16
    resources:
        mem_mb=32000,
        runtime=2880,
    benchmark:
        f"{OUT}/benchmarks/orthology.tsv"
    log:
        f"{OUT}/logs/07_orthology/orthology.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/orthology.py"


rule feature_files:
    """S4: GFF3 and GenBank of the analysis set; origin-spanning genes are written in
    each format's own way (plasmidann.features)."""
    input:
        annotation=f"{OUT}/06_annotation_tables/plasmid_annotation.tsv",
        fasta=f"{OUT}/01_analysis_set/analysis_set.fna",
        master=config["input"]["master_table"],
    output:
        gff3=f"{OUT}/06_annotation_tables/plasmid_annotation.gff3",
        genbank=f"{OUT}/06_annotation_tables/plasmid_annotation.gbk",
    resources:
        mem_mb=32000,
        runtime=480,
    benchmark:
        f"{OUT}/benchmarks/feature_files.tsv"
    log:
        f"{OUT}/logs/06_annotation_tables/feature_files.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/feature_files.py"


rule label_databases:
    """S4d: every unique protein against the plasmid-specific label databases.

    TADB, BacMet, oriTDB, mobileOG-db, dbAPIS and Anti-CRISPRdb by DIAMOND at PlasAnn's
    identity and coverage tiers (thresholds from the paper; PlasAnn's database and labels
    are not used), CARD protein homolog models by Perfect / Strict, and AMRFinderPlus by
    its own rules. One long table that protein_labels merges; a database that was not
    searched is NOT_RUN in the status table, never an empty result. See plasmidann.labeldb.
    """
    input:
        faa=f"{OUT}/03_dereplication/unique_proteins.faa",
        # No search starts before every tool and database has been confirmed present.
        preflight=f"{OUT}/05_annotation_cascade/preflight.tsv",
    output:
        tsv=f"{OUT}/08_protein_labels/protein_labels_plasmid.tsv",
        status=f"{OUT}/08_protein_labels/label_databases_status.tsv",
    params:
        labels=config["labels"],
        amrfinder=config["amrfinder"],
    # DIAMOND and AMRFinderPlus both thread; a third of the allocation lets the stage run
    # beside the cascade rather than queue behind it.
    threads: 32
    resources:
        # 0.61 GB peak on the test set (64 s wall, 260 CPU-s at 8 threads). DIAMOND's
        # memory is set by its block size, not by the query count, so 32 GB is a
        # scheduling figure with headroom, not a measurement at 3.5 M proteins.
        mem_mb=32000,
        # Extrapolated from the test set, linear in the query count: ~0.05 CPU-s per
        # protein, ~50 CPU-h at 3.5 M proteins, ~2 h on 32 threads. 24 h is the ceiling.
        runtime=1440,
    benchmark:
        f"{OUT}/benchmarks/label_databases.tsv"
    log:
        f"{OUT}/logs/08_protein_labels/label_databases.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/label_databases.py"


rule protein_labels:
    """Every functional label every tool produced, one long table, no role assigned."""
    input:
        hits=expand(f"{OUT}/05_annotation_cascade/{{tier}}/hits.tsv", tier=TIER_IDS),
        orthology=f"{OUT}/07_orthology/orthology.tsv",
        pfam_dat=config["references"]["pfam_dat"],
        # Search-cluster members take their representative's labels.
        selection=f"{OUT}/03_dereplication/selection.tsv",
        # The plasmid label databases (label_databases), merged in as their own kinds.
        labels_plasmid=f"{OUT}/08_protein_labels/protein_labels_plasmid.tsv",
    output:
        tsv=f"{OUT}/08_protein_labels/protein_labels.tsv",
    params:
        pfam_version=config["references"]["pfam_version"],
        swissprot_version=config["references"]["swissprot_version"],
        nr_version=config["references"]["nr_version"],
        eggnog_version=config["references"]["eggnog_version"],
        pharokka_db_version=config["references"]["pharokka_db_version"],
    resources:
        mem_mb=32000,
        runtime=240,
    benchmark:
        f"{OUT}/benchmarks/protein_labels.tsv"
    log:
        f"{OUT}/logs/08_protein_labels/protein_labels.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/protein_labels.py"
