# =====================================================================================
# The annotation cascade
#
# Each tier searches whatever the previous tier could not explain, and hands on the
# residue. Tier order is authority order.
#
# The one structural change from v1: narrowing uses `narrow_at` (permissive), while
# `min_explained` is applied post hoc at cascade_resolve. Using one number for both made
# the threshold unsweepable, because a protein withheld at T2 has no T4 result.
# =====================================================================================

rule preflight:
    """Fail in seconds if a tool or database is missing, rather than after 45 hours.

    A dependency of every tier and of the artefact screen, so no search can start until
    every configured tool and database has been confirmed present. v1 lost a 45-hour job
    to a missing DIAMOND binary, twice.
    """
    output:
        f"{OUT}/05_annotation_cascade/preflight.tsv",
    params:
        tiers=cascade["tiers"],
        artefact=cascade["artefact_screen"],
        structure=targets["structure"],
        orthology=targets["orthology"],
        foldseek_db=config.get("foldseek_db", "data/refs/foldseek/pdb"),
        prostt5=config.get("prostt5_model", "data/refs/foldseek/prostt5"),
    conda:
        "../envs/plasmidann.yaml"
    resources:
        mem_mb=2000,
        runtime=10,
    benchmark:
        f"{OUT}/benchmarks/preflight.tsv"
    log:
        f"{OUT}/logs/05_annotation_cascade/preflight.log",
    script:
        "../scripts/preflight.py"


rule check_hmmer_z:
    """S2z: refuse a declared -Z that no longer describes the protein set, as soon as
    dereplication has produced it - before the artefact screen or any tier searches."""
    input:
        unique=f"{OUT}/03_dereplication/unique_proteins.faa",
    output:
        f"{OUT}/03_dereplication/hmmer_z_checked.tsv",
    params:
        hmmer_z=cascade["hmmer_z"],
        n_controls=(config["controls"]["positive"]["n"]
                    + config["controls"]["negative"]["n_shuffled"]
                    + config["controls"]["negative"]["n_reverse_complement"]),
    conda:
        "../envs/plasmidann.yaml"
    resources:
        mem_mb=2000,
        runtime=30,
    log:
        f"{OUT}/logs/03_dereplication/check_hmmer_z.log",
    script:
        "../scripts/check_hmmer_z.py"


rule sweep_cohort:
    """Proteins that bypass narrowing entirely, so the threshold's cost is measurable.

    Roughly 2% of the compute buys the ability to state what narrow_at cost, with a
    confidence interval. Everything else in the run has no counterfactual.
    """
    input:
        faa=f"{OUT}/03_dereplication/cascade_input.faa",
        # Every tier depends on this rule, so none searches with an unconfirmed -Z.
        hmmer_z=f"{OUT}/03_dereplication/hmmer_z_checked.tsv",
    output:
        f"{OUT}/05_annotation_cascade/sweep_cohort.txt",
    params:
        fraction=cascade["sweep_cohort_fraction"],
        seed=config["seed"],
    conda:
        "../envs/plasmidann.yaml"
    resources:
        mem_mb=8000,
        runtime=30,
    benchmark:
        f"{OUT}/benchmarks/sweep_cohort.tsv"
    log:
        f"{OUT}/logs/05_annotation_cascade/sweep_cohort.log",
    script:
        "../scripts/sweep_cohort.py"


rule tier_search:
    """Search one tier, then hand the next tier only what stayed unexplained.

    One job per tier, with every core the run has. A search against a streamed database
    - DIAMOND reads the whole of nr per invocation - has a fixed cost that sharding the
    query set multiplies by the shard count, which is why the cascade is not sharded.
    """
    input:
        faa=tier_query,
        spans=tier_spans,
        # Earlier tiers' hits, where this tier skips what they named (skip_if_named_by).
        named=tier_named,
        sweep=f"{OUT}/05_annotation_cascade/sweep_cohort.txt",
        preflight=f"{OUT}/05_annotation_cascade/preflight.tsv",
    output:
        hits=f"{OUT}/05_annotation_cascade/{{tier}}/hits.tsv",
        unresolved=f"{OUT}/05_annotation_cascade/{{tier}}/unresolved.faa",
        spans=f"{OUT}/05_annotation_cascade/{{tier}}/spans.tsv",
    params:
        spec=lambda wc: TIER_BY_ID[wc.tier],
        # narrow_at, NOT min_explained. See the module docstring in tier_search.py.
        narrow_at=cascade["narrow_at"],
        hmmer_z=cascade["hmmer_z"],
        max_target_seqs=cascade["max_target_seqs"],
    threads: workflow.cores
    conda:
        "../envs/plasmidann.yaml"
    resources:
        # PER TIER, and MEASURED rather than declared, on full nr (T5 is now ClusteredNR,
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
        # nr on full nr: ~45 min fixed per pass and ~1.7 s per query on 96 threads
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
    """Assign functional_class per protein from every tier's hits.

    min_explained is applied here, post hoc, on a table where every protein in the
    interesting band has been seen by every tier - which is what makes it sweepable.
    """
    input:
        hits=expand(f"{OUT}/05_annotation_cascade/{{tier}}/hits.tsv", tier=TIER_IDS),
        # The last tier's spans carry the cumulative explained fraction for every protein
        # the cascade ever saw, because each tier writes forward everything it inherited.
        spans=f"{OUT}/05_annotation_cascade/{TIER_IDS[-1]}/spans.tsv",
        faa=f"{OUT}/03_dereplication/cascade_input.faa",
        # The proteins that skipped the cascade, so the table covers every protein.
        ps=f"{OUT}/03_dereplication/plasmidscope_proteins.tsv",
        # Search-cluster members, and the proteins no selected family holds (S2s).
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
    """S4b: COG and KEGG terms for the proteins the cascade named.

    Not a dark-hunting tier. S8 asks what a dark ORF's neighbours do, and the cascade
    answers in free text, which cannot be aggregated into pathways. Annotating the known
    fraction is what makes the unknown fraction interpretable.
    """
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
    """S4: GFF3 and GenBank alongside the TSV.

    160,375 ORFs were reconstructed across the origin of a circular plasmid and carry
    start > end. GFF3 forbids that and GenBank has dedicated syntax for it, so the same
    gene is written two different ways - see plasmidann.features.

    Takes the analysis-set FASTA, not the corpus. Reading the corpus made the scope of
    these two files a property of a config path rather than of the analysis set: on the
    100-plasmid test configuration it wrote 208,245 GenBank records.
    """
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


rule protein_labels:
    """S4c: every functional label every tool produced, in one long table.

    The substrate for grouping proteins into replication, mobilisation and conjugation.
    Long rather than wide because the vocabulary is open: Pfam-A 38.2 alone has 30,134
    families, of which 67 mention replication in their description and 42 mention
    conjugation - against the 16 and 15 named by the curated list this replaces.

    Nothing here assigns a biological role. The grouping is derived from the labels
    observed, which is what makes it describable in a methods section.
    """
    input:
        hits=expand(f"{OUT}/05_annotation_cascade/{{tier}}/hits.tsv", tier=TIER_IDS),
        orthology=f"{OUT}/07_orthology/orthology.tsv",
        pfam_dat=config["references"]["pfam_dat"],
        # Search-cluster members take their representative's labels.
        selection=f"{OUT}/03_dereplication/selection.tsv",
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
