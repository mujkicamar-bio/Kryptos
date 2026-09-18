# =====================================================================================
# S3: the annotation cascade
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
        f"{OUT}/s3/preflight.tsv",
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
    log:
        f"{OUT}/logs/s3/preflight.log",
    script:
        "../scripts/preflight.py"


rule sweep_cohort:
    """Proteins that bypass narrowing entirely, so the threshold's cost is measurable.

    Roughly 2% of the compute buys the ability to state what narrow_at cost, with a
    confidence interval. Everything else in the run has no counterfactual.
    """
    input:
        f"{OUT}/s2/cascade_input.faa",
    output:
        f"{OUT}/s3/sweep_cohort.txt",
    params:
        fraction=cascade["sweep_cohort_fraction"],
        seed=config["seed"],
        # Checked here, against the real protein count, before any tier runs.
        hmmer_z=cascade["hmmer_z"],
    conda:
        "../envs/plasmidann.yaml"
    resources:
        mem_mb=8000,
        runtime=30,
    log:
        f"{OUT}/logs/s3/sweep_cohort.log",
    script:
        "../scripts/sweep_cohort.py"


rule shard_cascade_input:
    """S3: split the query set into independently searchable, resumable units.

    Without this the deepest tier is a single job of several days over 3.5M queries against
    a 375 GB database, and any failure in it loses all of that work. Every other long stage
    was already sharded; the cascade was the one that was not.
    """
    input:
        f"{OUT}/s2/cascade_input.faa",
    output:
        expand(f"{OUT}/s3/input/{{cshard}}.faa", cshard=CASCADE_SHARDS),
    resources:
        mem_mb=8000,
        runtime=60,
    log:
        f"{OUT}/logs/s3/shard_cascade_input.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/shard_cascade_input.py"


rule tier_search:
    """Search one tier of one shard, then hand the next tier only what stayed unexplained.

    The unit of work is (tier, shard). A failure costs one shard of one tier rather than
    the tier, and the shards of a tier run concurrently across the allocation.
    """
    input:
        faa=tier_query,
        spans=tier_spans,
        sweep=f"{OUT}/s3/sweep_cohort.txt",
        preflight=f"{OUT}/s3/preflight.tsv",
    output:
        hits=f"{OUT}/s3/{{tier}}/{{cshard}}/hits.tsv",
        unresolved=f"{OUT}/s3/{{tier}}/{{cshard}}/unresolved.faa",
        spans=f"{OUT}/s3/{{tier}}/{{cshard}}/spans.tsv",
    params:
        spec=lambda wc: TIER_BY_ID[wc.tier],
        # narrow_at, NOT min_explained. See the module docstring in tier_search.py.
        narrow_at=cascade["narrow_at"],
        hmmer_z=cascade["hmmer_z"],
        max_target_seqs=cascade["max_target_seqs"],
    threads: 4
    conda:
        "../envs/plasmidann.yaml"
    resources:
        # PER TIER, not one number for all 256 jobs. A T1 shard reads a 2.2 GB pressed HMM
        # library; a T4 shard streams 350 GB of nr and DIAMOND holds a large block of it in
        # memory. Giving both the same 16 GB is how a tier gets OOM-killed on day three,
        # after every shard has already cleared T1-T3.
        #
        # The T4 numbers are a deliberate over-allocation, not a measurement:
        # workflow/bench_nr.sbatch exists to measure the real figure and needs a project
        # allocation. Over-allocating costs concurrency; under-allocating costs the run.
        mem_mb=lambda wc: 64000 if TIER_BY_ID[wc.tier]["db"].endswith("nr.dmnd") else 16000,
        runtime=lambda wc: 2880 if TIER_BY_ID[wc.tier]["db"].endswith("nr.dmnd") else 1440,
    log:
        f"{OUT}/logs/s3/{{tier}}/{{cshard}}.log",
    script:
        "../scripts/tier_search.py"


rule cascade_resolve:
    """Assign functional_class per protein from every tier's hits.

    min_explained is applied here, post hoc, on a table where every protein in the
    interesting band has been seen by every tier - which is what makes it sweepable.
    """
    input:
        hits=expand(f"{OUT}/s3/{{tier}}/{{cshard}}/hits.tsv",
                    tier=TIER_IDS, cshard=CASCADE_SHARDS),
        # The last tier's spans carry the cumulative explained fraction for every protein
        # the cascade ever saw, because each tier writes forward everything it inherited.
        # One file per shard, and a protein appears in exactly one of them.
        spans=expand(f"{OUT}/s3/{TIER_IDS[-1]}/{{cshard}}/spans.tsv",
                     cshard=CASCADE_SHARDS),
        faa=f"{OUT}/s2/cascade_input.faa",
    output:
        f"{OUT}/s3/protein_annotation.tsv",
    params:
        thresholds=cascade,
        tier_order=TIER_IDS,
    resources:
        mem_mb=32000,
        runtime=480,
    log:
        f"{OUT}/logs/s3/resolve.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/cascade_resolve.py"


rule annotate_plasmids:
    """S4: the primary deliverable - every ORF with its annotation and provenance."""
    input:
        prot=f"{OUT}/s3/protein_annotation.tsv",
        index=f"{OUT}/s1/orf_index.tsv",
        map=f"{OUT}/s2/protein_map.tsv",
        artefact=f"{OUT}/s2b/artefact_flags.tsv",
    output:
        f"{OUT}/s4/plasmid_annotation.tsv",
    resources:
        mem_mb=24000,
        runtime=240,
    log:
        f"{OUT}/logs/s4/annotate.log",
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
        prot=f"{OUT}/s3/protein_annotation.tsv",
        faa=f"{OUT}/s2/unique_proteins.faa",
    output:
        f"{OUT}/s4b/orthology.tsv",
    params:
        orthology=targets["orthology"],
    threads: 16
    resources:
        mem_mb=32000,
        runtime=2880,
    log:
        f"{OUT}/logs/s4b/orthology.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/orthology.py"


rule feature_files:
    """S4: GFF3 and GenBank alongside the TSV.

    160,375 ORFs were reconstructed across the origin of a circular plasmid and carry
    start > end. GFF3 forbids that and GenBank has dedicated syntax for it, so the same
    gene is written two different ways - see plasmidann.features.

    Takes the SHARDS, not the corpus FASTA. Reading the corpus made the scope of these two
    files a property of a config path rather than of the input the run was given: on the
    100-plasmid test configuration it wrote 208,245 GenBank records.
    """
    input:
        annotation=f"{OUT}/s4/plasmid_annotation.tsv",
        shards=[SHARD_PATHS[s] for s in SHARDS],
        master=config["input"]["master_table"],
    output:
        gff3=f"{OUT}/s4/plasmid_annotation.gff3",
        genbank=f"{OUT}/s4/plasmid_annotation.gbk",
    resources:
        mem_mb=32000,
        runtime=480,
    log:
        f"{OUT}/logs/s4/feature_files.log",
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
        hits=expand(f"{OUT}/s3/{{tier}}/{{cshard}}/hits.tsv",
                    tier=TIER_IDS, cshard=CASCADE_SHARDS),
        orthology=f"{OUT}/s4b/orthology.tsv",
        pfam_dat=config["references"]["pfam_dat"],
    output:
        tsv=f"{OUT}/s4c/protein_labels.tsv",
    params:
        pfam_version=config["references"]["pfam_version"],
        swissprot_version=config["references"]["swissprot_version"],
        nr_version=config["references"]["nr_version"],
        eggnog_version=config["references"]["eggnog_version"],
    resources:
        mem_mb=32000,
        runtime=240,
    log:
        f"{OUT}/logs/s4c/protein_labels.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/protein_labels.py"
