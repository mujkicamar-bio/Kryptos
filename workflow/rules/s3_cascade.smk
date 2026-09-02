rule tier_search:
    """Search one tier, then hand the next tier only what stayed unnamed."""
    input:
        faa=tier_query,
    output:
        hits=f"{OUT}/s3/{{tier}}/hits.tsv",
        unresolved=f"{OUT}/s3/{{tier}}/unresolved.faa",
    params:
        spec=lambda wc: TIER_BY_ID[wc.tier],
    threads: 16
    conda:
        "../envs/search.yaml"
    log:
        f"{OUT}/logs/s3/{{tier}}.log",
    script:
        "../scripts/tier_search.py"


rule cascade_resolve:
    """Assign functional_class per protein from all tier hits (design 6.1-6.3)."""
    input:
        hits=expand(f"{OUT}/s3/{{tier}}/hits.tsv", tier=TIER_IDS),
        faa=f"{OUT}/s2/unique_proteins.faa",
    output:
        f"{OUT}/s3/protein_annotation.tsv",
    log:
        f"{OUT}/logs/s3/resolve.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/cascade_resolve.py"


rule annotate_plasmids:
    """S4: the primary deliverable - every ORF with its annotation and provenance."""
    input:
        prot=f"{OUT}/s3/protein_annotation.tsv",
        index=f"{OUT}/s1/orf_index.tsv",
        map=f"{OUT}/s2/protein_map.tsv",
    output:
        f"{OUT}/s4/plasmid_annotation.tsv",
    log:
        f"{OUT}/logs/s4/annotate.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/annotate_plasmids.py"
