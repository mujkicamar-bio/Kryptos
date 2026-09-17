"""Column contracts for every table in the relational store (spec §4.1).

Declaring the columns in one place is what makes the store auditable: a reviewer can read
this file and know exactly what each stage is allowed to produce, and validate_frame turns
a silent schema drift into an immediate failure.

Only the tables Plan 1 writes are declared here. Later plans extend TABLES; they must not
redefine an existing table's columns without a spec revision.
"""

TABLES = {
    # Spec §62.1 - one row per plasmid in the analysis set.
    "plasmids": [
        "plasmid_id", "accession", "source_database", "length", "topology",
        "host", "host_taxonomy", "habitat", "mob_class", "mob_cluster",
        "plasmid_type", "bioproject", "record_class",
    ],
    # Spec §62.2, §8.3 - one row per predicted ORF occurrence.
    "orf_occurrences": [
        "orf_occurrence_id", "plasmid_id", "protein_id",
        "start", "end", "strand", "nucleotide_length", "protein_length",
        "cds_nucleotide_sequence", "partial", "start_type",
        "origin_spanning",
        # QC, spec §9. Flags, never deletions.
        "antifam_hit", "antifam_model", "antifam_score", "antifam_evalue",
        "opposite_strand_overlap_fraction", "opposite_strand_overlap_partner",
        "discovery_excluded_short", "discovery_excluded_partial",
        "discovery_excluded_antifam", "discovery_eligible", "exclusion_reason",
    ],
    # Spec §62.3 - one row per unique protein sequence.
    "proteins": [
        "protein_id", "protein_sequence", "protein_length",
        "occurrence_count", "plasmid_count",
        "has_complete_occurrence", "discovery_eligible", "record_class",
    ],
    # Spec §62.3 addendum - one row per (protein, source, label). LONG, not wide: the
    # label vocabulary is open (Pfam-A 38.2 alone has 30,134 families, and nr product
    # names are unbounded), so it cannot be columns. This table is the substrate from
    # which functional categories - replication, mobilisation, conjugation - are derived.
    # The hand-curated family list it replaces named 15 of the 42 Pfam families whose
    # description mentions conjugation, and assigned every role with no source.
    "protein_labels": [
        "protein_id", "source", "tier", "kind", "label", "accession",
        "evidence_evalue", "evidence_coverage", "database", "database_version",
    ],
    # Spec §58 - positive and negative controls travel with the real data and are excluded
    # from biological occurrence-integrity checks.
    "controls": [
        "protein_id", "record_class", "control_category",
        "source_protein_id", "expected_outcome",
    ],
    # Spec §6 - one row per production run.
    "run_manifest": [
        "pipeline_version", "git_commit", "input_dataset_hash",
        "plasmid_count", "orf_count", "unique_protein_count",
        "config_hash", "database_versions", "software_versions",
        "shard_counts", "started_at", "finished_at",
    ],
}


def columns(table):
    """Ordered column names for one table."""
    if table not in TABLES:
        raise KeyError(f"no schema declared for table {table!r}")
    return list(TABLES[table])


def validate_frame(table, frame):
    """Raise if an Arrow table does not match its declared contract.

    Both directions are checked. A missing column breaks consumers downstream; an
    undeclared one means a stage is emitting something no schema documents, which is how
    columns with two different meanings end up sharing a name.
    """
    declared = set(columns(table))
    present = set(frame.column_names)
    missing = declared - present
    undeclared = present - declared
    if missing:
        raise ValueError(f"{table}: missing declared columns {sorted(missing)}")
    if undeclared:
        raise ValueError(f"{table}: undeclared columns {sorted(undeclared)}")
