"""Family identifiers."""


def family_id(resolution, representative_protein_id):
    """Resolution plus the representative's seq_id, so adding an unrelated protein
    elsewhere in the dataset cannot renumber this family."""
    return f"{resolution}:{representative_protein_id}"
