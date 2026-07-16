"""Output column definitions and flatteners."""

from __future__ import annotations

from flexres.models import ResidueComparison, SkippedComparison, SkippedResidue

RESIDUE_COMPARISON_COLUMNS = [
    "reference_pdb_id", "reference_model", "reference_chain", "reference_residue_name",
    "reference_residue_number", "reference_insertion_code", "uniprot_id",
    "uniprot_residue_number", "comparison_pdb_id", "comparison_model", "comparison_chain",
    "comparison_residue_name", "comparison_residue_number", "comparison_insertion_code",
    "experimental_method", "resolution_angstrom", "aligned_ca_count", "global_ca_rmsd_angstrom",
    "ca_displacement_angstrom", "backbone_rmsd_angstrom", "side_chain_rmsd_angstrom", "in_reference_pocket",
    "distance_to_reference_ligand_angstrom", "comparison_release_date", "comparison_polymer_entity_id",
    "reference_sequence_coverage", "comparison_sequence_coverage", "sequence_identity", "alternate_location_used",
]

SKIPPED_COMPARISON_COLUMNS = ["pdb_id", "model", "chain", "stage", "reason", "details"]
SKIPPED_RESIDUE_COLUMNS = [
    "reference_pdb_id", "reference_chain", "reference_residue_number", "uniprot_residue_number",
    "comparison_pdb_id", "comparison_chain", "comparison_model", "reason", "missing_reference_atoms",
    "missing_comparison_atoms",
]


def flatten_comparison(row: ResidueComparison) -> dict:
    return {
        "reference_pdb_id": row.reference.pdb_id,
        "reference_model": row.reference.model_id,
        "reference_chain": row.reference.chain_id,
        "reference_residue_name": row.reference.residue_name,
        "reference_residue_number": row.reference.residue_number,
        "reference_insertion_code": row.reference.insertion_code,
        "uniprot_id": row.reference.uniprot_id,
        "uniprot_residue_number": row.reference.uniprot_residue_number,
        "comparison_pdb_id": row.comparison.pdb_id,
        "comparison_model": row.comparison.model_id,
        "comparison_chain": row.comparison.chain_id,
        "comparison_residue_name": row.comparison.residue_name,
        "comparison_residue_number": row.comparison.residue_number,
        "comparison_insertion_code": row.comparison.insertion_code,
        "experimental_method": row.experimental_method,
        "resolution_angstrom": row.resolution_angstrom,
        "aligned_ca_count": row.aligned_ca_count,
        "global_ca_rmsd_angstrom": row.global_ca_rmsd_angstrom,
        "ca_displacement_angstrom": row.ca_displacement_angstrom,
        "backbone_rmsd_angstrom": row.backbone_rmsd_angstrom,
        "side_chain_rmsd_angstrom": row.side_chain_rmsd_angstrom,
        "in_reference_pocket": row.in_reference_pocket,
        "distance_to_reference_ligand_angstrom": row.distance_to_reference_ligand_angstrom,
        "comparison_release_date": row.comparison_release_date,
        "comparison_polymer_entity_id": row.comparison_polymer_entity_id,
        "reference_sequence_coverage": row.reference_sequence_coverage,
        "comparison_sequence_coverage": row.comparison_sequence_coverage,
        "sequence_identity": row.sequence_identity,
        "alternate_location_used": row.alternate_location_used,
    }


def flatten_skipped_comparison(row: SkippedComparison) -> dict:
    return row.model_dump()


def flatten_skipped_residue(row: SkippedResidue) -> dict:
    return row.model_dump()
