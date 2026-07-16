"""Validated internal data models."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field


class FlexresModel(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)


class ResidueKey(FlexresModel):
    pdb_id: str
    model_id: int
    chain_id: str
    residue_number: int
    insertion_code: str = ""
    residue_name: str
    uniprot_id: str
    uniprot_residue_number: int


class AtomCoordinates(FlexresModel):
    n: list[float] | None = None
    ca: list[float] | None = None
    c: list[float] | None = None
    o: list[float] | None = None


class ResidueObservation(FlexresModel):
    key: ResidueKey
    atoms: AtomCoordinates
    all_heavy_atoms: list[list[float]] = Field(default_factory=list)
    side_chain_heavy_atoms: dict[str, list[float]] = Field(default_factory=dict)
    alternate_location_used: str | None = None
    distance_to_reference_ligand_angstrom: float | None = None


class ChainObservation(FlexresModel):
    pdb_id: str
    model_id: int
    chain_id: str
    uniprot_id: str
    residues: dict[int, ResidueObservation]
    experimental_method: str | None = None
    resolution_angstrom: float | None = None
    release_date: str | None = None
    polymer_entity_id: str | None = None
    author_chain_id: str | None = None
    sequence_identity: float | None = None
    reference_sequence_coverage: float | None = None
    comparison_sequence_coverage: float | None = None


class AlignmentResult(FlexresModel):
    rotation: list[list[float]]
    translation: list[float]
    global_ca_rmsd_angstrom: float
    aligned_ca_count: int
    alignment_rmsd_after_outlier_rejection: float | None = None
    aligned_atom_count_after_outlier_rejection: int | None = None
    alignment_cycles: int | None = None
    alignment_rmsd_before_outlier_rejection: float | None = None
    aligned_atom_count_before_outlier_rejection: int | None = None
    raw_alignment_score: float | None = None
    aligned_residue_count: int | None = None
    aligned_residue_pairs: list[dict[str, str | int]] = Field(default_factory=list)
    backend: Literal["pymol", "kabsch", "biopython"] = "kabsch"


class ResidueComparison(FlexresModel):
    reference: ResidueKey
    comparison: ResidueKey
    ca_displacement_angstrom: float | None
    backbone_rmsd_angstrom: float | None
    side_chain_rmsd_angstrom: float | None
    in_reference_pocket: bool
    distance_to_reference_ligand_angstrom: float | None = None
    experimental_method: str | None = None
    resolution_angstrom: float | None = None
    aligned_ca_count: int
    global_ca_rmsd_angstrom: float
    comparison_release_date: str | None = None
    comparison_polymer_entity_id: str | None = None
    reference_sequence_coverage: float | None = None
    comparison_sequence_coverage: float | None = None
    sequence_identity: float | None = None
    alternate_location_used: str | None = None


class SkippedComparison(FlexresModel):
    pdb_id: str
    model: int | None = None
    chain: str | None = None
    stage: str
    reason: str
    details: str = ""


class SkippedResidue(FlexresModel):
    reference_pdb_id: str
    reference_chain: str
    reference_residue_number: int | None = None
    uniprot_residue_number: int | None = None
    comparison_pdb_id: str | None = None
    comparison_chain: str | None = None
    comparison_model: int | None = None
    reason: str
    missing_reference_atoms: str = ""
    missing_comparison_atoms: str = ""


class AnalysisSettings(FlexresModel):
    pdb_id: str
    chain_id: str
    output_dir: Path
    all_residues: bool = False
    ligand_name: str | None = None
    ligand_chain: str | None = None
    ligand_residue: int | None = None
    ligand_insertion_code: str = ""
    pocket_radius: float = 8.0
    include_reference_entry: bool = True
    cache_dir: Path | None = None
    force_redownload: bool = False
    max_workers: int = 4
    alignment_backend: Literal["kabsch", "pymol", "biopython"] = "pymol"
    pymol_executable: str = "pymol"
    alignment_rmsd_cutoff: float = 3.0
    output_formats: set[str] = Field(default_factory=lambda: {"csv", "json"})
    exclude_mutations: bool = True
    resume: bool = False
    log_level: str = "INFO"


class ComparisonTarget(FlexresModel):
    pdb_id: str
    chain_id: str
    entity_id: str | None = None
    author_chain_id: str | None = None
    sequence_identity: float | None = None
    sequence_coverage: float | None = None
    experimental_method: str | None = None
    resolution_angstrom: float | None = None
    release_date: str | None = None


JsonDict = dict[str, Any]
