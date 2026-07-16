from __future__ import annotations

import math

import numpy as np
import pytest

from flexres.analysis.alignment import align_chain_biopython, align_chain_kabsch
from flexres.analysis.comparisons import compare_residues
from flexres.analysis.metrics import backbone_rmsd, ca_displacement, side_chain_rmsd
from flexres.models import AtomCoordinates, ChainObservation, ResidueKey, ResidueObservation


def residue(resno: int, ca: list[float], shift_backbone: float = 0.0) -> ResidueObservation:
    key = ResidueKey(
        pdb_id="TEST",
        model_id=1,
        chain_id="A",
        residue_number=resno,
        residue_name="ALA",
        uniprot_id="P1",
        uniprot_residue_number=resno,
    )
    x, y, z = ca
    atoms = AtomCoordinates(
        n=[x - 0.4 + shift_backbone, y, z],
        ca=ca,
        c=[x + 0.4 + shift_backbone, y, z],
        o=[x + 0.8 + shift_backbone, y, z],
    )
    cb = [x, y + 1.0, z]
    return ResidueObservation(key=key, atoms=atoms, all_heavy_atoms=[atoms.n, atoms.ca, atoms.c, atoms.o, cb], side_chain_heavy_atoms={"CB": cb})  # type: ignore[list-item]


def chain(coords: list[list[float]], pdb_id: str = "TEST") -> ChainObservation:
    return ChainObservation(
        pdb_id=pdb_id,
        model_id=1,
        chain_id="A",
        uniprot_id="P1",
        residues={idx + 1: residue(idx + 1, coord) for idx, coord in enumerate(coords)},
    )


def test_ca_displacement_and_backbone_rmsd() -> None:
    assert ca_displacement([0, 0, 0], [3, 4, 0]) == 5.0
    assert math.isclose(
        backbone_rmsd(
            {"N": [0, 0, 0], "CA": [1, 0, 0], "C": [2, 0, 0], "O": [3, 0, 0]},
            {"N": [1, 0, 0], "CA": [2, 0, 0], "C": [3, 0, 0], "O": [4, 0, 0]},
        ),
        1.0,
    )
    assert side_chain_rmsd({"CB": [0, 0, 0]}, {"CB": [0, 3, 4]}) == 5.0


def test_known_rigid_transform_aligns_to_zero() -> None:
    ref = chain([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    theta = np.pi / 2
    rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    comp_coords = [(np.asarray(c) @ rot + np.array([10, -4, 2])).tolist() for c in [[0, 0, 0], [1, 0, 0], [0, 1, 0]]]
    comp = chain(comp_coords, pdb_id="COMP")
    aligned, result = align_chain_kabsch(ref, comp)
    assert result.aligned_ca_count == 3
    assert result.global_ca_rmsd_angstrom < 1e-10
    rows, skipped = compare_residues(ref, aligned, result, {1, 2, 3}, {})
    assert skipped == []
    assert max(row.ca_displacement_angstrom for row in rows) < 1e-10  # type: ignore[arg-type]


def test_structural_variation_after_alignment_is_localized() -> None:
    ref = chain([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])
    comp = chain([[0, 0, 0], [1, 0, 0], [0, 1, 0], [3, 1, 0]], pdb_id="COMP")
    aligned, result = align_chain_kabsch(ref, comp)
    rows, _ = compare_residues(ref, aligned, result, {1, 2, 3, 4}, {})
    by_res = {row.reference.uniprot_residue_number: row for row in rows}
    assert by_res[4].ca_displacement_angstrom > by_res[1].ca_displacement_angstrom


def test_biopython_backend_aligns_and_keeps_outlier_residue_for_reporting() -> None:
    pytest.importorskip("Bio")
    ref = chain([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])
    comp = chain([[0, 0, 0], [1, 0, 0], [0, 1, 0], [4, 1, 0]], pdb_id="COMP")
    aligned, result = align_chain_biopython(ref, comp, outlier_cutoff_angstrom=1.0)
    assert result.backend == "biopython"
    assert result.aligned_ca_count == 3
    assert len(result.aligned_residue_pairs) == 4
    rows, _ = compare_residues(ref, aligned, result, {1, 2, 3, 4}, {})
    by_res = {row.reference.uniprot_residue_number: row for row in rows}
    assert by_res[4].ca_displacement_angstrom > 2.0
