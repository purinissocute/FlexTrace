"""Residue-level metric calculations."""

from __future__ import annotations

import numpy as np

from flexres.structures.atoms import BACKBONE_ATOMS


def ca_displacement(reference_ca: list[float], comparison_ca: list[float]) -> float:
    """Euclidean C-alpha displacement in angstroms."""
    return float(np.linalg.norm(np.asarray(comparison_ca, dtype=float) - np.asarray(reference_ca, dtype=float)))


def backbone_rmsd(reference_atoms: dict[str, list[float]], comparison_atoms: dict[str, list[float]]) -> float:
    """Backbone N-CA-C-O RMSD after global chain alignment."""
    diffs = []
    for atom in BACKBONE_ATOMS:
        ref = np.asarray(reference_atoms[atom], dtype=float)
        comp = np.asarray(comparison_atoms[atom], dtype=float)
        diffs.append(float(np.sum((comp - ref) ** 2)))
    return float(np.sqrt(sum(diffs) / len(BACKBONE_ATOMS)))


def side_chain_rmsd(reference_atoms: dict[str, list[float]], comparison_atoms: dict[str, list[float]]) -> float:
    """Side-chain heavy-atom RMSD after global chain alignment."""
    common_atoms = sorted(set(reference_atoms) & set(comparison_atoms))
    if not common_atoms:
        raise ValueError("no common side-chain heavy atoms")
    diffs = []
    for atom in common_atoms:
        ref = np.asarray(reference_atoms[atom], dtype=float)
        comp = np.asarray(comparison_atoms[atom], dtype=float)
        diffs.append(float(np.sum((comp - ref) ** 2)))
    return float(np.sqrt(sum(diffs) / len(common_atoms)))


def missing_backbone_atoms(atom_values: object) -> list[str]:
    return [name for name in BACKBONE_ATOMS if getattr(atom_values, name.lower()) is None]
