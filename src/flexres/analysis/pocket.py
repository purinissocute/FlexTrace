"""Reference ligand pocket detection."""

from __future__ import annotations

from flexres.models import ChainObservation
from flexres.structures.atoms import min_distance


def detect_pocket_residues(
    reference: ChainObservation,
    ligand_heavy_atoms: list[list[float]],
    radius_angstrom: float,
) -> dict[int, float]:
    """Return UniProt residue numbers within heavy-atom distance of ligand."""
    pocket: dict[int, float] = {}
    for uniprot_resno, residue in reference.residues.items():
        if not residue.all_heavy_atoms:
            continue
        distance = min_distance(residue.all_heavy_atoms, ligand_heavy_atoms)
        if distance <= radius_angstrom:
            pocket[uniprot_resno] = distance
    return dict(sorted(pocket.items()))
