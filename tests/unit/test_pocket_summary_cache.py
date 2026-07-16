from __future__ import annotations

from flexres.analysis.pocket import detect_pocket_residues
from flexres.analysis.summaries import summarize_residue_comparisons
from flexres.cache.manager import CacheManager
from flexres.models import AtomCoordinates, ChainObservation, ResidueKey, ResidueObservation


def obs(resno: int, atom: list[float]) -> ResidueObservation:
    return ResidueObservation(
        key=ResidueKey(
            pdb_id="TEST",
            model_id=1,
            chain_id="A",
            residue_number=resno,
            residue_name="ALA",
            uniprot_id="P1",
            uniprot_residue_number=resno,
        ),
        atoms=AtomCoordinates(ca=atom),
        all_heavy_atoms=[atom],
    )


def test_pocket_detection_uses_heavy_atom_distance() -> None:
    chain = ChainObservation(
        pdb_id="TEST",
        model_id=1,
        chain_id="A",
        uniprot_id="P1",
        residues={1: obs(1, [0, 0, 0]), 2: obs(2, [10, 0, 0])},
    )
    pocket = detect_pocket_residues(chain, [[1, 0, 0]], 8.0)
    assert set(pocket) == {1}


def test_summary_sorting_and_statistics() -> None:
    rows = [
        {"uniprot_id": "P1", "uniprot_residue_number": 1, "reference_residue_name": "ALA", "reference_residue_number": 1, "in_reference_pocket": True, "comparison_pdb_id": "PDB2", "comparison_chain": "A", "comparison_model": 1, "side_chain_rmsd_angstrom": 1.0},
        {"uniprot_id": "P1", "uniprot_residue_number": 1, "reference_residue_name": "ALA", "reference_residue_number": 1, "in_reference_pocket": True, "comparison_pdb_id": "PDB3", "comparison_chain": "B", "comparison_model": 1, "side_chain_rmsd_angstrom": 3.0},
        {"uniprot_id": "P1", "uniprot_residue_number": 2, "reference_residue_name": "GLY", "reference_residue_number": 2, "in_reference_pocket": True, "comparison_pdb_id": "PDB2", "comparison_chain": "A", "comparison_model": 1, "side_chain_rmsd_angstrom": 9.0},
    ]
    summary = summarize_residue_comparisons(rows)
    assert summary[0]["uniprot_residue_number"] == 2
    assert summary[1]["comparison_count"] == 2
    assert "in_reference_pocket" not in summary[1]
    assert summary[1]["PDB2_A_side_chain_rmsd_angstrom"] == 1.0
    assert summary[1]["PDB3_B_side_chain_rmsd_angstrom"] == 3.0
    assert summary[1]["mean_side_chain_rmsd_angstrom"] == 2.0
    assert summary[1]["max_side_chain_rmsd_angstrom"] == 3.0
    assert summary[1]["p90_side_chain_rmsd_angstrom"] > 2.0


def test_cache_atomic_json_roundtrip(tmp_path) -> None:
    cache = CacheManager(tmp_path)
    path = cache.search_path("P1")
    cache.atomic_write_json(path, {"ok": True})
    assert cache.is_valid_file(path)
    assert cache.read_json(path) == {"ok": True}
