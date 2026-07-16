from __future__ import annotations

from flexres.analysis.comparisons import compare_residues
from flexres.models import AlignmentResult, AnalysisSettings, AtomCoordinates, ChainObservation, ResidueKey, ResidueObservation
from flexres.output.schemas import flatten_comparison


def make_obs(resno: int, atoms: AtomCoordinates, side_chain: dict[str, list[float]] | None = None) -> ResidueObservation:
    return ResidueObservation(
        key=ResidueKey(
            pdb_id="PDB1",
            model_id=1,
            chain_id="A",
            residue_number=resno,
            residue_name="ALA",
            uniprot_id="P1",
            uniprot_residue_number=resno,
        ),
        atoms=atoms,
        all_heavy_atoms=[],
        side_chain_heavy_atoms=side_chain or {"CB": [1, 1, 0]},
    )


def test_missing_backbone_atom_is_skipped_not_zero() -> None:
    ref = ChainObservation(
        pdb_id="PDB1",
        model_id=1,
        chain_id="A",
        uniprot_id="P1",
        residues={1: make_obs(1, AtomCoordinates(n=[0, 0, 0], ca=[1, 0, 0], c=[2, 0, 0], o=[3, 0, 0]))},
    )
    comp = ChainObservation(
        pdb_id="PDB2",
        model_id=1,
        chain_id="A",
        uniprot_id="P1",
        residues={1: make_obs(1, AtomCoordinates(n=[0, 0, 0], ca=[1, 0, 0], c=[2, 0, 0]))},
    )
    alignment = AlignmentResult(rotation=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], translation=[0, 0, 0], global_ca_rmsd_angstrom=0, aligned_ca_count=1)
    rows, skipped = compare_residues(ref, comp, alignment, {1}, {})
    assert rows == []
    assert skipped[0].reason == "missing_backbone_atoms"
    assert skipped[0].missing_comparison_atoms == "O"


def test_flatten_comparison_has_required_columns() -> None:
    ref = ChainObservation(
        pdb_id="PDB1",
        model_id=1,
        chain_id="A",
        uniprot_id="P1",
        residues={1: make_obs(1, AtomCoordinates(n=[0, 0, 0], ca=[1, 0, 0], c=[2, 0, 0], o=[3, 0, 0]))},
    )
    comp = ChainObservation(
        pdb_id="PDB2",
        model_id=1,
        chain_id="A",
        uniprot_id="P1",
        residues={1: make_obs(1, AtomCoordinates(n=[0, 0, 0], ca=[2, 0, 0], c=[2, 0, 0], o=[3, 0, 0]))},
    )
    alignment = AlignmentResult(rotation=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], translation=[0, 0, 0], global_ca_rmsd_angstrom=0, aligned_ca_count=1)
    rows, _ = compare_residues(ref, comp, alignment, {1}, {1: 4.0})
    flat = flatten_comparison(rows[0])
    assert flat["reference_pdb_id"] == "PDB1"
    assert flat["comparison_pdb_id"] == "PDB1"
    assert flat["side_chain_rmsd_angstrom"] == 0.0
    assert flat["in_reference_pocket"] is True


def test_mutated_comparison_residue_is_skipped_not_whole_chain() -> None:
    ref = ChainObservation(
        pdb_id="PDB1",
        model_id=1,
        chain_id="A",
        uniprot_id="P1",
        residues={
            1: make_obs(1, AtomCoordinates(n=[0, 0, 0], ca=[1, 0, 0], c=[2, 0, 0], o=[3, 0, 0])),
            2: make_obs(2, AtomCoordinates(n=[0, 1, 0], ca=[1, 1, 0], c=[2, 1, 0], o=[3, 1, 0])),
        },
    )
    comp = ChainObservation(
        pdb_id="PDB2",
        model_id=1,
        chain_id="A",
        uniprot_id="P1",
        residues={
            1: make_obs(1, AtomCoordinates(n=[0, 0, 0], ca=[1, 0, 0], c=[2, 0, 0], o=[3, 0, 0])),
            2: make_obs(2, AtomCoordinates(n=[0, 1, 0], ca=[1, 1, 0], c=[2, 1, 0], o=[3, 1, 0])),
        },
    )
    alignment = AlignmentResult(rotation=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], translation=[0, 0, 0], global_ca_rmsd_angstrom=0, aligned_ca_count=2)
    rows, skipped = compare_residues(ref, comp, alignment, {1, 2}, {}, {2: "A2V"})
    assert [row.reference.uniprot_residue_number for row in rows] == [1]
    assert skipped[0].reason == "comparison_sequence_mismatch_mutation"
    assert skipped[0].uniprot_residue_number == 2
    assert skipped[0].missing_comparison_atoms == "A2V"


def test_default_alignment_rmsd_cutoff_is_three_angstroms(tmp_path) -> None:
    settings = AnalysisSettings(pdb_id="TEST", chain_id="A", output_dir=tmp_path, all_residues=True)
    assert settings.alignment_backend == "pymol"
    assert settings.alignment_rmsd_cutoff == 3.0
