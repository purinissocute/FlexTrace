from __future__ import annotations

from flexres.mapping.sequence import canonical_mismatch_positions, canonical_mismatches
from flexres.models import AtomCoordinates, ChainObservation, ResidueKey, ResidueObservation


def obs(resno: int, residue_name: str) -> ResidueObservation:
    return ResidueObservation(
        key=ResidueKey(
            pdb_id="TEST",
            model_id=1,
            chain_id="A",
            residue_number=resno,
            residue_name=residue_name,
            uniprot_id="P1",
            uniprot_residue_number=resno,
        ),
        atoms=AtomCoordinates(ca=[0, 0, 0]),
    )


def test_canonical_mismatches_detects_point_mutation() -> None:
    chain = ChainObservation(
        pdb_id="TEST",
        model_id=1,
        chain_id="A",
        uniprot_id="P1",
        residues={1: obs(1, "ALA"), 2: obs(2, "VAL")},
    )
    assert canonical_mismatches(chain, "AI") == ["I2V"]
    assert canonical_mismatch_positions(chain, "AI") == [(2, "I2V")]


def test_canonical_mismatches_treats_mse_as_methionine() -> None:
    chain = ChainObservation(
        pdb_id="TEST",
        model_id=1,
        chain_id="A",
        uniprot_id="P1",
        residues={1: obs(1, "MSE")},
    )
    assert canonical_mismatches(chain, "M") == []
