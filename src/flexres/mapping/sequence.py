"""Sequence identity checks against UniProt residue positions."""

from __future__ import annotations

from flexres.models import ChainObservation
from flexres.structures.atoms import AA3_TO_1


def canonical_mismatches(chain: ChainObservation, canonical_sequence: str, limit: int = 20) -> list[str]:
    """Return mapped residues that do not match the canonical UniProt sequence."""
    return [mismatch for _, mismatch in canonical_mismatch_positions(chain, canonical_sequence, limit)]


def canonical_mismatch_positions(chain: ChainObservation, canonical_sequence: str, limit: int = 20) -> list[tuple[int, str]]:
    """Return UniProt residue numbers and labels that differ from the canonical sequence."""
    mismatches: list[str] = []
    positioned: list[tuple[int, str]] = []
    for uniprot_resno, residue in sorted(chain.residues.items()):
        if uniprot_resno < 1 or uniprot_resno > len(canonical_sequence):
            continue
        observed = AA3_TO_1.get(residue.key.residue_name.upper())
        expected = canonical_sequence[uniprot_resno - 1]
        if observed is None or observed == expected:
            continue
        mismatch = f"{expected}{uniprot_resno}{observed}"
        mismatches.append(mismatch)
        positioned.append((uniprot_resno, mismatch))
        if len(positioned) >= limit:
            break
    return positioned
