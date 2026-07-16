"""Residue mapping transformations."""

from __future__ import annotations


def invert_uniprot_mapping(mapping: dict[tuple[int, str], tuple[str, int]], accession: str) -> dict[int, tuple[int, str]]:
    """Return UniProt residue number to PDB residue ID mapping for one accession."""
    result: dict[int, tuple[int, str]] = {}
    for pdb_key, (mapped_accession, uniprot_resno) in mapping.items():
        if mapped_accession == accession:
            result[uniprot_resno] = pdb_key
    return result
