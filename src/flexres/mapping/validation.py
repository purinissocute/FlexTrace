"""Mapping validation."""

from __future__ import annotations

from flexres.exceptions import AmbiguousUniProtMappingError, UniProtMappingError


def require_single_uniprot(accessions: set[str], pdb_id: str, chain_id: str) -> str:
    if not accessions:
        raise UniProtMappingError(f"no UniProt mapping is available for {pdb_id} chain {chain_id}")
    if len(accessions) > 1:
        values = ", ".join(sorted(accessions))
        raise AmbiguousUniProtMappingError(f"{pdb_id} chain {chain_id} maps to multiple UniProt accessions: {values}")
    return next(iter(accessions))
