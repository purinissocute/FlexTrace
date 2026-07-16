"""Reference-chain validation and mapping."""

from __future__ import annotations

from pathlib import Path

from flexres.api.sifts import SIFTSClient
from flexres.exceptions import ChainNotFoundError, ReferenceStructureError
from flexres.mapping.validation import require_single_uniprot
from flexres.structures.chains import chain_is_protein, get_chain, get_model
from flexres.structures.parser import read_structure


def validate_reference_chain(path: Path, pdb_id: str, chain_id: str) -> None:
    structure = read_structure(path)
    model = get_model(structure, 1)
    try:
        chain = get_chain(model, chain_id)
    except ChainNotFoundError:
        raise ChainNotFoundError(f"chain {chain_id} was not found in reference PDB {pdb_id}") from None
    if not chain_is_protein(chain):
        raise ReferenceStructureError(f"chain {chain_id} in {pdb_id} is not a protein chain")


def reference_uniprot(sifts: SIFTSClient, pdb_id: str, chain_id: str, force: bool = False) -> str:
    accessions = sifts.chain_uniprot_accessions(pdb_id, chain_id, force=force)
    return require_single_uniprot(accessions, pdb_id, chain_id)
