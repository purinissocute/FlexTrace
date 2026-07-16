"""Chain selection and validation."""

from __future__ import annotations

from flexres.exceptions import ChainNotFoundError
from flexres.structures.atoms import is_amino_acid


def get_model(structure: object, model_id: int = 1) -> object:
    models = list(structure)
    if model_id < 1 or model_id > len(models):
        raise ChainNotFoundError(f"model {model_id} was not found")
    return models[model_id - 1]


def get_chain(model: object, chain_id: str) -> object:
    for chain in model:
        if getattr(chain, "name", "") == chain_id:
            return chain
    raise ChainNotFoundError(f"chain {chain_id} was not found")


def chain_is_protein(chain: object) -> bool:
    return any(is_amino_acid(residue) for residue in chain)


def model_count(structure: object) -> int:
    return len(list(structure))
