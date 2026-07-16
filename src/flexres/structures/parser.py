"""Parse mmCIF files into analysis-ready observations."""

from __future__ import annotations

from pathlib import Path

import gemmi

from flexres.exceptions import LigandNotFoundError, StructureParseError
from flexres.models import AtomCoordinates, ChainObservation, ResidueKey, ResidueObservation
from flexres.structures.atoms import BACKBONE_ATOMS, choose_altloc_atoms, clean_gemmi_text, coord_list, is_amino_acid, is_hydrogen
from flexres.structures.chains import chain_is_protein, get_chain, get_model


def read_structure(path: Path) -> gemmi.Structure:
    try:
        return gemmi.read_structure(str(path))
    except Exception as exc:  # noqa: BLE001
        raise StructureParseError(f"failed to parse {path}: {exc}") from exc


def residue_number_and_icode(residue: object) -> tuple[int, str]:
    seqid = getattr(residue, "seqid")
    return int(seqid.num), clean_gemmi_text(getattr(seqid, "icode", ""))


def residue_label_number(residue: object) -> int | None:
    label_seq = getattr(residue, "label_seq", None)
    return int(label_seq) if label_seq is not None else None


def parse_chain_observation(
    path: Path,
    pdb_id: str,
    chain_id: str,
    uniprot_id: str,
    residue_map: dict[tuple[int, str], tuple[str, int]],
    model_id: int = 1,
    experimental_method: str | None = None,
    resolution_angstrom: float | None = None,
) -> ChainObservation:
    structure = read_structure(path)
    model = get_model(structure, model_id)
    chain = get_chain(model, chain_id)
    if not chain_is_protein(chain):
        raise StructureParseError(f"chain {chain_id} in {pdb_id} is not a protein chain")
    residues: dict[int, ResidueObservation] = {}
    for residue in chain:
        if not is_amino_acid(residue):
            continue
        pdb_resno, icode = residue_number_and_icode(residue)
        label_resno = residue_label_number(residue)
        mapping = (
            residue_map.get((pdb_resno, icode))
            or residue_map.get((pdb_resno, ""))
            or (residue_map.get((label_resno, "")) if label_resno is not None else None)
        )
        if mapping is None:
            continue
        mapped_accession, uniprot_resno = mapping
        if mapped_accession != uniprot_id:
            continue
        chosen, altloc_used = choose_altloc_atoms(list(residue))
        atom_coords = {name.lower(): coord_list(chosen[name]) if name in chosen else None for name in BACKBONE_ATOMS}
        heavy = [coord_list(atom) for atom in chosen.values() if not is_hydrogen(atom)]
        side_chain = {
            name: coord_list(atom)
            for name, atom in chosen.items()
            if name not in BACKBONE_ATOMS and not is_hydrogen(atom)
        }
        key = ResidueKey(
            pdb_id=pdb_id.upper(),
            model_id=model_id,
            chain_id=chain_id,
            residue_number=pdb_resno,
            insertion_code=icode,
            residue_name=str(getattr(residue, "name", "")).strip().upper(),
            uniprot_id=uniprot_id,
            uniprot_residue_number=uniprot_resno,
        )
        residues[uniprot_resno] = ResidueObservation(
            key=key,
            atoms=AtomCoordinates(**atom_coords),
            all_heavy_atoms=heavy,
            side_chain_heavy_atoms=side_chain,
            alternate_location_used=altloc_used,
        )
    return ChainObservation(
        pdb_id=pdb_id.upper(),
        model_id=model_id,
        chain_id=chain_id,
        uniprot_id=uniprot_id,
        residues=residues,
        experimental_method=experimental_method,
        resolution_angstrom=resolution_angstrom,
    )


def find_ligand_heavy_atoms(
    path: Path,
    chain_id: str,
    ligand_name: str,
    ligand_residue: int,
    insertion_code: str = "",
    model_id: int = 1,
) -> list[list[float]]:
    structure = read_structure(path)
    model = get_model(structure, model_id)
    candidates: list[tuple[str, int, str, str, list[list[float]]]] = []
    for chain in model:
        for residue in chain:
            name = str(getattr(residue, "name", "")).strip().upper()
            if is_amino_acid(residue) or name in {"HOH", "WAT"}:
                continue
            resno, icode = residue_number_and_icode(residue)
            chosen, _ = choose_altloc_atoms(list(residue))
            heavy = [coord_list(atom) for atom in chosen.values() if not is_hydrogen(atom)]
            candidates.append((name, resno, icode, getattr(chain, "name", ""), heavy))
    matches = [
        heavy
        for name, resno, icode, cid, heavy in candidates
        if name == ligand_name.upper()
        and cid == chain_id
        and resno == ligand_residue
        and (insertion_code == "" or icode == insertion_code)
    ]
    if len(matches) == 1:
        return matches[0]
    similar = [
        f"- {name}, chain {cid}, residue {resno}{icode}".rstrip()
        for name, resno, icode, cid, _ in candidates
        if name == ligand_name.upper() or resno == ligand_residue
    ]
    detail = "\n".join(similar[:20]) if similar else "No ligand-like residues were found."
    raise LigandNotFoundError(
        f"ligand {ligand_name} {chain_id} {ligand_residue}{insertion_code} was not found uniquely.\n\n"
        f"Possible matching ligands:\n{detail}"
    )
