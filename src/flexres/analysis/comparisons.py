"""Assemble residue comparisons after global alignment."""

from __future__ import annotations

from flexres.analysis.metrics import backbone_rmsd, ca_displacement, missing_backbone_atoms, side_chain_rmsd
from flexres.models import AlignmentResult, ChainObservation, ResidueComparison, ResidueObservation, SkippedResidue


def compare_residues(
    reference: ChainObservation,
    comparison: ChainObservation,
    alignment: AlignmentResult,
    residue_numbers: set[int],
    pocket_distances: dict[int, float],
    skipped_mutation_residues: dict[int, str] | None = None,
) -> tuple[list[ResidueComparison], list[SkippedResidue]]:
    rows: list[ResidueComparison] = []
    skipped: list[SkippedResidue] = []
    skipped_mutation_residues = skipped_mutation_residues or {}
    for uniprot_resno in sorted(residue_numbers):
        ref = reference.residues.get(uniprot_resno)
        comp = comparison.residues.get(uniprot_resno)
        if uniprot_resno in skipped_mutation_residues:
            skipped.append(
                _skip(
                    reference,
                    comparison,
                    uniprot_resno,
                    "comparison_sequence_mismatch_mutation",
                    missing_comp=[skipped_mutation_residues[uniprot_resno]],
                    ref_resno=ref.key.residue_number if ref is not None else None,
                )
            )
            continue
        if ref is None:
            skipped.append(_skip(reference, comparison, uniprot_resno, "missing_uniprot_mapping"))
            continue
        if comp is None:
            skipped.append(_skip(reference, comparison, uniprot_resno, "residue_not_present_in_comparison", ref_resno=ref.key.residue_number))
            continue
        missing_ref = missing_backbone_atoms(ref.atoms)
        missing_comp = missing_backbone_atoms(comp.atoms)
        if ref.atoms.ca is None or comp.atoms.ca is None:
            skipped.append(_skip(reference, comparison, uniprot_resno, "missing_ca", missing_ref, missing_comp, ref.key.residue_number))
            continue
        if missing_ref or missing_comp:
            skipped.append(
                _skip(reference, comparison, uniprot_resno, "missing_backbone_atoms", missing_ref, missing_comp, ref.key.residue_number)
            )
            continue
        if not ref.side_chain_heavy_atoms or not comp.side_chain_heavy_atoms:
            skipped.append(_skip(reference, comparison, uniprot_resno, "missing_side_chain_heavy_atoms", ref_resno=ref.key.residue_number))
            continue
        ref_atoms = {"N": ref.atoms.n, "CA": ref.atoms.ca, "C": ref.atoms.c, "O": ref.atoms.o}
        comp_atoms = {"N": comp.atoms.n, "CA": comp.atoms.ca, "C": comp.atoms.c, "O": comp.atoms.o}
        try:
            side_chain_value = side_chain_rmsd(ref.side_chain_heavy_atoms, comp.side_chain_heavy_atoms)
        except ValueError:
            skipped.append(_skip(reference, comparison, uniprot_resno, "no_common_side_chain_heavy_atoms", ref_resno=ref.key.residue_number))
            continue
        rows.append(
            ResidueComparison(
                reference=ref.key,
                comparison=comp.key,
                ca_displacement_angstrom=ca_displacement(ref.atoms.ca, comp.atoms.ca),
                backbone_rmsd_angstrom=backbone_rmsd(ref_atoms, comp_atoms),  # type: ignore[arg-type]
                side_chain_rmsd_angstrom=side_chain_value,
                in_reference_pocket=uniprot_resno in pocket_distances,
                distance_to_reference_ligand_angstrom=pocket_distances.get(uniprot_resno),
                experimental_method=comparison.experimental_method,
                resolution_angstrom=comparison.resolution_angstrom,
                aligned_ca_count=alignment.aligned_ca_count,
                global_ca_rmsd_angstrom=alignment.global_ca_rmsd_angstrom,
                comparison_release_date=comparison.release_date,
                comparison_polymer_entity_id=comparison.polymer_entity_id,
                reference_sequence_coverage=comparison.reference_sequence_coverage,
                comparison_sequence_coverage=comparison.comparison_sequence_coverage,
                sequence_identity=comparison.sequence_identity,
                alternate_location_used=comp.alternate_location_used or ref.alternate_location_used,
            )
        )
    return sorted(rows, key=lambda r: (-(r.side_chain_rmsd_angstrom or -1.0), r.comparison.pdb_id, r.comparison.chain_id, r.reference.uniprot_residue_number)), skipped


def compare_pymol_aligned_residues(
    reference: ChainObservation,
    comparison: ChainObservation,
    alignment: AlignmentResult,
    residue_numbers: set[int],
    pocket_distances: dict[int, float],
) -> tuple[list[ResidueComparison], list[SkippedResidue]]:
    """Compare residues paired by PyMOL alignment rather than SIFTS residue numbering."""
    ref_by_resno = {_residue_auth_key(residue): residue for residue in reference.residues.values()}
    comp_by_resno = {_residue_auth_key(residue): residue for residue in comparison.residues.values()}
    rows: list[ResidueComparison] = []
    skipped: list[SkippedResidue] = []
    seen: set[int] = set()
    for pair in alignment.aligned_residue_pairs:
        try:
            ref_resno = int(str(pair["reference_residue_number"]))
            comp_resno = int(str(pair["comparison_residue_number"]))
        except (KeyError, ValueError):
            continue
        ref = ref_by_resno.get(ref_resno)
        comp = comp_by_resno.get(comp_resno)
        if ref is None or comp is None:
            continue
        uniprot_resno = ref.key.uniprot_residue_number
        if uniprot_resno not in residue_numbers or uniprot_resno in seen:
            continue
        seen.add(uniprot_resno)
        if ref.key.residue_name != comp.key.residue_name:
            skipped.append(_skip(reference, comparison, uniprot_resno, "comparison_sequence_mismatch_mutation", missing_comp=[f"{ref.key.residue_name}{uniprot_resno}{comp.key.residue_name}"], ref_resno=ref.key.residue_number))
            continue
        row, skip = _compare_residue(reference, comparison, alignment, uniprot_resno, ref, comp, pocket_distances)
        if row is not None:
            rows.append(row)
        elif skip is not None:
            skipped.append(skip)
    return sorted(rows, key=lambda r: (-(r.side_chain_rmsd_angstrom or -1.0), r.comparison.pdb_id, r.comparison.chain_id, r.reference.uniprot_residue_number)), skipped


def _compare_residue(
    reference: ChainObservation,
    comparison: ChainObservation,
    alignment: AlignmentResult,
    uniprot_resno: int,
    ref: ResidueObservation,
    comp: ResidueObservation,
    pocket_distances: dict[int, float],
) -> tuple[ResidueComparison | None, SkippedResidue | None]:
    missing_ref = missing_backbone_atoms(ref.atoms)
    missing_comp = missing_backbone_atoms(comp.atoms)
    if ref.atoms.ca is None or comp.atoms.ca is None:
        return None, _skip(reference, comparison, uniprot_resno, "missing_ca", missing_ref, missing_comp, ref.key.residue_number)
    if missing_ref or missing_comp:
        return None, _skip(reference, comparison, uniprot_resno, "missing_backbone_atoms", missing_ref, missing_comp, ref.key.residue_number)
    if not ref.side_chain_heavy_atoms or not comp.side_chain_heavy_atoms:
        return None, _skip(reference, comparison, uniprot_resno, "missing_side_chain_heavy_atoms", ref_resno=ref.key.residue_number)
    ref_atoms = {"N": ref.atoms.n, "CA": ref.atoms.ca, "C": ref.atoms.c, "O": ref.atoms.o}
    comp_atoms = {"N": comp.atoms.n, "CA": comp.atoms.ca, "C": comp.atoms.c, "O": comp.atoms.o}
    try:
        side_chain_value = side_chain_rmsd(ref.side_chain_heavy_atoms, comp.side_chain_heavy_atoms)
    except ValueError:
        return None, _skip(reference, comparison, uniprot_resno, "no_common_side_chain_heavy_atoms", ref_resno=ref.key.residue_number)
    comp_key = comp.key.model_copy(update={"uniprot_id": ref.key.uniprot_id, "uniprot_residue_number": ref.key.uniprot_residue_number})
    return ResidueComparison(
        reference=ref.key,
        comparison=comp_key,
        ca_displacement_angstrom=ca_displacement(ref.atoms.ca, comp.atoms.ca),
        backbone_rmsd_angstrom=backbone_rmsd(ref_atoms, comp_atoms),  # type: ignore[arg-type]
        side_chain_rmsd_angstrom=side_chain_value,
        in_reference_pocket=uniprot_resno in pocket_distances,
        distance_to_reference_ligand_angstrom=pocket_distances.get(uniprot_resno),
        experimental_method=comparison.experimental_method,
        resolution_angstrom=comparison.resolution_angstrom,
        aligned_ca_count=alignment.aligned_ca_count,
        global_ca_rmsd_angstrom=alignment.global_ca_rmsd_angstrom,
        comparison_release_date=comparison.release_date,
        comparison_polymer_entity_id=comparison.polymer_entity_id,
        reference_sequence_coverage=comparison.reference_sequence_coverage,
        comparison_sequence_coverage=comparison.comparison_sequence_coverage,
        sequence_identity=comparison.sequence_identity,
        alternate_location_used=comp.alternate_location_used or ref.alternate_location_used,
    ), None


def _residue_auth_key(residue: ResidueObservation) -> int:
    return residue.key.residue_number


def _skip(
    reference: ChainObservation,
    comparison: ChainObservation,
    uniprot_resno: int,
    reason: str,
    missing_ref: list[str] | None = None,
    missing_comp: list[str] | None = None,
    ref_resno: int | None = None,
) -> SkippedResidue:
    return SkippedResidue(
        reference_pdb_id=reference.pdb_id,
        reference_chain=reference.chain_id,
        reference_residue_number=ref_resno,
        uniprot_residue_number=uniprot_resno,
        comparison_pdb_id=comparison.pdb_id,
        comparison_chain=comparison.chain_id,
        comparison_model=comparison.model_id,
        reason=reason,
        missing_reference_atoms=",".join(missing_ref or []),
        missing_comparison_atoms=",".join(missing_comp or []),
    )
