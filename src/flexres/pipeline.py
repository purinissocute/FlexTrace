"""End-to-end analysis pipeline."""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from flexres.analysis.alignment import align_chain_biopython, align_chain_kabsch, align_chain_pymol
from flexres.analysis.comparisons import compare_pymol_aligned_residues, compare_residues
from flexres.analysis.pocket import detect_pocket_residues
from flexres.analysis.summaries import summarize_residue_comparisons
from flexres.api.rcsb import RCSBClient, validate_pdb_id
from flexres.api.sifts import SIFTSClient
from flexres.api.uniprot import UniProtClient
from flexres.cache.manager import CacheManager
from flexres.exceptions import AlignmentError, FlexresError
from flexres.mapping.reference import reference_uniprot, validate_reference_chain
from flexres.mapping.sequence import canonical_mismatch_positions
from flexres.models import AnalysisSettings, ComparisonTarget, SkippedComparison, SkippedResidue
from flexres.output.csv_writer import write_csv
from flexres.output.json_writer import write_json
from flexres.output.manifest import build_manifest
from flexres.output.schemas import (
    RESIDUE_COMPARISON_COLUMNS,
    SKIPPED_COMPARISON_COLUMNS,
    SKIPPED_RESIDUE_COLUMNS,
    flatten_comparison,
    flatten_skipped_comparison,
    flatten_skipped_residue,
)
from flexres.structures.parser import find_ligand_heavy_atoms, parse_chain_observation, read_structure

LOGGER = logging.getLogger("flexres")


def run_analysis(settings: AnalysisSettings) -> dict:
    settings.pdb_id = validate_pdb_id(settings.pdb_id)
    settings.output_dir.mkdir(parents=True, exist_ok=True)
    _configure_logging(settings.output_dir / "flexres.log", settings.log_level)
    cache = CacheManager(settings.cache_dir)
    if settings.cache_dir is None:
        settings.cache_dir = cache.root
    rcsb = RCSBClient(cache)
    sifts = SIFTSClient(cache)
    uniprot = UniProtClient(cache)

    LOGGER.info("Downloading reference %s", settings.pdb_id)
    ref_path = rcsb.download_mmcif(settings.pdb_id, force=settings.force_redownload)
    validate_reference_chain(ref_path, settings.pdb_id, settings.chain_id)
    uniprot_id = reference_uniprot(sifts, settings.pdb_id, settings.chain_id, force=settings.force_redownload)
    ref_map = sifts.chain_residue_map(settings.pdb_id, settings.chain_id, force=settings.force_redownload)
    reference = parse_chain_observation(ref_path, settings.pdb_id, settings.chain_id, uniprot_id, ref_map)
    canonical_sequence = uniprot.canonical_sequence(uniprot_id, force=settings.force_redownload) if settings.exclude_mutations else None
    reference_mismatches: list[tuple[int, str]] = []
    if canonical_sequence is not None:
        reference_mismatches = canonical_mismatch_positions(reference, canonical_sequence)
        if reference_mismatches:
            LOGGER.warning(
                "Reference chain differs from canonical UniProt sequence at mapped residue(s): %s",
                ", ".join(label for _, label in reference_mismatches),
            )

    if settings.all_residues:
        pocket_distances: dict[int, float] = {}
        residue_numbers = set(reference.residues)
    else:
        if not (settings.ligand_name and settings.ligand_chain and settings.ligand_residue is not None):
            raise FlexresError("pocket mode requires --ligand-name, --ligand-chain and --ligand-residue")
        ligand_atoms = find_ligand_heavy_atoms(
            ref_path,
            settings.ligand_chain,
            settings.ligand_name,
            settings.ligand_residue,
            settings.ligand_insertion_code,
        )
        pocket_distances = detect_pocket_residues(reference, ligand_atoms, settings.pocket_radius)
        residue_numbers = set(pocket_distances)
    reference_mismatch_residues = {resno for resno, _ in reference_mismatches}
    if reference_mismatch_residues:
        residue_numbers -= reference_mismatch_residues

    targets = rcsb.search_uniprot(uniprot_id, force=settings.force_redownload)
    if settings.include_reference_entry:
        targets.append(ComparisonTarget(pdb_id=settings.pdb_id, chain_id=settings.chain_id))
    targets = _dedupe_targets(targets)

    comparisons = []
    skipped_comparisons: list[SkippedComparison] = []
    skipped_residues = [
        SkippedResidue(
            reference_pdb_id=settings.pdb_id,
            reference_chain=settings.chain_id,
            reference_residue_number=reference.residues[resno].key.residue_number if resno in reference.residues else None,
            uniprot_residue_number=resno,
            reason="reference_sequence_mismatch_mutation",
            missing_reference_atoms=label,
        )
        for resno, label in reference_mismatches
    ]
    alignments = []

    with ThreadPoolExecutor(max_workers=max(1, settings.max_workers)) as executor:
        futures = {
            executor.submit(_analyse_target, settings, cache.root, target, reference, uniprot_id, canonical_sequence, residue_numbers, pocket_distances): target
            for target in targets
        }
        for future in as_completed(futures):
            target = futures[future]
            try:
                result = future.result()
            except Exception as exc:  # noqa: BLE001
                skipped_comparisons.append(
                    SkippedComparison(
                        pdb_id=target.pdb_id,
                        chain=target.chain_id,
                        stage="comparison",
                        reason="unsupported_structure",
                        details=str(exc),
                    )
                )
                continue
            comparisons.extend(result["comparisons"])
            skipped_residues.extend(result["skipped_residues"])
            skipped_comparisons.extend(result["skipped_comparisons"])
            alignments.extend(result["alignments"])

    comparison_rows = sorted(
        [flatten_comparison(row) for row in comparisons if row.side_chain_rmsd_angstrom is not None],
        key=lambda row: (-row["side_chain_rmsd_angstrom"], row["comparison_pdb_id"], row["comparison_model"], row["comparison_chain"], row["uniprot_residue_number"]),
    )
    summary_rows = summarize_residue_comparisons(comparison_rows)
    skipped_comparison_rows = sorted([flatten_skipped_comparison(row) for row in skipped_comparisons], key=lambda r: (r["pdb_id"], r.get("model") or 0, r.get("chain") or "", r["reason"]))
    skipped_residue_rows = sorted([flatten_skipped_residue(row) for row in skipped_residues], key=lambda r: (r.get("comparison_pdb_id") or "", r.get("comparison_model") or 0, r.get("comparison_chain") or "", r.get("uniprot_residue_number") or -1, r["reason"]))

    if "csv" in settings.output_formats:
        write_csv(settings.output_dir / "residue_comparisons.csv", comparison_rows, RESIDUE_COMPARISON_COLUMNS)
        write_csv(settings.output_dir / "residue_summary.csv", summary_rows)
        write_csv(settings.output_dir / "skipped_comparisons.csv", skipped_comparison_rows, SKIPPED_COMPARISON_COLUMNS)
        write_csv(settings.output_dir / "skipped_residues.csv", skipped_residue_rows, SKIPPED_RESIDUE_COLUMNS)

    counts = {
        "number_of_pdb_entries_found": len({t.pdb_id for t in targets}),
        "number_of_chains_found": len(targets),
        "number_of_models_analysed": len({(r["comparison_pdb_id"], r["comparison_model"], r["comparison_chain"]) for r in comparison_rows}),
        "number_of_comparisons_completed": len({(r["comparison_pdb_id"], r["comparison_model"], r["comparison_chain"]) for r in comparison_rows}),
        "number_of_comparisons_skipped": len(skipped_comparison_rows),
    }
    manifest = build_manifest(settings, uniprot_id, counts)
    analysis_json = {
        "reference": {"pdb_id": settings.pdb_id, "chain": settings.chain_id, "model": 1},
        "uniprot_mapping": {"uniprot_id": uniprot_id},
        "pocket": {"enabled": not settings.all_residues, "radius_angstrom": settings.pocket_radius, "residues": sorted(pocket_distances)},
        "settings": settings.model_dump(mode="json"),
        "retrieval_summary": {"targets": [t.model_dump() for t in targets]},
        "comparison_summary": counts,
        "alignments": alignments,
        "residue_comparisons": comparison_rows,
        "residue_summary": summary_rows,
        "skipped_comparisons": skipped_comparison_rows,
        "skipped_residues": skipped_residue_rows,
    }
    if "json" in settings.output_formats:
        write_json(settings.output_dir / "analysis.json", analysis_json)
        write_json(settings.output_dir / "manifest.json", manifest)
    return analysis_json


def _analyse_target(
    settings: AnalysisSettings,
    cache_root: Path,
    target: ComparisonTarget,
    reference,
    uniprot_id: str,
    canonical_sequence: str | None,
    residue_numbers: set[int],
    pocket_distances: dict[int, float],
) -> dict:
    cache = CacheManager(cache_root)
    rcsb = RCSBClient(cache)
    sifts = SIFTSClient(cache)
    skipped_comparisons: list[SkippedComparison] = []
    skipped_residues = []
    comparisons = []
    alignments = []
    try:
        path = rcsb.download_mmcif(target.pdb_id, force=settings.force_redownload)
        mapping = sifts.chain_residue_map(target.pdb_id, target.chain_id, force=settings.force_redownload)
        read_structure(path)
    except Exception as exc:  # noqa: BLE001
        return {"comparisons": [], "skipped_residues": [], "skipped_comparisons": [SkippedComparison(pdb_id=target.pdb_id, chain=target.chain_id, stage="download_parse", reason="parse_failed", details=str(exc))], "alignments": []}
    for model_id in [1]:
        if target.pdb_id.upper() == settings.pdb_id.upper() and target.chain_id == settings.chain_id and model_id == 1:
            continue
        try:
            comp = parse_chain_observation(path, target.pdb_id, target.chain_id, uniprot_id, mapping, model_id=model_id)
            mutated_residues: dict[int, str] = {}
            if canonical_sequence is not None:
                mutated_residues = dict(canonical_mismatch_positions(comp, canonical_sequence, limit=10_000))
            comp = comp.model_copy(update={
                "polymer_entity_id": target.entity_id,
                "author_chain_id": target.author_chain_id,
                "sequence_identity": target.sequence_identity,
                "comparison_sequence_coverage": target.sequence_coverage,
                "experimental_method": target.experimental_method,
                "resolution_angstrom": target.resolution_angstrom,
                "release_date": target.release_date,
            })
            if settings.alignment_backend == "pymol":
                aligned, alignment = align_chain_pymol(
                    ref_path_for_worker(settings, cache_root),
                    path,
                    settings.chain_id,
                    target.chain_id,
                    lambda aligned_path: parse_chain_observation(aligned_path, target.pdb_id, target.chain_id, uniprot_id, mapping, model_id=model_id),
                    settings.pymol_executable,
                )
                aligned = aligned.model_copy(update={
                    "polymer_entity_id": target.entity_id,
                    "author_chain_id": target.author_chain_id,
                    "sequence_identity": target.sequence_identity,
                    "comparison_sequence_coverage": target.sequence_coverage,
                    "experimental_method": target.experimental_method,
                    "resolution_angstrom": target.resolution_angstrom,
                    "release_date": target.release_date,
                })
            elif settings.alignment_backend == "biopython":
                aligned, alignment = align_chain_biopython(reference, comp)
            else:
                aligned, alignment = align_chain_kabsch(reference, comp)
            if alignment.global_ca_rmsd_angstrom > settings.alignment_rmsd_cutoff:
                skipped_comparisons.append(
                    SkippedComparison(
                        pdb_id=target.pdb_id,
                        model=model_id,
                        chain=target.chain_id,
                        stage="alignment",
                        reason="poor_global_alignment",
                        details=(
                            f"global C-alpha RMSD {alignment.global_ca_rmsd_angstrom:.4f} A "
                            f"exceeds cutoff {settings.alignment_rmsd_cutoff:.4f} A"
                        ),
                    )
                )
                continue
            if settings.alignment_backend in {"pymol", "biopython"}:
                rows, skips = compare_pymol_aligned_residues(reference, aligned, alignment, residue_numbers, pocket_distances)
            else:
                rows, skips = compare_residues(reference, aligned, alignment, residue_numbers, pocket_distances, mutated_residues)
            comparisons.extend(rows)
            skipped_residues.extend(skips)
            alignments.append({"pdb_id": target.pdb_id, "chain": target.chain_id, "model": model_id, **alignment.model_dump()})
        except AlignmentError as exc:
            skipped_comparisons.append(SkippedComparison(pdb_id=target.pdb_id, model=model_id, chain=target.chain_id, stage="alignment", reason="alignment_failed", details=str(exc)))
        except Exception as exc:  # noqa: BLE001
            skipped_comparisons.append(SkippedComparison(pdb_id=target.pdb_id, model=model_id, chain=target.chain_id, stage="comparison", reason="unsupported_structure", details=str(exc)))
    return {"comparisons": comparisons, "skipped_residues": skipped_residues, "skipped_comparisons": skipped_comparisons, "alignments": alignments}


def _dedupe_targets(targets: list[ComparisonTarget]) -> list[ComparisonTarget]:
    seen = set()
    result = []
    for target in sorted(targets, key=lambda t: (t.pdb_id, t.chain_id, t.entity_id or "")):
        key = (target.pdb_id.upper(), target.chain_id)
        if key in seen:
            continue
        seen.add(key)
        result.append(target)
    return result


def ref_path_for_worker(settings: AnalysisSettings, cache_root: Path) -> Path:
    return CacheManager(cache_root).structure_path(settings.pdb_id)


def _configure_logging(path: Path, level: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=path,
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        force=True,
    )
