"""Whole-chain structural alignment."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from flexres.exceptions import AlignmentError
from flexres.models import AlignmentResult, ChainObservation, ResidueObservation
from flexres.structures.atoms import AA3_TO_1


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Return rotation, translation, and RMSD to superpose mobile onto target."""
    if mobile.shape != target.shape or mobile.ndim != 2 or mobile.shape[1] != 3 or mobile.shape[0] == 0:
        raise AlignmentError("cannot align empty or mismatched coordinate arrays")
    mobile_centroid = mobile.mean(axis=0)
    target_centroid = target.mean(axis=0)
    m0 = mobile - mobile_centroid
    t0 = target - target_centroid
    covariance = m0.T @ t0
    v, _, wt = np.linalg.svd(covariance)
    d = np.sign(np.linalg.det(v @ wt))
    correction = np.diag([1.0, 1.0, d])
    rotation = v @ correction @ wt
    translation = target_centroid - mobile_centroid @ rotation
    transformed = mobile @ rotation + translation
    rmsd = float(np.sqrt(np.mean(np.sum((transformed - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def common_ca_arrays(reference: ChainObservation, comparison: ChainObservation) -> tuple[np.ndarray, np.ndarray, list[int]]:
    common = sorted(set(reference.residues) & set(comparison.residues))
    ref_coords: list[list[float]] = []
    comp_coords: list[list[float]] = []
    used: list[int] = []
    for uniprot_resno in common:
        ref_ca = reference.residues[uniprot_resno].atoms.ca
        comp_ca = comparison.residues[uniprot_resno].atoms.ca
        if ref_ca is None or comp_ca is None:
            continue
        ref_coords.append(ref_ca)
        comp_coords.append(comp_ca)
        used.append(uniprot_resno)
    if not used:
        raise AlignmentError("no common mapped C-alpha atoms")
    return np.asarray(ref_coords, dtype=float), np.asarray(comp_coords, dtype=float), used


def align_chain_kabsch(reference: ChainObservation, comparison: ChainObservation) -> tuple[ChainObservation, AlignmentResult]:
    """Align comparison to reference using common mapped C-alpha atoms."""
    target, mobile, used = common_ca_arrays(reference, comparison)
    rotation, translation, rmsd = kabsch(mobile, target)
    transformed_residues = {
        resno: transform_residue(residue, rotation, translation)
        for resno, residue in comparison.residues.items()
    }
    aligned = comparison.model_copy(update={"residues": transformed_residues})
    result = AlignmentResult(
        rotation=rotation.tolist(),
        translation=translation.tolist(),
        global_ca_rmsd_angstrom=rmsd,
        aligned_ca_count=len(used),
        aligned_residue_count=len(used),
        backend="kabsch",
    )
    return aligned, result


def align_chain_biopython(
    reference: ChainObservation,
    comparison: ChainObservation,
    outlier_cutoff_angstrom: float = 2.0,
    report_ca_cutoff_angstrom: float = 5.0,
    max_cycles: int = 5,
) -> tuple[ChainObservation, AlignmentResult]:
    """Align comparison to reference using exact residue pairs and Bio.SVDSuperimposer."""
    try:
        from Bio.SVDSuperimposer import SVDSuperimposer
    except ImportError as exc:
        raise AlignmentError("Biopython is required for --alignment-backend biopython") from exc

    candidate_pairs = _author_residue_pairs(reference, comparison)
    fit_pairs = _fit_pairs(candidate_pairs)
    if len(fit_pairs) < 3:
        candidate_pairs = _mapped_residue_pairs(reference, comparison)
        fit_pairs = _fit_pairs(candidate_pairs)
    if len(fit_pairs) < 3:
        candidate_pairs = _sequence_same_residue_pairs(reference, comparison)
        fit_pairs = _fit_pairs(candidate_pairs)
    if not fit_pairs:
        raise AlignmentError("no same-residue C-alpha atoms after sequence alignment")

    kept_indices = np.arange(len(fit_pairs))
    target_all = np.asarray([pair[0].atoms.ca for pair in fit_pairs], dtype=float)
    mobile_all = np.asarray([pair[1].atoms.ca for pair in fit_pairs], dtype=float)
    rotation, translation, rmsd_before = _svd_fit(SVDSuperimposer, mobile_all, target_all)
    rmsd_after = rmsd_before
    cycles = 0
    for cycle in range(1, max_cycles + 1):
        target = np.asarray([fit_pairs[index][0].atoms.ca for index in kept_indices], dtype=float)
        mobile = np.asarray([fit_pairs[index][1].atoms.ca for index in kept_indices], dtype=float)
        rotation, translation, rmsd_after = _svd_fit(SVDSuperimposer, mobile, target)
        transformed = mobile @ rotation + translation
        distances = np.sqrt(np.sum((transformed - target) ** 2, axis=1))
        next_indices = kept_indices[distances <= outlier_cutoff_angstrom]
        if len(next_indices) < 3 or len(next_indices) == len(kept_indices):
            break
        kept_indices = next_indices
        cycles = cycle

    transformed_residues = {
        resno: transform_residue(residue, rotation, translation)
        for resno, residue in comparison.residues.items()
    }
    aligned = comparison.model_copy(update={"residues": transformed_residues})
    residue_pairs = [
        _aligned_pair(reference.chain_id, comparison.chain_id, ref, comp)
        for ref, comp in candidate_pairs
        if _ca_distance_after_transform(ref, comp, rotation, translation) <= report_ca_cutoff_angstrom
    ]
    result = AlignmentResult(
        rotation=rotation.tolist(),
        translation=translation.tolist(),
        global_ca_rmsd_angstrom=float(rmsd_after),
        aligned_ca_count=int(len(kept_indices)),
        alignment_rmsd_after_outlier_rejection=float(rmsd_after),
        aligned_atom_count_after_outlier_rejection=int(len(kept_indices)),
        alignment_cycles=cycles,
        alignment_rmsd_before_outlier_rejection=float(rmsd_before),
        aligned_atom_count_before_outlier_rejection=int(len(fit_pairs)),
        aligned_residue_count=len(residue_pairs),
        aligned_residue_pairs=residue_pairs,
        backend="biopython",
    )
    return aligned, result


def _chain_sequence(chain: ChainObservation) -> tuple[list[ResidueObservation], str]:
    residues = sorted(chain.residues.values(), key=lambda residue: (residue.key.residue_number, residue.key.insertion_code))
    return residues, "".join(_one_letter(residue) for residue in residues)


def _author_residue_pairs(reference: ChainObservation, comparison: ChainObservation) -> list[tuple[ResidueObservation, ResidueObservation]]:
    comp_by_author = {(_author_key(residue)): residue for residue in comparison.residues.values()}
    pairs = []
    for ref in sorted(reference.residues.values(), key=lambda residue: (residue.key.residue_number, residue.key.insertion_code)):
        comp = comp_by_author.get(_author_key(ref))
        if comp is not None:
            pairs.append((ref, comp))
    return pairs


def _mapped_residue_pairs(reference: ChainObservation, comparison: ChainObservation) -> list[tuple[ResidueObservation, ResidueObservation]]:
    return [(reference.residues[resno], comparison.residues[resno]) for resno in sorted(set(reference.residues) & set(comparison.residues))]


def _fit_pairs(candidate_pairs: list[tuple[ResidueObservation, ResidueObservation]]) -> list[tuple[ResidueObservation, ResidueObservation]]:
    return [
        (ref, comp)
        for ref, comp in candidate_pairs
        if _one_letter(ref) == _one_letter(comp) and _one_letter(ref) != "X" and ref.atoms.ca is not None and comp.atoms.ca is not None
    ]


def _sequence_same_residue_pairs(reference: ChainObservation, comparison: ChainObservation) -> list[tuple[ResidueObservation, ResidueObservation]]:
    try:
        from Bio.Align import PairwiseAligner
    except ImportError as exc:
        raise AlignmentError("Biopython is required for sequence alignment") from exc
    ref_residues, ref_sequence = _chain_sequence(reference)
    comp_residues, comp_sequence = _chain_sequence(comparison)
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(ref_sequence, comp_sequence)
    if not alignments:
        raise AlignmentError("no sequence alignment could be computed")
    pairs = []
    for ref_block, comp_block in zip(*alignments[0].aligned):
        for ref_index, comp_index in zip(range(ref_block[0], ref_block[1]), range(comp_block[0], comp_block[1])):
            ref = ref_residues[ref_index]
            comp = comp_residues[comp_index]
            if _one_letter(ref) == _one_letter(comp) and _one_letter(ref) != "X":
                pairs.append((ref, comp))
    return pairs


def _author_key(residue: ResidueObservation) -> tuple[int, str]:
    return residue.key.residue_number, residue.key.insertion_code


def _one_letter(residue: ResidueObservation) -> str:
    return AA3_TO_1.get(residue.key.residue_name.upper(), "X")


def _svd_fit(superimposer_cls, mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    sup = superimposer_cls()
    sup.set(target, mobile)
    sup.run()
    rotation, translation = sup.get_rotran()
    return np.asarray(rotation, dtype=float), np.asarray(translation, dtype=float), float(sup.get_rms())


def _aligned_pair(reference_chain_id: str, comparison_chain_id: str, ref: ResidueObservation, comp: ResidueObservation) -> dict[str, str | int]:
    return {
        "reference_chain": reference_chain_id,
        "reference_residue_number": ref.key.residue_number,
        "reference_residue_name": ref.key.residue_name,
        "comparison_chain": comparison_chain_id,
        "comparison_residue_number": comp.key.residue_number,
        "comparison_residue_name": comp.key.residue_name,
    }


def _ca_distance_after_transform(ref: ResidueObservation, comp: ResidueObservation, rotation: np.ndarray, translation: np.ndarray) -> float:
    if ref.atoms.ca is None or comp.atoms.ca is None:
        return float("inf")
    transformed = np.asarray(comp.atoms.ca, dtype=float) @ rotation + translation
    return float(np.linalg.norm(transformed - np.asarray(ref.atoms.ca, dtype=float)))


def align_chain_pymol(
    reference_path: Path,
    comparison_path: Path,
    reference_chain_id: str,
    comparison_chain_id: str,
    parse_aligned_comparison,
    pymol_executable: str = "pymol",
) -> tuple[ChainObservation, AlignmentResult]:
    """Align comparison to reference with PyMOL cmd.align and parse transformed coordinates."""
    pymol_executable = _resolve_pymol_executable(pymol_executable)
    with tempfile.TemporaryDirectory(prefix="flexres_pymol_") as temp_dir:
        temp = Path(temp_dir)
        script_path = temp / "align.py"
        aligned_path = temp / "aligned_mobile.cif"
        metrics_path = temp / "metrics.json"
        script_path.write_text(
            _pymol_script(
                _path_for_pymol(reference_path, pymol_executable),
                _path_for_pymol(comparison_path, pymol_executable),
                reference_chain_id,
                comparison_chain_id,
                _path_for_pymol(aligned_path, pymol_executable),
                _path_for_pymol(metrics_path, pymol_executable),
            ),
            encoding="utf-8",
        )
        script_arg = _path_for_pymol(script_path, pymol_executable)
        try:
            completed = subprocess.run(
                [pymol_executable, "-cq", script_arg],
                check=False,
                capture_output=True,
                text=True,
                timeout=300,
            )
        except FileNotFoundError as exc:
            raise AlignmentError(f"PyMOL executable was not found: {pymol_executable}") from exc
        except subprocess.TimeoutExpired as exc:
            raise AlignmentError("PyMOL alignment timed out") from exc
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout or "").strip()
            raise AlignmentError(f"PyMOL alignment failed: {detail}")
        if not aligned_path.exists() or not metrics_path.exists():
            raise AlignmentError("PyMOL did not produce aligned coordinates")
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
        aligned = parse_aligned_comparison(aligned_path)
        result = AlignmentResult(
            rotation=np.eye(3).tolist(),
            translation=[0.0, 0.0, 0.0],
            global_ca_rmsd_angstrom=float(metrics["rmsd_after_refinement"]),
            aligned_ca_count=int(metrics["aligned_atom_count_after_refinement"]),
            alignment_rmsd_after_outlier_rejection=float(metrics["rmsd_after_refinement"]),
            aligned_atom_count_after_outlier_rejection=int(metrics["aligned_atom_count_after_refinement"]),
            alignment_cycles=int(metrics["cycles"]),
            alignment_rmsd_before_outlier_rejection=float(metrics["rmsd_before_refinement"]),
            aligned_atom_count_before_outlier_rejection=int(metrics["aligned_atom_count_before_refinement"]),
            raw_alignment_score=float(metrics["raw_alignment_score"]),
            aligned_residue_count=int(metrics["aligned_residue_count"]),
            aligned_residue_pairs=metrics.get("aligned_residue_pairs", []),
            backend="pymol",
        )
        return aligned, result


def _path_for_pymol(path: Path, pymol_executable: str) -> str:
    """Return a path string PyMOL can open, converting to Windows paths for .exe launchers."""
    if not pymol_executable.lower().endswith(".exe"):
        return str(path)
    try:
        converted = subprocess.run(
            ["wslpath", "-w", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return str(path)
    return converted.stdout.strip()


def _resolve_pymol_executable(pymol_executable: str) -> str:
    if Path(pymol_executable).parts and (Path(pymol_executable).is_absolute() or "/" in pymol_executable or "\\" in pymol_executable):
        return pymol_executable
    found = shutil.which(pymol_executable)
    if found:
        return found
    env_bin = Path(sys.executable).parent / pymol_executable
    if env_bin.exists():
        return str(env_bin)
    return pymol_executable


def _pymol_script(
    reference_path: Path,
    comparison_path: Path,
    reference_chain_id: str,
    comparison_chain_id: str,
    aligned_path: Path,
    metrics_path: Path,
) -> str:
    return f"""
from __future__ import annotations

import json
from pymol import cmd

cmd.load({str(reference_path)!r}, "ref")
cmd.load({str(comparison_path)!r}, "mobile")
result = cmd.align(
    "mobile and chain {comparison_chain_id} and polymer and name CA",
    "ref and chain {reference_chain_id} and polymer and name CA",
    quiet=1,
    object="aln",
    reset=1,
)
ref_atoms = {{}}
mobile_atoms = {{}}
cmd.iterate("ref and chain {reference_chain_id} and polymer and name CA", "ref_atoms[index]=(chain, resi, resn)", space={{"ref_atoms": ref_atoms}})
cmd.iterate("mobile and chain {comparison_chain_id} and polymer and name CA", "mobile_atoms[index]=(chain, resi, resn)", space={{"mobile_atoms": mobile_atoms}})
pairs = []
for column in cmd.get_raw_alignment("aln"):
    if len(column) != 2:
        continue
    first, second = column
    atoms = {{first[0]: first[1], second[0]: second[1]}}
    if "ref" not in atoms or "mobile" not in atoms:
        continue
    ref_atom = ref_atoms.get(atoms["ref"])
    mobile_atom = mobile_atoms.get(atoms["mobile"])
    if ref_atom is None or mobile_atom is None:
        continue
    pairs.append({{
        "reference_chain": ref_atom[0],
        "reference_residue_number": ref_atom[1],
        "reference_residue_name": ref_atom[2],
        "comparison_chain": mobile_atom[0],
        "comparison_residue_number": mobile_atom[1],
        "comparison_residue_name": mobile_atom[2],
    }})
cmd.save({str(aligned_path)!r}, "mobile")
metrics = {{
    "rmsd_after_refinement": result[0],
    "aligned_atom_count_after_refinement": result[1],
    "cycles": result[2],
    "rmsd_before_refinement": result[3],
    "aligned_atom_count_before_refinement": result[4],
    "raw_alignment_score": result[5],
    "aligned_residue_count": result[6],
    "aligned_residue_pairs": pairs,
}}
with open({str(metrics_path)!r}, "w", encoding="utf-8") as handle:
    json.dump(metrics, handle)
cmd.quit()
"""


def transform_point(point: list[float] | None, rotation: np.ndarray, translation: np.ndarray) -> list[float] | None:
    if point is None:
        return None
    return (np.asarray(point, dtype=float) @ rotation + translation).astype(float).tolist()


def transform_residue(residue: ResidueObservation, rotation: np.ndarray, translation: np.ndarray) -> ResidueObservation:
    atoms = residue.atoms.model_copy(
        update={
            "n": transform_point(residue.atoms.n, rotation, translation),
            "ca": transform_point(residue.atoms.ca, rotation, translation),
            "c": transform_point(residue.atoms.c, rotation, translation),
            "o": transform_point(residue.atoms.o, rotation, translation),
        }
    )
    heavy = [transform_point(point, rotation, translation) for point in residue.all_heavy_atoms]
    side_chain = {
        name: transformed
        for name, point in residue.side_chain_heavy_atoms.items()
        if (transformed := transform_point(point, rotation, translation)) is not None
    }
    return residue.model_copy(update={"atoms": atoms, "all_heavy_atoms": [p for p in heavy if p is not None], "side_chain_heavy_atoms": side_chain})
