"""Atom and residue coordinate utilities."""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Iterable

import numpy as np

BACKBONE_ATOMS = ("N", "CA", "C", "O")
AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "MSE",
}
AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
    "MSE": "M",
}


def clean_gemmi_text(value: object) -> str:
    """Normalize Gemmi blank/null text markers to an empty string."""
    text = str(value or "").replace("\x00", "").strip()
    return "" if text in {".", "?"} else text


def is_hydrogen(atom: object) -> bool:
    element = getattr(getattr(atom, "element", None), "name", "")
    name = getattr(atom, "name", "").strip()
    return element == "H" or name.startswith("H")


def is_amino_acid(residue: object) -> bool:
    return getattr(residue, "name", "").strip().upper() in AA3


def coord_list(atom: object) -> list[float]:
    pos = getattr(atom, "pos")
    return [float(pos.x), float(pos.y), float(pos.z)]


def choose_altloc_atoms(atoms: Iterable[object]) -> tuple[dict[str, object], str | None]:
    """Choose one deterministic atom for each atom name."""
    grouped: dict[str, list[object]] = defaultdict(list)
    for atom in atoms:
        grouped[getattr(atom, "name", "").strip()].append(atom)
    chosen: dict[str, object] = {}
    altloc_used: str | None = None
    for name, candidates in grouped.items():
        def key(atom: object) -> tuple[float, int, str]:
            alt = clean_gemmi_text(getattr(atom, "altloc", ""))
            alt_rank = 2 if alt == "" else 1 if alt == "A" else 0
            return (float(getattr(atom, "occ", 0.0) or 0.0), alt_rank, "".join(chr(255 - ord(c)) for c in alt))

        selected = sorted(candidates, key=key, reverse=True)[0]
        chosen[name] = selected
        alt = clean_gemmi_text(getattr(selected, "altloc", ""))
        if alt:
            altloc_used = alt if altloc_used is None else min(altloc_used, alt)
    return chosen, altloc_used


def distance(a: list[float], b: list[float]) -> float:
    return float(np.linalg.norm(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)))


def min_distance(points_a: Iterable[list[float]], points_b: Iterable[list[float]]) -> float:
    best = math.inf
    for a in points_a:
        for b in points_b:
            best = min(best, distance(a, b))
    return float(best)
