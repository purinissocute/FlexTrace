from __future__ import annotations

from flexres.structures.parser import residue_label_number


class Residue:
    label_seq = 4


def test_residue_label_number_returns_gemmi_label_seq() -> None:
    assert residue_label_number(Residue()) == 4
