from __future__ import annotations

from pathlib import Path

import pytest

from flexres.analysis.alignment import _pymol_script, align_chain_pymol
from flexres.exceptions import AlignmentError


def test_pymol_script_uses_cmd_align() -> None:
    script = _pymol_script(
        Path("ref.cif"),
        Path("mobile.cif"),
        "A",
        "B",
        Path("aligned.cif"),
        Path("metrics.json"),
    )
    assert "cmd.align(" in script
    assert "mobile and chain B and polymer and name CA" in script
    assert "ref and chain A and polymer and name CA" in script
    assert "cmd.save(" in script


def test_align_chain_pymol_reports_missing_executable(tmp_path) -> None:
    with pytest.raises(AlignmentError, match="PyMOL executable was not found"):
        align_chain_pymol(
            tmp_path / "ref.cif",
            tmp_path / "mobile.cif",
            "A",
            "A",
            lambda _path: None,
            "definitely-not-a-pymol-binary",
        )
