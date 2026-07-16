from __future__ import annotations

import pytest

from flexres.api.rcsb import validate_pdb_id


def test_validate_pdb_id_normalizes() -> None:
    assert validate_pdb_id("2q8f") == "2Q8F"


@pytest.mark.parametrize("bad", ["", "123", "12345", "AB-C"])
def test_validate_pdb_id_rejects_bad_values(bad: str) -> None:
    with pytest.raises(ValueError):
        validate_pdb_id(bad)
