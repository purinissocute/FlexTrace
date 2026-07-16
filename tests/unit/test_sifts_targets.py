from __future__ import annotations

from flexres.api.sifts import SIFTSClient


def test_chain_residue_map_prefers_author_numbering() -> None:
    client = SIFTSClient.__new__(SIFTSClient)
    client.mappings = lambda pdb_id, force=False: {  # type: ignore[method-assign]
        "8c1d": {
            "UniProt": {
                "O14965": {
                    "mappings": [
                        {
                            "chain_id": "A",
                            "struct_asym_id": "A",
                            "unp_start": 126,
                            "start": {"author_residue_number": None, "residue_number": 3},
                            "end": {"author_residue_number": 388, "residue_number": 265},
                        }
                    ]
                }
            }
        }
    }
    mapping = client.chain_residue_map("8C1D", "A")
    assert mapping[(126, "")] == ("O14965", 126)
    assert mapping[(199, "")] == ("O14965", 199)
    assert mapping[(388, "")] == ("O14965", 388)
