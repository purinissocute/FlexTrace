"""Local cache management with atomic writes."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any

from platformdirs import user_cache_dir


class CacheManager:
    """Manage cached API responses and downloaded structures."""

    def __init__(self, root: Path | None = None) -> None:
        self.root = Path(root) if root else Path(user_cache_dir("flexres"))
        self.structures = self.root / "structures"
        self.mappings = self.root / "mappings"
        self.searches = self.root / "searches"
        self.metadata = self.root / "metadata"
        for path in (self.structures, self.mappings, self.searches, self.metadata):
            path.mkdir(parents=True, exist_ok=True)

    def structure_path(self, pdb_id: str) -> Path:
        return self.structures / f"{pdb_id.lower()}.cif"

    def mapping_path(self, pdb_id: str) -> Path:
        return self.mappings / f"{pdb_id.lower()}.json"

    def search_path(self, uniprot_id: str) -> Path:
        return self.searches / f"{uniprot_id}.json"

    def uniprot_sequence_path(self, uniprot_id: str) -> Path:
        return self.metadata / f"uniprot_{uniprot_id.upper()}.fasta"

    @staticmethod
    def is_valid_file(path: Path) -> bool:
        return path.exists() and path.is_file() and path.stat().st_size > 0

    @staticmethod
    def atomic_write_bytes(path: Path, content: bytes) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
            handle.write(content)
            temp = Path(handle.name)
        temp.replace(path)

    @staticmethod
    def atomic_write_json(path: Path, payload: Any) -> None:
        data = json.dumps(payload, indent=2, sort_keys=True).encode("utf-8")
        CacheManager.atomic_write_bytes(path, data)

    @staticmethod
    def read_json(path: Path) -> Any:
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
