"""PDBe/SIFTS residue mapping client."""

from __future__ import annotations

from typing import Any

import requests

from flexres.cache.manager import CacheManager
from flexres.config import PDBe_API, USER_AGENT
from flexres.api.retry import retry, stop_after_attempt, wait_exponential


class SIFTSClient:
    """Fetch residue-level UniProt mappings from PDBe."""

    def __init__(self, cache: CacheManager, timeout: float = 30.0) -> None:
        self.cache = cache
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": USER_AGENT})

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=8))
    def _get_json(self, url: str) -> dict[str, Any]:
        response = self.session.get(url, timeout=self.timeout)
        response.raise_for_status()
        return response.json()

    def mappings(self, pdb_id: str, force: bool = False) -> dict[str, Any]:
        pdb_id = pdb_id.lower()
        path = self.cache.mapping_path(pdb_id)
        if not force and self.cache.is_valid_file(path):
            return self.cache.read_json(path)
        payload = self._get_json(f"{PDBe_API}/mappings/uniprot/{pdb_id}")
        self.cache.atomic_write_json(path, payload)
        return payload

    def chain_residue_map(self, pdb_id: str, chain_id: str, force: bool = False) -> dict[tuple[int, str], tuple[str, int]]:
        payload = self.mappings(pdb_id, force=force)
        entry = payload.get(pdb_id.lower(), {}).get("UniProt", {})
        result: dict[tuple[int, str], tuple[str, int]] = {}
        for accession, acc_payload in entry.items():
            for mapping in acc_payload.get("mappings", []):
                if mapping.get("chain_id") != chain_id and mapping.get("struct_asym_id") != chain_id:
                    continue
                start = mapping.get("start", {})
                end = mapping.get("end", {})
                label_start = int(start.get("residue_number"))
                label_end = int(end.get("residue_number"))
                pdb_start = _author_or_label_start(start, end, label_start, label_end)
                pdb_end = int(end.get("author_residue_number") or label_end)
                uni_start = int(mapping.get("unp_start"))
                step = 1 if pdb_end >= pdb_start else -1
                for offset, pdb_resno in enumerate(range(pdb_start, pdb_end + step, step)):
                    result[(pdb_resno, "")] = (accession, uni_start + offset)
        return result

    def chain_uniprot_accessions(self, pdb_id: str, chain_id: str, force: bool = False) -> set[str]:
        payload = self.mappings(pdb_id, force=force)
        entry = payload.get(pdb_id.lower(), {}).get("UniProt", {})
        accessions: set[str] = set()
        for accession, acc_payload in entry.items():
            for mapping in acc_payload.get("mappings", []):
                if mapping.get("chain_id") == chain_id or mapping.get("struct_asym_id") == chain_id:
                    accessions.add(accession)
        return accessions

def _author_or_label_start(start: dict[str, Any], end: dict[str, Any], label_start: int, label_end: int) -> int:
    author_start = start.get("author_residue_number")
    if author_start is not None:
        return int(author_start)
    author_end = end.get("author_residue_number")
    if author_end is not None:
        label_length = abs(label_end - label_start)
        return int(author_end) - label_length if label_end >= label_start else int(author_end) + label_length
    return label_start
