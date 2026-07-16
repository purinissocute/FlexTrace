"""RCSB API and file download helpers."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import requests

from flexres.cache.manager import CacheManager
from flexres.config import RCSB_DATA_API, RCSB_FILES, RCSB_SEARCH_API, USER_AGENT
from flexres.exceptions import StructureDownloadError
from flexres.models import ComparisonTarget
from flexres.api.retry import retry, stop_after_attempt, wait_exponential

PDB_ID_RE = re.compile(r"^[A-Za-z0-9]{4}$")


def validate_pdb_id(pdb_id: str) -> str:
    """Return a normalized PDB ID or raise ValueError."""
    if not PDB_ID_RE.match(pdb_id):
        raise ValueError(f"invalid PDB ID {pdb_id!r}; expected four alphanumeric characters")
    return pdb_id.upper()


class RCSBClient:
    """Small cache-aware RCSB client."""

    def __init__(self, cache: CacheManager, timeout: float = 30.0) -> None:
        self.cache = cache
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": USER_AGENT})

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=8))
    def _get(self, url: str) -> requests.Response:
        response = self.session.get(url, timeout=self.timeout)
        response.raise_for_status()
        return response

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=8))
    def _post_json(self, url: str, payload: dict[str, Any]) -> requests.Response:
        response = self.session.post(url, json=payload, timeout=self.timeout)
        response.raise_for_status()
        return response

    def download_mmcif(self, pdb_id: str, force: bool = False) -> Path:
        pdb_id = validate_pdb_id(pdb_id)
        path = self.cache.structure_path(pdb_id)
        if not force and self.cache.is_valid_file(path):
            return path
        url = f"{RCSB_FILES}/{pdb_id}.cif"
        try:
            response = self._get(url)
        except Exception as exc:  # noqa: BLE001
            raise StructureDownloadError(f"failed to download mmCIF for {pdb_id}: {exc}") from exc
        if not response.content:
            raise StructureDownloadError(f"downloaded mmCIF for {pdb_id} was empty")
        self.cache.atomic_write_bytes(path, response.content)
        return path

    def entry_metadata(self, pdb_id: str) -> dict[str, Any]:
        pdb_id = validate_pdb_id(pdb_id)
        return self._get(f"{RCSB_DATA_API}/entry/{pdb_id}").json()

    def polymer_entity_metadata(self, pdb_id: str, entity_id: str) -> dict[str, Any]:
        pdb_id = validate_pdb_id(pdb_id)
        return self._get(f"{RCSB_DATA_API}/polymer_entity/{pdb_id}/{entity_id}").json()

    def search_uniprot(self, uniprot_id: str, force: bool = False) -> list[ComparisonTarget]:
        path = self.cache.search_path(uniprot_id)
        if not force and self.cache.is_valid_file(path):
            payload = self.cache.read_json(path)
        else:
            query = {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uniprot_id,
                    },
                },
                "return_type": "polymer_entity",
                "request_options": {"paginate": {"start": 0, "rows": 10000}},
            }
            payload = self._post_json(RCSB_SEARCH_API, query).json()
            self.cache.atomic_write_json(path, payload)
        targets: list[ComparisonTarget] = []
        for row in payload.get("result_set", []):
            identifier = row.get("identifier", "")
            if "_" not in identifier:
                continue
            pdb_id, entity_id = identifier.split("_", 1)
            try:
                meta = self.polymer_entity_metadata(pdb_id, entity_id)
            except Exception:
                meta = {}
            asym_ids = (
                meta.get("rcsb_polymer_entity_container_identifiers", {}).get("asym_ids")
                or meta.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids")
                or []
            )
            auth_ids = meta.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids") or asym_ids
            for idx, chain_id in enumerate(asym_ids):
                targets.append(
                    ComparisonTarget(
                        pdb_id=pdb_id.upper(),
                        chain_id=chain_id,
                        entity_id=entity_id,
                        author_chain_id=auth_ids[idx] if idx < len(auth_ids) else None,
                    )
                )
        return sorted(targets, key=lambda t: (t.pdb_id, t.chain_id, t.entity_id or ""))
