"""UniProt API helpers."""

from __future__ import annotations

import requests

from flexres.api.retry import retry, stop_after_attempt, wait_exponential
from flexres.cache.manager import CacheManager
from flexres.config import UNIPROT_API, USER_AGENT


def normalize_uniprot_accession(accession: str) -> str:
    """Normalize a UniProt accession for file names and comparisons."""
    return accession.strip().upper()


class UniProtClient:
    """Small cache-aware UniProt client."""

    def __init__(self, cache: CacheManager, timeout: float = 30.0) -> None:
        self.cache = cache
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": USER_AGENT})

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=8))
    def _get_text(self, url: str) -> str:
        response = self.session.get(url, timeout=self.timeout)
        response.raise_for_status()
        return response.text

    def canonical_sequence(self, accession: str, force: bool = False) -> str:
        accession = normalize_uniprot_accession(accession)
        path = self.cache.uniprot_sequence_path(accession)
        if not force and self.cache.is_valid_file(path):
            return _parse_fasta(path.read_text(encoding="utf-8"))
        payload = self._get_text(f"{UNIPROT_API}/{accession}.fasta")
        sequence = _parse_fasta(payload)
        self.cache.atomic_write_bytes(path, payload.encode("utf-8"))
        return sequence


def _parse_fasta(payload: str) -> str:
    lines = [line.strip() for line in payload.splitlines() if line.strip() and not line.startswith(">")]
    return "".join(lines).upper()
