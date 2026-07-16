"""Download facade."""

from __future__ import annotations

from pathlib import Path

from flexres.api.rcsb import RCSBClient


def download_mmcif(client: RCSBClient, pdb_id: str, force: bool = False) -> Path:
    """Download one mmCIF file through the configured RCSB client."""
    return client.download_mmcif(pdb_id, force=force)
