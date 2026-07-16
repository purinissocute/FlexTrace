"""Run manifest generation."""

from __future__ import annotations

import importlib.metadata
import platform
from datetime import datetime, timezone
from typing import Any

from flexres import __version__
from flexres.config import PDBe_API, RCSB_DATA_API, RCSB_FILES, RCSB_SEARCH_API
from flexres.models import AnalysisSettings


def dependency_versions() -> dict[str, str]:
    names = ["gemmi", "numpy", "pandas", "scipy", "requests", "pydantic", "typer", "rich", "tenacity", "platformdirs"]
    versions = {"python": platform.python_version()}
    for name in names:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "not installed"
    return versions


def build_manifest(
    settings: AnalysisSettings,
    uniprot_id: str,
    counts: dict[str, int],
) -> dict[str, Any]:
    return {
        "package_version": __version__,
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "command_line_arguments": settings.model_dump(mode="json"),
        "reference_pdb_id": settings.pdb_id,
        "reference_chain": settings.chain_id,
        "uniprot_accession": uniprot_id,
        "api_endpoints_used": {
            "rcsb_data": RCSB_DATA_API,
            "rcsb_search": RCSB_SEARCH_API,
            "rcsb_files": RCSB_FILES,
            "pdbe": PDBe_API,
        },
        **counts,
        "cache_directory": str(settings.cache_dir) if settings.cache_dir else None,
        "software_dependency_versions": dependency_versions(),
    }
