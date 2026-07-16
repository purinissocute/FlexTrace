"""Command-line interface."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import typer
from rich.console import Console

from flexres.exceptions import FlexresError
from flexres.models import AnalysisSettings
from flexres.pipeline import run_analysis

app = typer.Typer(help="Identify residues with experimentally observed side-chain conformational variation.")
console = Console()


@app.callback()
def main() -> None:
    """Flexres command group."""


@app.command()
def analyse(
    pdb: str = typer.Option(..., "--pdb", help="Reference PDB ID."),
    chain: str = typer.Option(..., "--chain", help="Reference protein chain ID."),
    output: Path = typer.Option(..., "--output", help="Output directory."),
    all_residues: bool = typer.Option(False, "--all-residues", help="Analyse all comparable residues."),
    ligand_name: str | None = typer.Option(None, "--ligand-name", help="Reference ligand residue name."),
    ligand_chain: str | None = typer.Option(None, "--ligand-chain", help="Reference ligand chain ID."),
    ligand_residue: int | None = typer.Option(None, "--ligand-residue", help="Reference ligand residue number."),
    ligand_insertion_code: str = typer.Option("", "--ligand-insertion-code", help="Reference ligand insertion code."),
    pocket_radius: float = typer.Option(8.0, "--pocket-radius", help="Pocket heavy-atom radius in angstroms."),
    include_reference_entry: bool = typer.Option(True, "--include-reference-entry/--exclude-reference-entry"),
    cache_dir: Path | None = typer.Option(None, "--cache-dir", help="Cache directory."),
    force_redownload: bool = typer.Option(False, "--force-redownload", help="Ignore cached API responses and mmCIF files."),
    max_workers: int = typer.Option(4, "--max-workers", min=1),
    alignment_backend: Literal["pymol", "kabsch", "biopython"] = typer.Option("pymol", "--alignment-backend", help="Structural alignment backend."),
    pymol_executable: str = typer.Option("pymol", "--pymol-executable", help="PyMOL executable path when using --alignment-backend pymol."),
    alignment_rmsd_cutoff: float = typer.Option(3.0, "--alignment-rmsd-cutoff", help="Skip comparison chains with global C-alpha alignment RMSD above this value."),
    output_format: str = typer.Option("csv,json", "--output-format", help="Comma-separated output formats: csv,json."),
    exclude_mutations: bool = typer.Option(True, "--exclude-mutations/--include-mutations", help="Require mapped residues to match the canonical UniProt sequence."),
    log_level: str = typer.Option("INFO", "--log-level"),
    resume: bool = typer.Option(False, "--resume", help="Reserved for resumable runs."),
) -> None:
    """Run residue variability analysis."""
    ligand_complete = ligand_name is not None and ligand_chain is not None and ligand_residue is not None
    if all_residues == ligand_complete:
        raise typer.BadParameter("choose exactly one of --all-residues or complete ligand options")
    formats = {item.strip().lower() for item in output_format.split(",") if item.strip()}
    unknown = formats - {"csv", "json"}
    if unknown:
        raise typer.BadParameter(f"unsupported output format(s): {', '.join(sorted(unknown))}")
    settings = AnalysisSettings(
        pdb_id=pdb,
        chain_id=chain,
        output_dir=output,
        all_residues=all_residues,
        ligand_name=ligand_name,
        ligand_chain=ligand_chain,
        ligand_residue=ligand_residue,
        ligand_insertion_code=ligand_insertion_code,
        pocket_radius=pocket_radius,
        include_reference_entry=include_reference_entry,
        cache_dir=cache_dir,
        force_redownload=force_redownload,
        max_workers=max_workers,
        alignment_backend=alignment_backend,
        pymol_executable=pymol_executable,
        alignment_rmsd_cutoff=alignment_rmsd_cutoff,
        output_formats=formats,
        exclude_mutations=exclude_mutations,
        resume=resume,
        log_level=log_level,
    )
    try:
        run_analysis(settings)
    except FlexresError as exc:
        console.print(f"[red]Error:[/red] {exc}")
        raise typer.Exit(1) from exc
    except Exception as exc:  # noqa: BLE001
        console.print(f"[red]Unexpected error:[/red] {exc}")
        raise typer.Exit(1) from exc
    console.print(f"Wrote flexres outputs to {output}")


if __name__ == "__main__":
    app()
