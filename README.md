# FlexRes

FlexRes identifies residues with observed side-chain conformational differences across experimentally determined structures of the same UniProt protein.



## Install

Recommended conda environment, including PyMOL:

```bash
conda env create -f environment.yml
conda activate flexres
```

Pip install with the PyMOL extra:

```bash
python -m pip install "flexres[pymol]"
```

Development install:

```bash
python -m pip install -e ".[test]"
```

PyMOL is the default structural alignment backend. The package also includes an experimental Biopython/SVDSuperimposer backend and a deterministic UniProt-mapped C-alpha Kabsch backend for testing or fallback use in standard Python environments.

## Command-line use

Pocket mode:

```bash
flexres analyse \
  --pdb 8C1D \
  --chain A \
  --ligand-name T5L \
  --ligand-chain A \
  --ligand-residue 402 \
  --pocket-radius 8 \
  --output results/8C1D_A_pocket
```

Whole-protein mode:

```bash
flexres analyse \
  --pdb 8C1D \
  --chain A \
  --all-residues \
  --output results/8C1D_A
```

Use PyMOL for structural superposition:

```bash
flexres analyse \
  --pdb 8C1D \
  --chain A \
  --all-residues \
  --output results/8C1D_A_pymol
```

Use the experimental Biopython backend:

```bash
flexres analyse \
  --pdb 8C1D \
  --chain A \
  --all-residues \
  --alignment-backend biopython \
  --output results/8C1D_A_biopython
```

Exactly one analysis selection is required: `--all-residues`, or the complete ligand specification `--ligand-name`, `--ligand-chain`, and `--ligand-residue`.

## Method

Flexres maps the selected reference chain to one UniProt accession using PDBe/SIFTS-style residue mappings. Candidate structures are searched by the same UniProt accession, not by homologous-sequence expansion. By default, mapped residues must also match the canonical UniProt sequence. Mutated reference or comparison residues are reported and excluded from residue-level side-chain RMSD calculations, while the rest of the structure is still analysed. Residues are compared only when both structures map to the same UniProt accession and UniProt residue number. PDB residue numbers alone are not used as residue identity. PyMOL performs the default structural superposition; UniProt mapping is still used for target selection and residue identity.

In pocket mode, the reference ligand is selected by ligand residue name, ligand chain, and ligand residue number. The reference pocket is every amino-acid residue in the selected reference chain with at least one protein heavy atom within the selected radius of any ligand heavy atom. Hydrogens are ignored. The resulting UniProt residue set is analysed across all comparison chains, even when the comparison structure lacks the ligand.

For each comparison chain and model, the comparison chain is globally aligned to the reference chain. Metrics are then calculated from the aligned comparison coordinates:

- `ca_displacement_angstrom`: Euclidean distance between corresponding C-alpha atoms after global alignment. This is the per-residue CA-only RMSD value.
- `backbone_rmsd_angstrom`: RMSD over N, CA, C, and O without a second residue-level superposition.
- `side_chain_rmsd_angstrom`: RMSD over corresponding side-chain heavy atoms without a second residue-level superposition. This is the primary ranked metric.

Comparison chains with global C-alpha alignment RMSD above `--alignment-rmsd-cutoff` are skipped and reported in `skipped_comparisons.csv`; the default cutoff is 3.0 angstroms.

Missing residues or atoms are not imputed and are not assigned zero. They are reported in `skipped_residues.csv`.

Multiple mapped chains from the same PDB entry are reported separately because they are separate experimental observations. Only coordinate model 1 is analysed for each chain.

## Outputs

The analysis directory contains:

- `residue_comparisons.csv`: ranked residue observations, sorted by descending side-chain RMSD.
- `residue_summary.csv`: final per-residue side-chain RMSD table, with one column per comparison PDB/chain plus mean, standard deviation, maximum, and p90, sorted by descending maximum side-chain RMSD.
- `skipped_comparisons.csv`: skipped chains/models and reasons.
- `skipped_residues.csv`: omitted residue metrics and missing atom details.
- `analysis.json`: full precision structured results.
- `manifest.json`: run timestamp, settings, API endpoints, counts, cache path, and dependency versions.
- `flexres.log`: validation, retrieval, alignment, skip, and output log messages.

No output labels residues as flexible, inflexible, recommended, or not recommended.

## Development

```bash
python -m pytest -q
```

Live API tests should be marked separately and are not required for the default suite:

```bash
python -m pytest -q -m live
```

## APIs and cache

Flexres uses RCSB for mmCIF downloads and UniProt-to-PDB search, and PDBe/SIFTS-style mappings for chain and residue mapping. Cached files are stored under the platform cache directory by default, normally `~/.cache/flexres/`, with separate folders for structures, mappings, searches, and metadata. Use `--cache-dir` to override and `--force-redownload` to refresh cached content.

