"""Residue-level summary tables."""

from __future__ import annotations

import pandas as pd


def summarize_residue_comparisons(rows: list[dict]) -> list[dict]:
    if not rows:
        return []
    frame = pd.DataFrame(rows)
    frame["comparison_label"] = frame.apply(_comparison_label, axis=1)
    group_cols = [
        "uniprot_id",
        "uniprot_residue_number",
        "reference_residue_name",
        "reference_residue_number",
    ]
    grouped = frame.groupby(group_cols, dropna=False)["side_chain_rmsd_angstrom"]
    summary = grouped.agg(
        comparison_count="count",
        mean_side_chain_rmsd_angstrom="mean",
        sd_side_chain_rmsd_angstrom="std",
        max_side_chain_rmsd_angstrom="max",
        p90_side_chain_rmsd_angstrom=lambda s: s.quantile(0.9),
    ).reset_index()
    wide = frame.pivot_table(
        index=group_cols,
        columns="comparison_label",
        values="side_chain_rmsd_angstrom",
        aggfunc="mean",
        dropna=False,
    ).reset_index()
    wide.columns.name = None
    summary = summary.merge(wide, on=group_cols, how="left")
    comparison_cols = sorted(
        column for column in summary.columns if column not in {
            *group_cols,
            "comparison_count",
            "mean_side_chain_rmsd_angstrom",
            "sd_side_chain_rmsd_angstrom",
            "max_side_chain_rmsd_angstrom",
            "p90_side_chain_rmsd_angstrom",
        }
    )
    summary = summary[
        [
            *group_cols,
            "comparison_count",
            *comparison_cols,
            "mean_side_chain_rmsd_angstrom",
            "sd_side_chain_rmsd_angstrom",
            "max_side_chain_rmsd_angstrom",
            "p90_side_chain_rmsd_angstrom",
        ]
    ]
    summary = summary.sort_values(
        ["max_side_chain_rmsd_angstrom", "uniprot_residue_number"],
        ascending=[False, True],
        kind="mergesort",
    )
    return summary.where(pd.notnull(summary), None).to_dict(orient="records")


def _comparison_label(row: pd.Series) -> str:
    pdb_id = row["comparison_pdb_id"]
    chain = row["comparison_chain"]
    return f"{pdb_id}_{chain}_side_chain_rmsd_angstrom"
