"""CSV serialization."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def write_csv(path: Path, rows: list[dict], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    if columns is not None:
        for column in columns:
            if column not in frame:
                frame[column] = None
        frame = frame[columns]
    float_cols = frame.select_dtypes(include=["float"]).columns
    frame[float_cols] = frame[float_cols].round(4)
    frame.to_csv(path, index=False)
