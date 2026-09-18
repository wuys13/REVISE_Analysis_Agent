"""Moran result figures from saved result tables."""
from __future__ import annotations

from pathlib import Path
import pandas as pd

from .impact import _pyplot, _save


def plot_moran_distribution(moran: pd.DataFrame, destination: Path, *, title: str = "Moran I") -> Path:
    if "moran" not in moran:
        raise ValueError("Moran table requires a moran column")
    values = pd.to_numeric(moran["moran"], errors="coerce").dropna()
    plt = _pyplot()
    figure, axis = plt.subplots()
    if values.empty:
        axis.text(.5, .5, "No computable Moran I values", transform=axis.transAxes,
                  ha="center", va="center", color="#667")
    else:
        axis.hist(values, bins=min(40, max(5, len(values))))
    axis.set_xlabel("Moran I")
    axis.set_title(title)
    return _save(figure, destination)
