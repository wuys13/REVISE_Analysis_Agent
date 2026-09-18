"""Pathway-score figures from saved score series."""
from __future__ import annotations

from pathlib import Path
import pandas as pd

from .impact import _pyplot, _save


def plot_pathway_scores(scores: pd.Series, destination: Path, *, title: str | None = None) -> Path:
    values = pd.to_numeric(scores, errors="coerce").dropna()
    figure, axis = _pyplot().subplots()
    if values.empty:
        axis.text(.5, .5, "No computable AUCell scores", transform=axis.transAxes,
                  ha="center", va="center", color="#667")
    else:
        axis.hist(values, bins=min(40, max(5, len(values))))
    axis.set_xlabel(scores.name or "AUCell score")
    axis.set_title(title or "Pathway activity distribution")
    return _save(figure, destination)
