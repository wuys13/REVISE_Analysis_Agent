"""Reusable impact figures from tabular workflow outputs."""
from __future__ import annotations

from pathlib import Path
import pandas as pd


def plot_partition_sizes(summary: pd.DataFrame, destination: Path) -> Path:
    """Plot saved partition summary counts by side and scope."""
    plt = _pyplot()
    if not {"side", "scope"}.issubset(summary):
        raise ValueError("partition summary requires side and scope columns")
    count_column = next((column for column in ("n_clusters", "n_units") if column in summary), None)
    if count_column is None:
        raise ValueError("partition summary needs n_units or n_clusters")
    pivot = summary.pivot(index="scope", columns="side", values=count_column)
    axis = pivot.plot.bar(rot=0)
    axis.set_ylabel(count_column.replace("_", " "))
    axis.set_title("Independent partitions")
    return _save(axis.figure, destination)


def plot_gain(gain: pd.DataFrame, destination: Path) -> Path:
    """Plot only valid windows already selected by the workflow."""
    plt = _pyplot()
    if "gain_neff" not in gain:
        raise ValueError("gain table requires gain_neff")
    values = pd.to_numeric(gain["gain_neff"], errors="coerce").dropna()
    figure, axis = _pyplot().subplots()
    if values.empty:
        _no_values(axis, "No computable valid-window gain values")
    else:
        axis.hist(values, bins=min(30, max(5, len(values))))
    axis.set_xlabel("SVC − Raw effective diversity")
    axis.set_title("Common valid spatial windows")
    return _save(axis.figure, destination)


def plot_anatomy_context(context: pd.DataFrame, destination: Path) -> Path:
    """Draw the saved full-Raw Level1 anatomy context in native coordinates."""
    category = "level1_region" if "level1_region" in context else "anatomy_class"
    required = {"x", "y", category}
    if missing := required - set(context.columns):
        raise ValueError(f"anatomy context missing columns: {sorted(missing)}")
    plt = _pyplot()
    colors = {"Tumor": "#c43c39", "Normal": "#39835b", "Interface": "#e89b20", "Other": "#87929e"}
    figure, axis = plt.subplots(figsize=(6, 5))
    for label, frame in context.groupby(category, dropna=False, sort=True):
        axis.scatter(frame["x"], frame["y"], s=7, alpha=.75,
                     color=colors.get(str(label), "#87929e"), label=str(label), rasterized=True)
    axis.set(title="Raw Level1 anatomy window context", xlabel="x", ylabel="y", aspect="equal")
    axis.legend(title="window context", markerscale=2)
    return _save(figure, destination)


def plot_membership_change(assignments: pd.DataFrame, destination: Path, *, title: str) -> Path:
    """Map the explicit Raw-defined cohort membership comparison only."""
    required = {"x", "y", "unit_changed"}
    if missing := required - set(assignments.columns):
        raise ValueError(f"membership map missing columns: {sorted(missing)}")
    plt = _pyplot()
    figure, axis = plt.subplots(figsize=(6, 5))
    changed = assignments["unit_changed"].astype(bool)
    axis.scatter(assignments.loc[~changed, "x"], assignments.loc[~changed, "y"], s=7,
                 color="#aeb7c2", alpha=.55, label="matched membership", rasterized=True)
    axis.scatter(assignments.loc[changed, "x"], assignments.loc[changed, "y"], s=9,
                 color="#b43c55", alpha=.9, label="changed membership", rasterized=True)
    axis.set(title=title, xlabel="Raw x", ylabel="Raw y", aspect="equal")
    axis.legend(markerscale=2)
    return _save(figure, destination)


def plot_window_field(windows: pd.DataFrame, destination: Path, *, value_column: str, title: str) -> Path:
    """Map a saved local-diversity or positive-gain field and its accepted mask."""
    required = {"window_x", "window_y", value_column}
    if missing := required - set(windows.columns):
        raise ValueError(f"window field missing columns: {sorted(missing)}")
    plt = _pyplot()
    values = pd.to_numeric(windows[value_column], errors="coerce")
    figure, axis = plt.subplots(figsize=(6, 5))
    if values.notna().any():
        scatter = axis.scatter(windows["window_x"], windows["window_y"], c=values, s=35,
                               cmap="viridis", alpha=.9, rasterized=True)
        figure.colorbar(scatter, ax=axis, label=value_column.replace("_", " "))
    else:
        _no_values(axis, f"No computable {value_column} values")
    if "in_region" in windows:
        selected = windows["in_region"].astype(bool)
        axis.scatter(windows.loc[selected, "window_x"], windows.loc[selected, "window_y"],
                     s=72, facecolors="none", edgecolors="#c43c39", linewidths=1.2,
                     label="stable region")
        if selected.any():
            axis.legend(markerscale=1)
    axis.set(title=title, xlabel="window x", ylabel="window y", aspect="equal")
    return _save(figure, destination)


def _no_values(axis, message: str) -> None:
    axis.text(.5, .5, message, transform=axis.transAxes, ha="center", va="center", color="#667")


def _pyplot():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("matplotlib is required to create figures") from exc
    return plt


def _save(figure, destination: Path) -> Path:
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(destination, dpi=160)
    _pyplot().close(figure)
    return destination
