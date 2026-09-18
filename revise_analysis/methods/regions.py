"""Spatial window helpers migrated from reconstruction-impact@e82dd13.

These helpers work on a single native coordinate/label side.  They do not
construct paired Raw/SVC objects or require their identifiers to overlap.
"""
from __future__ import annotations

from collections.abc import Sequence
import numpy as np
import pandas as pd


def effective_number(labels: Sequence[object]) -> float:
    """Return Shannon effective diversity, or NaN for no observed labels."""
    values = pd.Series(labels).dropna().astype(str)
    if values.empty:
        return float("nan")
    probabilities = values.value_counts(normalize=True).to_numpy()
    return float(np.exp(-(probabilities * np.log(probabilities)).sum()))


def assign_square_windows(coordinates: pd.DataFrame, *, window_side_length: float, origin: tuple[float, float] | None = None) -> pd.DataFrame:
    """Assign every native coordinate to one non-overlapping square window."""
    if not coordinates.index.is_unique:
        raise ValueError("Coordinates must have unique unit IDs")
    if not np.isfinite(window_side_length) or window_side_length <= 0:
        raise ValueError("window_side_length must be positive and finite")
    if not {"x", "y"}.issubset(coordinates.columns):
        raise KeyError("Coordinates must contain x and y columns")
    work = coordinates.loc[:, ["x", "y"]].copy().astype(float)
    if not np.isfinite(work.to_numpy()).all():
        raise ValueError("Coordinates must be finite")
    x0, y0 = origin if origin is not None else (float(work.x.min()), float(work.y.min()))
    work["window_x_index"] = np.floor((work.x - x0) / window_side_length).astype(int)
    work["window_y_index"] = np.floor((work.y - y0) / window_side_length).astype(int)
    work["window_id"] = work.window_x_index.astype(str) + "_" + work.window_y_index.astype(str)
    work["window_side_length"] = float(window_side_length)
    work["window_center_x"] = x0 + (work["window_x_index"] + .5) * window_side_length
    work["window_center_y"] = y0 + (work["window_y_index"] + .5) * window_side_length
    return work


def _labels_for(labels: pd.Series, ids: pd.Index) -> pd.Series:
    if not labels.index.is_unique or set(labels.index) != set(ids):
        raise ValueError("labels must have unique IDs and exactly match coordinates")
    labels = labels.reindex(ids)
    if labels.isna().any():
        raise ValueError("labels must not contain missing values")
    return labels.astype(str)


def compute_window_diversity(coordinates: pd.DataFrame, labels: pd.Series, *, window_side_length: float, origin: tuple[float, float] | None = None, min_units: int = 4, n_draws: int = 200, random_state: int = 42) -> pd.DataFrame:
    """Rarefy within each window and estimate one side's local diversity.

    Each eligible window contributes draws of size ``min_units``.  The returned
    estimates are therefore comparable across independently processed sides,
    without reusing unit IDs or draws between them.
    """
    if min_units < 1 or n_draws < 1:
        raise ValueError("min_units and n_draws must be at least one")
    assignments = assign_square_windows(coordinates, window_side_length=window_side_length, origin=origin)
    labels = _labels_for(labels, assignments.index)
    codes = pd.Series(pd.factorize(labels, sort=True)[0], index=labels.index)
    generator, rows = np.random.default_rng(random_state), []
    for window_id, frame in assignments.groupby("window_id", sort=True):
        n_units = int(frame.shape[0])
        row = {"window_id": str(window_id), "window_x": float(frame.window_center_x.iloc[0]),
               "window_y": float(frame.window_center_y.iloc[0]), "n_units": n_units,
               "valid_window": n_units >= min_units}
        if n_units < min_units:
            row.update({"k_obs": np.nan, "entropy": np.nan, "neff": np.nan, "evenness": np.nan})
            rows.append(row)
            continue
        values = codes.reindex(frame.index).to_numpy(dtype=int)
        metrics = []
        for _ in range(n_draws):
            draw = values[generator.choice(n_units, size=min_units, replace=False)]
            counts = np.bincount(draw)
            probabilities = counts[counts > 0] / min_units
            entropy = float(-(probabilities * np.log(probabilities)).sum())
            neff = float(np.exp(entropy))
            metrics.append((float(probabilities.size), entropy, neff, float(neff / probabilities.size)))
        means = np.asarray(metrics).mean(axis=0)
        row.update(dict(zip(("k_obs", "entropy", "neff", "evenness"), means, strict=True)))
        rows.append(row)
    return pd.DataFrame(rows)


def select_window_scale(coordinates: pd.DataFrame, *, candidate_window_sides: Sequence[float], min_parent_units: int = 4, origin: tuple[float, float] | None = None) -> tuple[dict, pd.DataFrame]:
    """Select the lower-tie occupancy knee from candidate square side lengths."""
    candidates = sorted({float(side) for side in candidate_window_sides})
    if min_parent_units < 1 or not candidates or any(not np.isfinite(x) or x <= 0 for x in candidates):
        raise ValueError("valid candidate sides and min_parent_units are required")
    rows = []
    for side in candidates:
        occupancy = assign_square_windows(coordinates, window_side_length=side, origin=origin).groupby("window_id").size()
        retained = occupancy[occupancy >= min_parent_units]
        rows.append({"window_side_length": side, "n_tissue_windows": int(occupancy.size), "n_valid_windows": int(retained.size), "valid_window_fraction": float(retained.size / occupancy.size) if occupancy.size else np.nan, "retained_parent_units": int(retained.sum()), "retained_parent_unit_fraction": float(retained.sum() / len(coordinates)) if len(coordinates) else np.nan})
    sensitivity = pd.DataFrame(rows)
    if not (sensitivity.n_valid_windows > 0).any():
        return {"status": "insufficient_support", "window_side_length": None, "min_parent_units": min_parent_units}, sensitivity
    x, y = np.log(sensitivity.window_side_length.to_numpy()), sensitivity.retained_parent_unit_fraction.to_numpy()
    denominator = float(np.hypot(y[-1] - y[0], x[-1] - x[0]))
    distances = np.zeros_like(x) if denominator == 0 else np.abs((y[-1] - y[0]) * x - (x[-1] - x[0]) * y + x[-1] * y[0] - y[-1] * x[0]) / denominator
    sensitivity["knee_distance"] = distances
    maximum = sensitivity.knee_distance.max()
    side = float(sensitivity.loc[np.isclose(sensitivity.knee_distance, maximum), "window_side_length"].min())
    return {"status": "ok", "window_side_length": side, "min_parent_units": min_parent_units}, sensitivity


def select_region_threshold(values: Sequence[float] | np.ndarray, *, n_bootstrap: int = 500, random_state: int = 42) -> tuple[dict, pd.DataFrame]:
    """Use the stable survival-curve breakpoint algorithm from the source method."""
    observed = np.asarray(values, dtype=float)
    observed = observed[np.isfinite(observed)]
    unavailable = lambda n=0: ({"status": "no_stable_threshold", "threshold": None, "n_windows": int(observed.size), "n_valid_bootstrap": n}, pd.DataFrame())
    if observed.size < 200 or n_bootstrap < 1:
        return unavailable()
    def breakpoint(sample):
        unique = np.unique(np.quantile(sample, np.linspace(0, 1, 81)))
        if unique.size < 5:
            return None
        survival = np.array([(sample >= value).mean() for value in unique])
        best = None
        for index in range(max(1, int(np.ceil(unique.size * .1))), min(unique.size - 1, int(np.floor(unique.size * .9))) + 1):
            left_x, left_y, right_x, right_y = unique[:index + 1], np.log(survival[:index + 1]), unique[index:], np.log(survival[index:])
            if left_x.size < 2 or right_x.size < 2:
                continue
            error = float(((left_y - np.polyval(np.polyfit(left_x, left_y, 1), left_x)) ** 2).sum() + ((right_y - np.polyval(np.polyfit(right_x, right_y, 1), right_x)) ** 2).sum())
            best = min(best, (error, float(unique[index]))) if best is not None else (error, float(unique[index]))
        return None if best is None else best[1]
    point = breakpoint(observed)
    if point is None:
        return unavailable()
    generator = np.random.default_rng(random_state)
    valid = np.asarray([x for x in (breakpoint(generator.choice(observed, observed.size, replace=True)) for _ in range(n_bootstrap)) if x is not None], dtype=float)
    audit = pd.DataFrame({"bootstrap_threshold": valid})
    if valid.size < int(np.ceil(.8 * n_bootstrap)):
        return {"status": "no_stable_threshold", "threshold": None, "n_windows": int(observed.size), "n_valid_bootstrap": int(valid.size)}, audit
    lower, upper = np.quantile(valid, [.025, .975])
    stable = np.ptp(observed) > 0 and upper - lower <= .25 * np.ptp(observed)
    return {"status": "ok" if stable else "no_stable_threshold", "threshold": float(np.median(valid)) if stable else None, "point_threshold": float(point), "ci_lower": float(lower), "ci_upper": float(upper), "n_windows": int(observed.size), "n_valid_bootstrap": int(valid.size)}, audit


def assign_anatomy_candidates(window_assignments: pd.DataFrame, level1_labels: pd.Series, *, tumor_label: str = "Tumor", normal_source_label: str = "Intestinal Epithelial") -> pd.DataFrame:
    """Classify every native window as Tumor, Normal, Interface, or Other."""
    if "window_id" not in window_assignments:
        raise KeyError("window_assignments must contain window_id")
    labels = _labels_for(level1_labels, window_assignments.index)
    rows = []
    for window_id, frame in window_assignments.groupby("window_id", sort=True):
        values = labels.reindex(frame.index)
        tumor_units, normal_units = int((values == tumor_label).sum()), int((values == normal_source_label).sum())
        if tumor_units and normal_units:
            region = "Interface"
        elif tumor_units:
            region = "Tumor"
        elif normal_units:
            region = "Normal"
        else:
            region = "Other"
        rows.append({"window_id": str(window_id), "window_x": float(frame.window_center_x.iloc[0]),
                     "window_y": float(frame.window_center_y.iloc[0]), "tumor_units": tumor_units,
                     "normal_units": normal_units, "tumor_candidate": int(bool(tumor_units)),
                     "normal_candidate": int(bool(normal_units)),
                     "union_candidate": int(bool(tumor_units)) + int(bool(normal_units)),
                     "level1_region": region})
    return pd.DataFrame(rows)


def flag_region_windows(window_metrics: pd.DataFrame, *, threshold: float, value_column: str, window_side_length: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply an accepted threshold and summarize its valid-window spatial extent.

    The caller selects ``value_column`` explicitly (for example ``neff`` or
    ``gain_neff``).  This avoids treating State diversity and reconstruction
    Gain as interchangeable scientific quantities.
    """
    required = {"window_id", "valid_window", "n_units", value_column}
    if missing := required - set(window_metrics.columns):
        raise KeyError(f"window_metrics missing columns: {sorted(missing)}")
    if not np.isfinite(threshold) or threshold <= 0 or not np.isfinite(window_side_length) or window_side_length <= 0:
        raise ValueError("threshold and window_side_length must be positive and finite")
    flagged = window_metrics.copy()
    flagged["in_region"] = flagged["valid_window"].astype(bool) & (pd.to_numeric(flagged[value_column], errors="coerce") >= threshold)
    valid, selected = flagged.loc[flagged.valid_window.astype(bool)], flagged.loc[flagged.in_region]
    area = float(window_side_length ** 2)
    valid_units, selected_units = int(valid.n_units.sum()), int(selected.n_units.sum())
    summary = pd.DataFrame([{"value_column": value_column, "threshold": float(threshold), "window_side_length": float(window_side_length), "n_valid_windows": int(valid.shape[0]), "n_region_windows": int(selected.shape[0]), "region_area": float(selected.shape[0] * area), "region_area_fraction": float(selected.shape[0] / valid.shape[0]) if not valid.empty else np.nan, "n_valid_units": valid_units, "n_region_units": selected_units, "region_unit_fraction": float(selected_units / valid_units) if valid_units else np.nan}])
    return summary, flagged
