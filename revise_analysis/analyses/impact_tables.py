"""Pure table transformations for the reconstruction-impact workflow.

These helpers deliberately operate on saved-table shaped data.  They make the
scientific joins auditable and keep the workflow from creating a global merged
Raw/SVC expression object.
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd


ANATOMY_CATEGORIES = ("Tumor", "Normal", "Interface", "Other", "Unknown")


def common_valid_delta(raw: pd.DataFrame, svc: pd.DataFrame, *,
                       raw_column: str = "neff", svc_column: str = "neff",
                       delta_column: str = "delta_neff") -> pd.DataFrame:
    """Join native window tables by physical window ID and retain valid pairs."""
    left, right = raw.set_index("window_id"), svc.set_index("window_id")
    common = left.index.intersection(right.index)
    if common.empty:
        return pd.DataFrame(columns=["window_id", "raw_neff", "svc_neff", delta_column])
    valid = (left.loc[common, "valid_window"].astype(bool)
             & right.loc[common, "valid_window"].astype(bool))
    common = common[valid.to_numpy()]
    result = pd.DataFrame({
        "window_id": common.astype(str),
        "raw_neff": pd.to_numeric(left.loc[common, raw_column], errors="coerce").to_numpy(),
        "svc_neff": pd.to_numeric(right.loc[common, svc_column], errors="coerce").to_numpy(),
        "raw_n_units": left.loc[common, "n_units"].to_numpy(),
        "svc_n_units": right.loc[common, "n_units"].to_numpy(),
        "window_x": right.loc[common, "window_x"].to_numpy(),
        "window_y": right.loc[common, "window_y"].to_numpy(),
    })
    result[delta_column] = result["svc_neff"] - result["raw_neff"]
    result["n_units"] = result[["raw_n_units", "svc_n_units"]].min(axis=1).astype(int)
    result["valid_window"] = True
    return result


def label_composition(unit_windows: pd.DataFrame, labels: pd.Series, *,
                      group_columns: Iterable[str]) -> pd.DataFrame:
    """Count reconstructed labels by explicit spatial grouping columns."""
    columns = list(group_columns)
    frame = unit_windows.copy()
    frame["reconstruction_label"] = labels.reindex(frame.index).astype("string")
    frame = frame.dropna(subset=["reconstruction_label"])
    if frame.empty:
        return pd.DataFrame(columns=[*columns, "reconstruction_label", "n_units", "fraction"])
    counts = (frame.groupby([*columns, "reconstruction_label"], dropna=False)
              .size().rename("n_units").reset_index())
    denominators = counts.groupby(columns, dropna=False)["n_units"].transform("sum")
    counts["fraction"] = counts["n_units"] / denominators
    counts["denominator_units"] = denominators
    return counts


def changed_unit_summary(changed: pd.DataFrame, *, group_columns: Iterable[str]) -> pd.DataFrame:
    """Summarize membership changes using only units actually compared."""
    columns = list(group_columns)
    if changed.empty:
        return pd.DataFrame(columns=[*columns, "n_compared", "n_changed", "changed_fraction"])
    frame = changed.dropna(subset=["unit_changed"]).copy()
    rows = []
    for values, group in frame.groupby(columns, dropna=False):
        values = values if isinstance(values, tuple) else (values,)
        flags = group["unit_changed"].astype(bool)
        rows.append({**dict(zip(columns, values, strict=True)),
                     "n_compared": int(flags.size), "n_changed": int(flags.sum()),
                     "changed_fraction": float(flags.mean())})
    return pd.DataFrame(rows)

def aggregate_scores(scores: pd.DataFrame, *, group_columns: Iterable[str]) -> pd.DataFrame:
    """Aggregate already-computed unit scores; missing scores remain missing."""
    columns = list(group_columns)
    required = {*columns, "program", "score"}
    if missing := required - set(scores):
        raise KeyError(f"score table missing columns: {sorted(missing)}")
    rows = []
    for values, group in scores.groupby([*columns, "program"], dropna=False):
        values = values if isinstance(values, tuple) else (values,)
        numeric = pd.to_numeric(group["score"], errors="coerce")
        valid = numeric.dropna()
        rows.append({**dict(zip([*columns, "program"], values, strict=True)),
                     "n_units": int(group.shape[0]), "n_scored": int(valid.size),
                     "score_mean": float(valid.mean()) if not valid.empty else np.nan,
                     "score_median": float(valid.median()) if not valid.empty else np.nan,
                     "missing_status": "ok" if valid.size == group.shape[0] else
                                       ("unavailable" if valid.empty else "partial")})
    return pd.DataFrame(rows)
