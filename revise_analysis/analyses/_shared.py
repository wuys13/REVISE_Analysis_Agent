"""Small, explicit helpers shared by analysis workflows."""
from __future__ import annotations

from pathlib import Path
from typing import Any
import hashlib
import json

import numpy as np
import pandas as pd


DEFAULT_SCOPES = ("All", "Fibroblast", "Mono_Macro", "T")
UNIT_ID_HASH_CONTRACT = "sha256-json-utf8-ordered-v1"


def normalize_labels(values: pd.Series) -> pd.Series:
    """Normalize categorical labels without converting missing values to text."""
    return values.astype("string").str.replace("/", "_", regex=False)


def normalize_label(value: Any) -> Any:
    """Normalize one configured label or scope while retaining scalar missingness."""
    if value is None or pd.isna(value):
        return value
    return str(value).replace("/", "_")


def analysis_parameters(parameters: dict | None, allowed: set[str]) -> dict:
    supplied = dict(parameters or {})
    unknown = sorted(set(supplied) - allowed)
    if unknown:
        raise ValueError(f"Unsupported analysis parameters: {', '.join(unknown)}")
    return supplied


def coordinates(adata: Any, spatial_key: str) -> pd.DataFrame:
    try:
        values = adata.obsm[spatial_key]
    except (AttributeError, KeyError) as exc:
        raise KeyError(f"Missing spatial coordinates at obsm[{spatial_key!r}]") from exc
    frame = pd.DataFrame(np.asarray(values), index=pd.Index(adata.obs_names, name="unit_id"))
    if frame.shape[1] < 2:
        raise ValueError("Spatial coordinates require at least two columns")
    return frame.iloc[:, :2].set_axis(["x", "y"], axis=1)


def deterministic_subset(adata: Any, sample_n_units: int | None, random_state: int) -> Any:
    if sample_n_units is None or sample_n_units >= adata.n_obs:
        return adata
    if sample_n_units < 1:
        raise ValueError("sample_n_units must be a positive integer")
    # Sampling each side independently avoids turning an optional pairing into a prerequisite.
    chosen = np.random.default_rng(random_state).choice(adata.n_obs, size=sample_n_units, replace=False)
    return adata[np.sort(chosen)].copy()


def unit_id_digest(values: Any) -> str:
    """Hash an ordered unit-ID cohort with an explicit, portable encoding."""
    payload = json.dumps([str(value) for value in values], ensure_ascii=False,
                         separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def save_table(table: pd.DataFrame | pd.Series, output_dir: Path, relative: str, outputs: dict) -> Path:
    path = output_dir / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(table, pd.Series):
        table.to_frame().to_csv(path)
    else:
        # Method tables such as Moran already carry their identifier as both
        # named index and a visible column; do not emit a second copy.
        table.to_csv(path, index=table.index.name not in table.columns)
    outputs[relative] = relative
    return path


def save_json(value: Any, output_dir: Path, relative: str, outputs: dict) -> Path:
    path = output_dir / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=_json_default), encoding="utf-8")
    outputs[relative] = relative
    return path


def _json_default(value: Any):
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, (Path,)):
        return str(value)
    if isinstance(value, pd.DataFrame):
        return value.to_dict(orient="records")
    if isinstance(value, pd.Series):
        return value.to_dict()
    raise TypeError(f"Not JSON serializable: {type(value).__name__}")


def unavailable(component: str, reason: str, target: list[dict]) -> None:
    target.append({"component": component, "reason": reason})
