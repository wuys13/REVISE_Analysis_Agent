"""Per-side Moran I workflow; each side builds its own neighbour graph."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from ._shared import analysis_parameters, deterministic_subset, save_table, unavailable


def compute_moran(adata: Any, *, spatial_key: str, n_neighbors: int = 6, transformation_state: str = "untransformed_nonnegative"):
    """Expose the reusable Moran calculation for notebooks."""
    from revise_analysis.methods.spatial import compute_moran as implementation
    return implementation(adata, spatial_key=spatial_key, n_neighbors=n_neighbors, transformation_state=transformation_state)


def run(sample: Any, output_dir: Path, parameters: dict | None = None) -> dict:
    """Save a gene-level Moran table for Raw and SVC independently."""
    supplied = analysis_parameters(parameters, {"moran_n_neighbors", "sample_n_units", "random_state"})
    effective = {
        "moran_n_neighbors": int(supplied.get("moran_n_neighbors", 6)),
        "sample_n_units": supplied.get("sample_n_units"),
        "random_state": int(supplied.get("random_state", 42)),
    }
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs, missing = {}, []
    for side in ("raw", "svc"):
        if reason := sample.expression_unavailable(side):
            unavailable(f"moran:{side}", reason, missing)
            continue
        try:
            adata = deterministic_subset(getattr(sample, side), effective["sample_n_units"], effective["random_state"] + (side == "svc"))
            table = compute_moran(adata, spatial_key=sample.spatial_key, n_neighbors=effective["moran_n_neighbors"], transformation_state=sample.expression(side)["scale"])
            save_table(table, output_dir, f"tables/moran_{side}.csv", outputs)
            record_uncomputed(table, side, missing)
        except (ImportError, KeyError, ValueError) as exc:
            unavailable(f"moran:{side}", str(exc), missing)
    _write_figures(output_dir, outputs, missing)
    status = "partial" if outputs and missing else ("succeeded" if outputs else "skipped")
    return {"status": status, "parameters": effective, "outputs": outputs, "unavailable": missing}


def _write_figures(output_dir: Path, outputs: dict, missing: list[dict]) -> None:
    try:
        from revise_analysis.plotting import plot_moran_distribution
        import pandas as pd
        for side in ("raw", "svc"):
            table = output_dir / f"tables/moran_{side}.csv"
            if table.exists():
                relative = f"figures/moran_{side}.png"
                plot_moran_distribution(pd.read_csv(table), output_dir / relative, title=f"{side.upper()} 原生 Moran I 分布")
                outputs[relative] = relative
    except ImportError as exc:
        unavailable("figures", str(exc), missing)


def record_uncomputed(table, side: str, missing: list[dict]) -> None:
    """Keep per-gene unavailable results visible in the workflow status."""
    if "status" not in table:
        return
    unavailable_rows = table.loc[~table["status"].isin(["computed", "ok"])]
    if not unavailable_rows.empty:
        reasons = unavailable_rows["reason"].value_counts().to_dict() if "reason" in table else {}
        unavailable(f"moran:{side}:genes", f"{len(unavailable_rows)}/{len(table)} genes not computable: {reasons}", missing)
