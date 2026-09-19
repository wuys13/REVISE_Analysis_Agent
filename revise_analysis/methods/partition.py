"""Independent Leiden partitions and the paired membership comparison.

The comparison implementation is migrated from
``revise.analysis.basic.partition_change`` at reconstruction-impact@e82dd13.
It is intentionally only used after matching observation identifiers.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Hashable, Sequence
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.sparse import issparse, SparseEfficiencyWarning
from .errors import PrerequisiteUnavailable
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score, f1_score, mutual_info_score, normalized_mutual_info_score


@dataclass(frozen=True)
class PartitionComparison:
    summary: pd.DataFrame
    mapping: pd.DataFrame
    contingency: pd.DataFrame
    assignments: pd.DataFrame


@dataclass(frozen=True)
class PreparedPartition:
    """A fixed cohort, feature set, PCA and neighbour graph for Leiden sweeps."""

    work: object
    embedding: pd.DataFrame | None
    summary: dict
    random_state: int
    transformation_state: str


def _validate_expression(adata) -> None:
    if adata.n_obs == 0 or adata.n_vars == 0 or adata.X is None:
        raise PrerequisiteUnavailable("Expression AnnData must have nonempty observations, genes and X")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise PrerequisiteUnavailable("Observation and gene identifiers must be unique")
    values = adata.X.data if issparse(adata.X) else np.asarray(adata.X)
    if not np.isfinite(values).all() or (values < 0).any():
        raise PrerequisiteUnavailable("Expression values must be finite and nonnegative")


def _as_labels(values: pd.Series | Sequence[Hashable], name: str) -> pd.Series:
    series = values.copy() if isinstance(values, pd.Series) else pd.Series(values)
    if not series.index.is_unique:
        raise ValueError(f"{name} labels must have unique unit IDs")
    if series.empty or series.isna().any():
        raise ValueError(f"{name} labels must be nonempty and contain no missing values")
    return series.astype(str)


def _label_key(value: str) -> tuple[int, float | str, str]:
    try:
        return (0, float(value), value)
    except ValueError:
        return (1, value, value)


def _vi(left: pd.Series, right: pd.Series) -> float:
    entropy = lambda s: float(-(s.value_counts(normalize=True).pipe(lambda p: p * np.log(p))).sum())
    return float(entropy(left) + entropy(right) - 2 * mutual_info_score(left, right))


def compare_partitions(raw_labels: pd.Series | Sequence[Hashable], reconstructed_labels: pd.Series | Sequence[Hashable], *, comparison_edge: str = "raw_to_recon_expression") -> PartitionComparison:
    """Hungarian comparison for two already paired label vectors."""
    raw, recon = _as_labels(raw_labels, "Raw"), _as_labels(reconstructed_labels, "Reconstructed")
    if set(raw.index) != set(recon.index):
        raise ValueError("Raw and reconstructed labels must contain the same unit IDs")
    recon = recon.reindex(raw.index)
    raw_categories, recon_categories = sorted(raw.unique(), key=_label_key), sorted(recon.unique(), key=_label_key)
    contingency = pd.crosstab(raw, recon).reindex(index=raw_categories, columns=recon_categories, fill_value=0)
    rows, cols = linear_sum_assignment(-contingency.to_numpy())
    pairs = {raw_categories[i]: recon_categories[j] for i, j in zip(rows, cols)}
    matched_recon = set(pairs.values())
    mapping_rows = []
    for raw_cluster in raw_categories:
        recon_cluster = pairs.get(raw_cluster)
        raw_n = int(contingency.loc[raw_cluster].sum())
        overlap = int(contingency.loc[raw_cluster, recon_cluster]) if recon_cluster is not None else 0
        recon_n = int(contingency[recon_cluster].sum()) if recon_cluster is not None else 0
        precision, recall = (overlap / recon_n if recon_n else np.nan), (overlap / raw_n if raw_n else np.nan)
        mapping_rows.append({"raw_cluster": raw_cluster, "recon_cluster": recon_cluster, "matched": recon_cluster is not None, "overlap_n": overlap, "raw_cluster_n": raw_n, "recon_cluster_n": recon_n, "precision": precision, "recall": recall, "f1": 2 * precision * recall / (precision + recall) if precision + recall else 0.0})
    for recon_cluster in recon_categories:
        if recon_cluster not in matched_recon:
            mapping_rows.append({"raw_cluster": pd.NA, "recon_cluster": recon_cluster, "matched": False, "overlap_n": 0, "raw_cluster_n": 0, "recon_cluster_n": int(contingency[recon_cluster].sum()), "precision": 0.0, "recall": np.nan, "f1": 0.0})
    mapping = pd.DataFrame(mapping_rows)
    mapped_raw = raw.map(pairs).fillna("__unmatched_raw_cluster__")
    changed = mapped_raw != recon
    assignments = pd.DataFrame({"unit_id": raw.index.astype(str), "raw_cluster": raw.to_numpy(), "matched_raw_cluster": mapped_raw.to_numpy(), "recon_cluster": recon.to_numpy(), "unit_changed": changed.to_numpy(), "comparison_edge": comparison_edge}, index=raw.index)
    labels_for_f1 = sorted(set(mapped_raw) | set(recon))
    accuracy = float(1 - changed.mean())
    macro_f1 = float(f1_score(recon, mapped_raw, labels=labels_for_f1, average="macro", zero_division=0))
    summary = pd.DataFrame([{"comparison_edge": comparison_edge, "n_units": int(raw.size), "n_raw_clusters": len(raw_categories), "n_recon_clusters": len(recon_categories), "matched_accuracy": accuracy, "st_unit_change_fraction": float(1 - accuracy), "matched_macro_f1": macro_f1, "balanced_cluster_change": float(1 - macro_f1), "ARI": float(adjusted_rand_score(raw, recon)), "AMI": float(adjusted_mutual_info_score(raw, recon)), "NMI": float(normalized_mutual_info_score(raw, recon)), "VI": _vi(raw, recon), "status": "ok"}])
    return PartitionComparison(summary, mapping, contingency, assignments)


def compare_membership(raw_labels: pd.Series | Sequence[Hashable], svc_labels: pd.Series | Sequence[Hashable]) -> PartitionComparison:
    """Compare membership only on shared identifiers; each side otherwise stays native."""
    raw, svc = _as_labels(raw_labels, "Raw"), _as_labels(svc_labels, "SVC")
    common = raw.index[raw.index.isin(svc.index)]
    if common.empty:
        summary = pd.DataFrame([{"comparison_edge": "raw_to_svc_membership", "n_units": 0, "status": "no_common_units"}])
        mapping = pd.DataFrame(columns=["raw_cluster", "recon_cluster", "matched", "overlap_n", "raw_cluster_n", "recon_cluster_n", "precision", "recall", "f1"])
        contingency = pd.DataFrame()
        assignments = pd.DataFrame(columns=["unit_id", "raw_cluster", "matched_raw_cluster", "recon_cluster", "unit_changed", "comparison_edge"])
        return PartitionComparison(summary, mapping, contingency, assignments)
    return compare_partitions(raw.loc[common], svc.reindex(common), comparison_edge="raw_to_svc_membership")


_TRANSFORMATION_ALIASES = {
    "untransformed": "untransformed_nonnegative",
    "untransformed_nonnegative": "untransformed_nonnegative",
    "log1p": "log1p_nonnegative",
    "log1p_nonnegative": "log1p_nonnegative",
}


def _transformation_state(value: str) -> str:
    try:
        return _TRANSFORMATION_ALIASES[str(value)]
    except KeyError as exc:
        raise ValueError(
            "transformation_state must be untransformed_nonnegative or log1p_nonnegative"
        ) from exc


def prepare_partition_graph(
    adata,
    *,
    n_top_genes: int = 2000,
    random_state: int = 42,
    min_genes: int = 0,
    min_cells: int = 0,
    transformation_state: str = "untransformed_nonnegative",
) -> PreparedPartition:
    """Prepare one immutable cohort/feature/PCA/kNN basis for Leiden runs.

    Untransformed nonnegative values are normalized and logged on their full
    gene axis before the dispersion-based ``seurat`` HVG method is applied.
    This state does not assert integer counts, so it must not imply the
    count-model ``seurat_v3`` method. Already log1p-transformed values use the
    same HVG method without a second transformation.
    """
    if (type(n_top_genes) is not int or n_top_genes < 1 or type(random_state) is not int
            or type(min_genes) is not int or min_genes < 0
            or type(min_cells) is not int or min_cells < 0):
        raise ValueError("n_top_genes/min_genes/min_cells/random_state are invalid")
    state = _transformation_state(transformation_state)
    _validate_expression(adata)
    if adata.n_obs < 3 or adata.n_vars < 2:
        raise PrerequisiteUnavailable("Leiden partitioning requires at least three units and two genes")
    try:
        import scanpy as sc
    except ImportError as exc:
        raise ImportError("scanpy is required for Leiden partitioning") from exc
    work = adata.copy()
    input_units, input_genes = work.n_obs, work.n_vars
    if min_genes:
        sc.pp.filter_cells(work, min_genes=min_genes)
    if min_cells:
        sc.pp.filter_genes(work, min_cells=min_cells)
    if work.n_obs < 3 or work.n_vars < 2:
        raise PrerequisiteUnavailable("QC retained fewer than three units or two genes")
    preprocessing = "declared_log1p_no_transform"
    if state == "untransformed_nonnegative":
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Some cells have zero counts")
            sc.pp.normalize_total(work, target_sum=1e4)
        sc.pp.log1p(work)
        preprocessing = "normalize_total_1e4_log1p"
    hvg_flavor = "seurat"
    if work.n_vars > n_top_genes:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)
            sc.pp.highly_variable_genes(work, n_top_genes=n_top_genes, flavor=hvg_flavor)
        if not work.var["highly_variable"].astype(bool).any():
            raise PrerequisiteUnavailable("HVG selection retained no genes")
        work = work[:, work.var["highly_variable"].to_numpy()].copy()
    n_pcs = min(40, work.n_obs - 1, work.n_vars - 1)
    if n_pcs < 1:
        raise PrerequisiteUnavailable("No usable principal components after filtering")
    sc.tl.pca(work, n_comps=n_pcs, svd_solver="arpack", random_state=random_state)
    sc.pp.neighbors(work, n_neighbors=min(10, work.n_obs - 1), n_pcs=n_pcs, random_state=random_state)
    embedding = None
    if "X_pca" in work.obsm:
        values = np.asarray(work.obsm["X_pca"])
        n_dimensions = min(2, values.shape[1])
        embedding = pd.DataFrame(values[:, :n_dimensions], index=work.obs_names, columns=[f"PC{i + 1}" for i in range(n_dimensions)])
    summary = {
        "status": "prepared", "n_units": int(work.n_obs), "input_units": int(input_units),
        "input_genes": int(input_genes), "n_genes": int(work.n_vars),
        "random_state": int(random_state), "transformation_state": state,
        "hvg_flavor": hvg_flavor, "preprocessing": preprocessing,
    }
    return PreparedPartition(work, embedding, summary, random_state, state)


def partition_prepared_graph(prepared: PreparedPartition, *, resolution: float = 0.5) -> dict:
    """Run Leiden on a copy of a prepared graph without changing the basis."""
    if not isinstance(prepared, PreparedPartition):
        raise TypeError("prepared must be a PreparedPartition")
    if not np.isfinite(resolution) or resolution <= 0:
        raise ValueError("resolution must be positive and finite")
    try:
        import scanpy as sc
    except ImportError as exc:
        raise ImportError("scanpy is required for Leiden partitioning") from exc
    work = prepared.work.copy()
    sc.tl.leiden(
        work, resolution=float(resolution), key_added="partition",
        random_state=prepared.random_state, flavor="igraph", n_iterations=2, directed=False,
    )
    labels = work.obs["partition"].astype(str).rename("partition")
    summary = dict(prepared.summary)
    summary.update({"status": "ok", "resolution": float(resolution), "n_clusters": int(labels.nunique())})
    embedding = None if prepared.embedding is None else prepared.embedding.copy()
    return {"labels": labels, "summary": summary, "embedding": embedding}


def select_matched_k_partition(
    prepared: PreparedPartition,
    *,
    candidate_resolutions: Sequence[float],
    target_k: int,
    main_resolution: float = 0.5,
) -> dict:
    """Select an exact or nearest-K Raw sensitivity partition deterministically."""
    if type(target_k) is not int or target_k < 1:
        raise ValueError("target_k must be a positive integer")
    if not np.isfinite(main_resolution) or main_resolution <= 0:
        raise ValueError("main_resolution must be positive and finite")
    resolutions = sorted({float(value) for value in candidate_resolutions})
    if not resolutions or any(not np.isfinite(value) or value <= 0 for value in resolutions):
        raise ValueError("candidate_resolutions must contain positive finite values")
    results, rows = {}, []
    for resolution in resolutions:
        result = partition_prepared_graph(prepared, resolution=resolution)
        actual_k = int(result["summary"]["n_clusters"])
        results[resolution] = result
        rows.append({
            "resolution": resolution, "target_k": target_k, "actual_k": actual_k,
            "k_delta": actual_k - target_k, "absolute_k_delta": abs(actual_k - target_k),
            "exact": actual_k == target_k,
        })
    candidates = pd.DataFrame(rows)
    exact = candidates.loc[candidates.exact]
    pool = exact if not exact.empty else candidates.loc[
        candidates.absolute_k_delta == candidates.absolute_k_delta.min()
    ]
    selected_index = min(
        pool.index,
        key=lambda index: (abs(float(pool.loc[index, "resolution"]) - main_resolution),
                           float(pool.loc[index, "resolution"])),
    )
    candidates["selected"] = False
    candidates.loc[selected_index, "selected"] = True
    selected_row = candidates.loc[selected_index]
    status = "exact" if bool(selected_row.exact) else "nearest"
    selection = {
        "status": status, "target_k": target_k,
        "actual_k": int(selected_row.actual_k), "k_delta": int(selected_row.k_delta),
        "resolution": float(selected_row.resolution), "main_resolution": float(main_resolution),
    }
    return {
        "selection": selection,
        "candidates": candidates,
        "selected": results[float(selected_row.resolution)],
    }


def compute_partition(
    adata,
    *,
    resolution: float = 0.5,
    n_top_genes: int = 2000,
    random_state: int = 42,
    min_genes: int = 0,
    min_cells: int = 0,
    transformation_state: str = "untransformed_nonnegative",
) -> dict:
    """Run an independent Leiden partition without altering the input AnnData."""
    if not np.isfinite(resolution) or resolution <= 0:
        raise ValueError("resolution must be positive and finite")
    prepared = prepare_partition_graph(
        adata, n_top_genes=n_top_genes, random_state=random_state,
        min_genes=min_genes, min_cells=min_cells,
        transformation_state=transformation_state,
    )
    return partition_prepared_graph(prepared, resolution=resolution)
