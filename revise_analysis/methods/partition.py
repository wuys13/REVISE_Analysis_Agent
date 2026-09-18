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
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score, f1_score, mutual_info_score, normalized_mutual_info_score


@dataclass(frozen=True)
class PartitionComparison:
    summary: pd.DataFrame
    mapping: pd.DataFrame
    contingency: pd.DataFrame
    assignments: pd.DataFrame


def _validate_expression(adata) -> None:
    if adata.n_obs == 0 or adata.n_vars == 0 or adata.X is None:
        raise ValueError("Expression AnnData must have nonempty observations, genes and X")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Observation and gene identifiers must be unique")
    values = adata.X.data if issparse(adata.X) else np.asarray(adata.X)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Expression values must be finite and nonnegative")


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


def compute_partition(adata, *, resolution: float = 0.5, n_top_genes: int = 2000, random_state: int = 42, min_genes: int = 0, min_cells: int = 0) -> dict:
    """Run an independent Leiden partition without altering the input AnnData."""
    if (not np.isfinite(resolution) or resolution <= 0 or type(n_top_genes) is not int or n_top_genes < 1
            or type(random_state) is not int or type(min_genes) is not int or min_genes < 0
            or type(min_cells) is not int or min_cells < 0):
        raise ValueError("resolution/n_top_genes/min_genes/min_cells are invalid")
    _validate_expression(adata)
    if adata.n_obs < 3 or adata.n_vars < 2:
        raise ValueError("Leiden partitioning requires at least three units and two genes")
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
        raise ValueError("QC retained fewer than three units or two genes")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Some cells have zero counts")
        sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    if work.n_vars > n_top_genes:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)
            sc.pp.highly_variable_genes(work, n_top_genes=n_top_genes, flavor="seurat_v3", check_values=False)
        work = work[:, work.var["highly_variable"].to_numpy()].copy()
    n_pcs = min(40, work.n_obs - 1, work.n_vars - 1)
    if n_pcs < 1:
        raise ValueError("No usable principal components after filtering")
    sc.tl.pca(work, n_comps=n_pcs, svd_solver="arpack")
    sc.pp.neighbors(work, n_neighbors=min(10, work.n_obs - 1), n_pcs=n_pcs, random_state=random_state)
    sc.tl.leiden(work, resolution=float(resolution), key_added="partition", random_state=random_state, flavor="igraph", n_iterations=2, directed=False)
    labels = work.obs["partition"].astype(str).rename("partition")
    embedding = None
    if "X_pca" in work.obsm:
        values = np.asarray(work.obsm["X_pca"])
        n_dimensions = min(2, values.shape[1])
        embedding = pd.DataFrame(values[:, :n_dimensions], index=work.obs_names, columns=[f"PC{i + 1}" for i in range(n_dimensions)])
    return {"labels": labels, "summary": {"status": "ok", "resolution": float(resolution), "n_units": int(work.n_obs), "n_clusters": int(labels.nunique()), "input_units": int(input_units), "input_genes": int(input_genes), "n_genes": int(work.n_vars), "random_state": int(random_state)}, "embedding": embedding}
