"""Native-coordinate spatial autocorrelation.

Moran's I calculation follows ``revise.analysis.paired_moran._side_moran``
at reconstruction-impact@e82dd13.  Graph construction is deliberately local
to one AnnData object, rather than implying a Raw/SVC correspondence.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse


def _validate_expression(adata) -> None:
    if adata.n_obs == 0 or adata.n_vars == 0 or adata.X is None:
        raise ValueError("Expression AnnData must have nonempty observations, genes and X")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Observation and gene identifiers must be unique")
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Expression values must be finite and nonnegative")


def _coordinates(adata, spatial_key: str) -> np.ndarray:
    if spatial_key not in adata.obsm:
        raise KeyError(f"Missing spatial coordinates: adata.obsm[{spatial_key!r}]")
    coordinates = np.asarray(adata.obsm[spatial_key], dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[0] != adata.n_obs or coordinates.shape[1] < 2:
        raise ValueError("Spatial coordinates must be n_obs by at least two dimensions")
    coordinates = coordinates[:, :2]
    if not np.isfinite(coordinates).all():
        raise ValueError("Spatial coordinates must be finite")
    return coordinates


def native_knn_graph(coordinates: np.ndarray, *, n_neighbors: int = 6) -> sparse.csr_matrix:
    """Build a deterministic symmetric binary kNN graph with no self edges."""
    if type(n_neighbors) is not int or n_neighbors < 1:
        raise ValueError("n_neighbors must be at least one")
    n = coordinates.shape[0]
    if n < 2:
        return sparse.csr_matrix((n, n), dtype=float)
    k = min(n_neighbors, n - 1)
    try:
        from sklearn.neighbors import NearestNeighbors
    except ImportError as exc:  # pragma: no cover - core project dependency
        raise ImportError("scikit-learn is required to construct the native kNN graph") from exc
    # This avoids materialising an n_obs × n_obs distance matrix. Search
    # using a two-dimensional KD tree bounds memory; the explicit
    # lexicographic ordering settles returned equal-distance ties by unit row.
    model = NearestNeighbors(n_neighbors=k + 1, algorithm="kd_tree", n_jobs=1)
    model.fit(coordinates)
    distances, candidates = model.kneighbors(coordinates, return_distance=True)
    neighbors = np.empty((n, k), dtype=int)
    for row in range(n):
        candidates_row = candidates[row]
        distances_row = distances[row]
        keep = candidates_row != row
        candidates_row, distances_row = candidates_row[keep], distances_row[keep]
        order = np.lexsort((candidates_row, distances_row))
        neighbors[row] = candidates_row[order[:k]]
    rows = np.repeat(np.arange(n), k)
    graph = sparse.csr_matrix((np.ones(rows.size), (rows, neighbors.ravel())), shape=(n, n))
    graph = graph.maximum(graph.T).tocsr()
    graph.setdiag(0)
    graph.eliminate_zeros()
    graph.data[:] = 1.0
    return graph


def compute_moran(
    adata,
    *,
    spatial_key: str = "spatial",
    n_neighbors: int = 6,
    min_units: int = 3,
    block_size: int = 64,
    transformation_state: str = "untransformed_nonnegative",
) -> pd.DataFrame:
    """Compute per-gene Moran's I on a side-native graph.

    Untransformed values are normalized and log1p transformed on a copy.
    Already log1p-transformed values are used as declared, preventing a second
    normalization or log transform.  Unknown scale is rejected for callers to
    record as an unmet method prerequisite.
    """
    if type(n_neighbors) is not int or type(min_units) is not int or type(block_size) is not int or n_neighbors < 1 or min_units < 1 or block_size < 1:
        raise ValueError("min_units and block_size must be at least one")
    aliases = {
        "untransformed": "untransformed_nonnegative",
        "untransformed_nonnegative": "untransformed_nonnegative",
        "log1p": "log1p_nonnegative",
        "log1p_nonnegative": "log1p_nonnegative",
    }
    try:
        state = aliases[str(transformation_state)]
    except KeyError as exc:
        raise ValueError(
            "transformation_state must be untransformed_nonnegative or log1p_nonnegative"
        ) from exc
    _validate_expression(adata)
    coordinates = _coordinates(adata, spatial_key)
    graph = native_knn_graph(coordinates, n_neighbors=n_neighbors)
    reason = "too_few_units" if adata.n_obs < min_units else "no_edges" if graph.nnz == 0 else ""
    work = adata.copy()
    if state == "untransformed_nonnegative":
        try:
            import scanpy as sc
        except ImportError as exc:
            raise ImportError("scanpy is required to normalize expression for Moran I") from exc
        sc.pp.normalize_total(work, target_sum=1e4)
        sc.pp.log1p(work)
    rows, weight_sum = [], float(graph.sum())
    for start in range(0, work.n_vars, block_size):
        block = work.X[:, start:start + block_size]
        values = np.asarray(block.toarray() if sparse.issparse(block) else block, dtype=float)
        if not np.isfinite(values).all():
            raise ValueError("Expression values must be finite")
        centered = values - values.mean(axis=0)
        denominator = np.einsum("ij,ij->j", centered, centered)
        numerator = np.einsum("ij,ij->j", centered, graph @ centered)
        for offset, gene in enumerate(work.var_names[start:start + block_size]):
            if reason:
                value, status, why = np.nan, "insufficient_support", reason
            elif denominator[offset] == 0:
                value, status, why = np.nan, "not_computable", "constant_expression"
            else:
                value = float(work.n_obs * numerator[offset] / (weight_sum * denominator[offset]))
                if not np.isfinite(value):
                    raise ValueError("Moran computation produced a nonfinite value")
                status, why = "computed", ""
            rows.append({"gene_id": str(gene), "moran": value, "status": status, "reason": why, "n_units": int(work.n_obs), "n_edges": int(graph.nnz), "transformation_state": state})
    return pd.DataFrame(rows).set_index("gene_id", drop=False)
