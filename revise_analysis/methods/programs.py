"""Pathway scoring through the audited OmicVerse AUCell provider.

This is migrated from reconstruction-impact@e82dd13.  There is no substitute
rank-score implementation: provider availability is part of the result's
scientific provenance.
"""
from __future__ import annotations

from collections.abc import Sequence
import importlib
from importlib import metadata

import numpy as np
import pandas as pd
from scipy import sparse

from .errors import PrerequisiteUnavailable


def _require_omicverse():
    try:
        return importlib.import_module("omicverse")
    except (ImportError, ModuleNotFoundError) as exc:
        if isinstance(exc, ModuleNotFoundError) and exc.name != "omicverse":
            raise ImportError(f"OmicVerse import failed because dependency {exc.name!r} is unavailable") from exc
        raise ImportError("OmicVerse is required for AUCell scoring; install the pathway optional dependency.") from exc


def _provider_metadata(omicverse) -> dict[str, str]:
    provider = getattr(getattr(omicverse, "single", None), "geneset_aucell", None)
    if provider is None or not callable(provider):
        raise RuntimeError("OmicVerse single.geneset_aucell provider is unavailable")
    version = getattr(omicverse, "__version__", None)
    if not version:
        try:
            version = metadata.version("omicverse")
        except metadata.PackageNotFoundError:
            version = "unknown"
    return {"provider": "omicverse", "provider_version": str(version), "scorer": "single.geneset_aucell", "provider_module": str(getattr(provider, "__module__", "unknown"))}


def compute_pathway_scores(adata, genes: Sequence[str], *, score_name: str = "EMT", auc_threshold: float = 0.01, seed: int = 42) -> dict:
    """Score one gene set with OmicVerse AUCell on a copy of ``adata``.

    Returns score values indexed by native observation IDs, coverage and exact
    provider metadata.  Missing optional dependencies fail explicitly.
    """
    if adata.n_obs == 0 or adata.n_vars == 0 or adata.X is None or not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise PrerequisiteUnavailable("Expression AnnData must have nonempty unique axes and X")
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if not np.isfinite(values).all() or (values < 0).any():
        raise PrerequisiteUnavailable("Expression values must be finite and nonnegative")
    if not score_name or not np.isfinite(auc_threshold) or not 0 < auc_threshold <= 1 or type(seed) is not int:
        raise PrerequisiteUnavailable("score_name, auc_threshold and seed must be valid")
    requested = list(dict.fromkeys(str(gene) for gene in genes))
    if not requested:
        raise PrerequisiteUnavailable("genes must be nonempty")
    overlapping = [gene for gene in requested if gene in adata.var_names]
    coverage = {"requested_genes": requested, "n_requested": len(requested), "present_genes": overlapping, "n_present": len(overlapping), "n_missing": len(requested) - len(overlapping), "fraction_present": len(overlapping) / len(requested)}
    if not overlapping:
        raise PrerequisiteUnavailable("Gene set has no overlap with adata.var_names")
    omicverse = _require_omicverse()
    provider = _provider_metadata(omicverse)
    # OmicVerse 1.7.5 interprets this parameter as a detection quantile,
    # not as the final fraction of ranked genes. Match its declared choices.
    if auc_threshold not in (0.01, 0.05, 0.10, 0.50, 1.0):
        raise PrerequisiteUnavailable("OmicVerse AUC_threshold must be one of 0.01, 0.05, 0.1, 0.5, 1.0")
    detected = np.asarray((adata.X != 0).sum(axis=1)).ravel() if sparse.issparse(adata.X) else np.count_nonzero(adata.X, axis=1)
    effective_fraction = float(pd.Series(detected).quantile(auc_threshold) / adata.n_vars)
    rank_cutoff = int(round(effective_fraction * adata.n_vars))
    if rank_cutoff <= 1 or rank_cutoff >= adata.n_vars:
        raise PrerequisiteUnavailable("AUCell unavailable: detection-derived rank cutoff must be between 2 and n_genes-1; input is too sparse or uniformly detected")
    work = adata.copy()
    if not sparse.issparse(work.X):
        work.X = sparse.csr_matrix(work.X)
    omicverse.single.geneset_aucell(adata=work, geneset_name=score_name, geneset=overlapping, AUC_threshold=float(auc_threshold), seed=int(seed))
    score_key = f"{score_name}_aucell"
    if score_key not in work.obs:
        raise RuntimeError(f"AUCell did not create {score_key!r}")
    scores = pd.Series(work.obs[score_key].to_numpy(dtype=float), index=work.obs_names, name=score_name)
    if not np.isfinite(scores.to_numpy()).all():
        raise RuntimeError(f"AUCell returned non-finite scores in {score_key!r}")
    return {"scores": scores, "coverage": coverage, "metadata": {**provider, "score_name": score_name, "score_key": score_key, "auc_threshold": float(auc_threshold), "effective_rank_fraction": effective_fraction, "rank_cutoff": rank_cutoff, "seed": int(seed)}}
