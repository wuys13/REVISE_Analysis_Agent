"""Audited memory adapter for the pinned OmicVerse 1.7.5 AUCell ranker."""
from __future__ import annotations

import hashlib
import importlib
import inspect
from importlib import metadata
from contextlib import contextmanager
from functools import partial
import threading
from typing import Any, Callable

import numpy as np


SOURCE_HASH_CONTRACT = "sha256-inspect-getsource-utf8-v1"
PERMUTATION_HASH_CONTRACT = "sha256-int64-little-endian-c-order-v1"
ADAPTER_NAME = "omicverse_fast_rank_preallocated_v1"
EXPECTED_OMICVERSE_VERSION = "1.7.5"
EXPECTED_SOURCE_HASHES = {
    "geneset_aucell": "fc430a4d6fd4e22686ab844919b9129eda71cb0835048d22b2a965932a9f6ea1",
    "fast_rank": "f2760bd6586cec0a4d5aac8ae39be7313112158ca80414b21b8ce30c9929f8e9",
    "_rank_sparse_row": "e270e03460cb4f8b66fec1d1d21e74eeb370e97701812f53a3cae3f0ecec9c08",
}
_PATCH_LOCK = threading.Lock()


def function_source_sha256(function: Callable[..., Any]) -> str:
    """Hash the exact UTF-8 function text returned by ``inspect.getsource``."""
    return hashlib.sha256(inspect.getsource(function).encode("utf-8")).hexdigest()


def package_version(name: str) -> str:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return "unknown"


def preallocated_fast_rank(
    matrix,
    seed: int = 42,
    *,
    rank_sparse_row: Callable[[Any, int], np.ndarray],
    audit: dict[str, Any],
    progress: Callable[[Any], Any],
) -> np.ndarray:
    """Return the exact OmicVerse ranks while retaining one rank matrix."""
    rng = np.random.default_rng(seed)
    n_rows, n_cols = matrix.shape
    shuffle_order = rng.permutation(n_cols)
    shuffled = matrix[:, shuffle_order]
    output = np.empty((n_rows, n_cols), dtype=np.int32)
    for row in progress(range(n_rows)):
        row_rank = rank_sparse_row(shuffled.getrow(row), n_cols)
        output[row, shuffle_order] = row_rank
    permutation = shuffle_order.astype("<i8", copy=False).tobytes(order="C")
    audit.update({
        "call_count": int(audit.get("call_count", 0)) + 1,
        "seed": int(seed),
        "n_genes": int(n_cols),
        "rank_dtype": str(output.dtype),
        "permutation_sha256": hashlib.sha256(permutation).hexdigest(),
        "permutation_hash_contract": PERMUTATION_HASH_CONTRACT,
    })
    return output


def _provider_components(omicverse):
    scgsea = importlib.import_module("omicverse.single._scgsea")
    aucell = importlib.import_module("omicverse.single._aucell")
    provider = getattr(getattr(omicverse, "single", None), "geneset_aucell", None)
    original_fast_rank = getattr(scgsea, "fast_rank", None)
    canonical_fast_rank = getattr(aucell, "fast_rank", None)
    rank_sparse_row = getattr(aucell, "_rank_sparse_row", None)
    if str(getattr(omicverse, "__version__", "unknown")) != EXPECTED_OMICVERSE_VERSION:
        raise RuntimeError("OmicVerse 1.7.5 AUCell adapter rejected provider version drift")
    if provider is not getattr(scgsea, "geneset_aucell", None):
        raise RuntimeError("OmicVerse 1.7.5 AUCell adapter rejected public scorer object drift")
    if original_fast_rank is not canonical_fast_rank:
        raise RuntimeError("OmicVerse 1.7.5 AUCell adapter rejected fast_rank object drift")
    functions = {
        "geneset_aucell": provider,
        "fast_rank": original_fast_rank,
        "_rank_sparse_row": rank_sparse_row,
    }
    try:
        observed_hashes = {name: function_source_sha256(function)
                           for name, function in functions.items()}
    except (OSError, TypeError) as exc:
        raise RuntimeError(
            "OmicVerse 1.7.5 AUCell adapter could not audit provider source"
        ) from exc
    if observed_hashes != EXPECTED_SOURCE_HASHES:
        changed = sorted(name for name in EXPECTED_SOURCE_HASHES
                         if observed_hashes.get(name) != EXPECTED_SOURCE_HASHES[name])
        raise RuntimeError(
            f"OmicVerse 1.7.5 AUCell adapter rejected source drift: {changed}"
        )
    progress = original_fast_rank.__globals__.get("tqdm", lambda values: values)
    return scgsea, provider, original_fast_rank, rank_sparse_row, progress, observed_hashes


@contextmanager
def audited_fast_rank(omicverse):
    """Temporarily replace exactly one audited OmicVerse module global."""
    if not _PATCH_LOCK.acquire(blocking=False):
        raise RuntimeError("OmicVerse 1.7.5 AUCell fast_rank adapter is already active")
    scgsea = None
    original = None
    replacement = None
    audit: dict[str, Any] = {}
    body_failed = False
    try:
        (scgsea, provider, original, rank_sparse_row, progress,
         observed_hashes) = _provider_components(omicverse)
        audit.update({
            "name": ADAPTER_NAME,
            "active": True,
            "omicverse_version": EXPECTED_OMICVERSE_VERSION,
            "ctxcore_version": package_version("ctxcore"),
            "public_scorer": "omicverse.single.geneset_aucell",
            "source_hash_contract": SOURCE_HASH_CONTRACT,
            "original_source_hashes": observed_hashes,
            "adapter_source_hash": function_source_sha256(preallocated_fast_rank),
            "full_cohort_single_provider_call": True,
            "call_count": 0,
        })
        replacement = partial(
            preallocated_fast_rank,
            rank_sparse_row=rank_sparse_row,
            audit=audit,
            progress=progress,
        )
        scgsea.fast_rank = replacement
        try:
            yield audit
        except BaseException:
            body_failed = True
            raise
    finally:
        restore_error = None
        if scgsea is not None and original is not None:
            if getattr(scgsea, "fast_rank", None) is not replacement:
                restore_error = RuntimeError(
                    "OmicVerse 1.7.5 AUCell fast_rank changed while adapter was active"
                )
            scgsea.fast_rank = original
            audit["restored_original"] = True
        _PATCH_LOCK.release()
        if restore_error is not None and not body_failed:
            raise restore_error
