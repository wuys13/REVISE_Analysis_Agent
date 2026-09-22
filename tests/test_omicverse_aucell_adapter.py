import hashlib
import inspect
import linecache
from types import ModuleType, SimpleNamespace
import threading

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from revise_analysis.methods import omicverse_aucell_adapter as adapter
from revise_analysis.methods import programs


def _rank_sparse_row(row_sparse, n_cols):
    row = row_sparse.toarray().ravel()
    sort_idx = np.lexsort((np.arange(n_cols), -row))
    rank_vals = np.empty(n_cols, dtype=np.int32)
    rank_vals[sort_idx] = np.arange(n_cols)
    return rank_vals


def _original_fast_rank(matrix, seed=42):
    rng = np.random.default_rng(seed)
    n_rows, n_cols = matrix.shape
    shuffle_order = rng.permutation(n_cols)
    shuffled = matrix[:, shuffle_order]
    ranks = np.stack([_rank_sparse_row(shuffled.getrow(i), n_cols)
                      for i in range(n_rows)])
    unshuffled = np.empty_like(ranks)
    unshuffled[:, shuffle_order] = ranks
    return unshuffled


@pytest.mark.parametrize("seed", [42, 117])
def test_preallocated_rank_is_bitwise_equal_and_does_not_mutate_csr(seed):
    # Includes an all-zero row, positive ties and an explicitly stored zero.
    matrix = sparse.csr_matrix((
        np.array([2.0, 2.0, 0.0, 5.0, 1.0, 1.0, 3.0]),
        np.array([0, 2, 4, 1, 0, 3, 4]),
        np.array([0, 0, 3, 4, 7]),
    ), shape=(4, 6))
    before = (matrix.data.copy(), matrix.indices.copy(), matrix.indptr.copy())

    audit = {}
    actual = adapter.preallocated_fast_rank(
        matrix, seed=seed, rank_sparse_row=_rank_sparse_row, audit=audit,
        progress=lambda values: values,
    )

    assert actual.dtype == np.int32
    assert np.array_equal(actual, _original_fast_rank(matrix, seed=seed))
    assert all(np.array_equal(left, right) for left, right in zip(
        before, (matrix.data, matrix.indices, matrix.indptr), strict=True
    ))
    permutation = np.random.default_rng(seed).permutation(matrix.shape[1])
    expected_hash = hashlib.sha256(
        permutation.astype("<i8", copy=False).tobytes(order="C")
    ).hexdigest()
    assert audit == {
        "call_count": 1,
        "seed": seed,
        "n_genes": matrix.shape[1],
        "rank_dtype": "int32",
        "permutation_sha256": expected_hash,
        "permutation_hash_contract": "sha256-int64-little-endian-c-order-v1",
    }


def test_adapter_source_hash_uses_inspect_source_utf8_contract():
    expected = hashlib.sha256(
        inspect.getsource(adapter.preallocated_fast_rank).encode("utf-8")
    ).hexdigest()
    assert adapter.function_source_sha256(adapter.preallocated_fast_rank) == expected


def _fake_provider_modules():
    aucell = ModuleType("omicverse.single._aucell")
    aucell.np = np
    aucell_source = (inspect.getsource(_rank_sparse_row) + "\n"
                     + inspect.getsource(_original_fast_rank)
                     .replace("_original_fast_rank", "fast_rank"))
    _exec_module(aucell, aucell_source)
    scgsea = ModuleType("omicverse.single._scgsea")
    scgsea.fast_rank = aucell.fast_rank
    _exec_module(scgsea,
        "def geneset_aucell(adata, geneset_name, geneset, AUC_threshold=0.01, seed=42):\n"
        "    adata.ranks = fast_rank(adata.X, seed=seed)\n"
        "    if getattr(adata, 'fail', False):\n"
        "        raise RuntimeError('provider failed')\n"
        "    if hasattr(adata, 'obs'):\n"
        "        adata.obs[f'{geneset_name}_aucell'] = adata.ranks[:, 0].astype(float)\n",
    )
    omicverse = SimpleNamespace(
        __version__="1.7.5",
        single=SimpleNamespace(geneset_aucell=scgsea.geneset_aucell),
    )
    return omicverse, scgsea, aucell


def _exec_module(module, source):
    filename = f"<{module.__name__}-test-provider>"
    source = source.rstrip() + "\n"
    linecache.cache[filename] = (len(source), None, source.splitlines(True), filename)
    exec(compile(source, filename, "exec"), module.__dict__)


def _install_fake_provider(monkeypatch):
    omicverse, scgsea, aucell = _fake_provider_modules()
    modules = {
        "omicverse.single._scgsea": scgsea,
        "omicverse.single._aucell": aucell,
    }
    monkeypatch.setattr(adapter.importlib, "import_module", modules.__getitem__)
    monkeypatch.setattr(adapter, "EXPECTED_SOURCE_HASHES", {
        "geneset_aucell": adapter.function_source_sha256(scgsea.geneset_aucell),
        "fast_rank": adapter.function_source_sha256(aucell.fast_rank),
        "_rank_sparse_row": adapter.function_source_sha256(aucell._rank_sparse_row),
    })
    monkeypatch.setattr(adapter, "package_version", lambda _name: "0.2.2")
    return omicverse, scgsea, aucell


def test_audited_context_patches_one_provider_global_and_restores_identity(monkeypatch):
    omicverse, scgsea, aucell = _install_fake_provider(monkeypatch)
    original = scgsea.fast_rank
    data = SimpleNamespace(X=sparse.csr_matrix([[0, 2, 2], [1, 0, 3]]))

    with adapter.audited_fast_rank(omicverse) as audit:
        assert scgsea.fast_rank is not original
        scgsea.geneset_aucell(data, "EMT", ["g1"], seed=117)

    assert scgsea.fast_rank is original
    assert np.array_equal(data.ranks, aucell.fast_rank(data.X, seed=117))
    assert audit["name"] == "omicverse_fast_rank_preallocated_v1"
    assert audit["active"] is True
    assert audit["restored_original"] is True
    assert audit["omicverse_version"] == "1.7.5"
    assert audit["ctxcore_version"] == "0.2.2"
    assert audit["public_scorer"] == "omicverse.single.geneset_aucell"
    assert audit["source_hash_contract"] == "sha256-inspect-getsource-utf8-v1"
    assert audit["original_source_hashes"] == adapter.EXPECTED_SOURCE_HASHES
    assert audit["adapter_source_hash"] == adapter.function_source_sha256(
        adapter.preallocated_fast_rank
    )
    assert audit["call_count"] == 1
    assert audit["full_cohort_single_provider_call"] is True


def test_audited_context_restores_identity_after_provider_error(monkeypatch):
    omicverse, scgsea, _ = _install_fake_provider(monkeypatch)
    original = scgsea.fast_rank
    data = SimpleNamespace(X=sparse.csr_matrix([[0, 1], [2, 0]]), fail=True)
    with pytest.raises(RuntimeError, match="provider failed"):
        with adapter.audited_fast_rank(omicverse):
            scgsea.geneset_aucell(data, "EMT", ["g1"])
    assert scgsea.fast_rank is original


def test_pinned_omicverse_function_hashes_are_explicit():
    assert adapter.EXPECTED_SOURCE_HASHES == {
        "geneset_aucell": "fc430a4d6fd4e22686ab844919b9129eda71cb0835048d22b2a965932a9f6ea1",
        "fast_rank": "f2760bd6586cec0a4d5aac8ae39be7313112158ca80414b21b8ce30c9929f8e9",
        "_rank_sparse_row": "e270e03460cb4f8b66fec1d1d21e74eeb370e97701812f53a3cae3f0ecec9c08",
    }


def test_audited_context_rejects_nested_and_concurrent_use(monkeypatch):
    omicverse, scgsea, _ = _install_fake_provider(monkeypatch)
    original = scgsea.fast_rank
    errors = []
    with adapter.audited_fast_rank(omicverse):
        with pytest.raises(RuntimeError, match="already active"):
            with adapter.audited_fast_rank(omicverse):
                pass

        thread = threading.Thread(target=lambda: _capture_context_error(omicverse, errors))
        thread.start()
        thread.join()
    assert len(errors) == 1 and "already active" in str(errors[0])
    assert scgsea.fast_rank is original


def _capture_context_error(omicverse, errors):
    try:
        with adapter.audited_fast_rank(omicverse):
            pass
    except BaseException as exc:
        errors.append(exc)


@pytest.mark.parametrize("drift", ["version", "hash", "provider_object", "rank_object"])
def test_audited_context_fails_closed_on_provider_drift(monkeypatch, drift):
    omicverse, scgsea, _ = _install_fake_provider(monkeypatch)
    if drift == "version":
        omicverse.__version__ = "1.7.6"
    elif drift == "hash":
        monkeypatch.setattr(adapter, "EXPECTED_SOURCE_HASHES", {
            **adapter.EXPECTED_SOURCE_HASHES, "fast_rank": "0" * 64,
        })
    elif drift == "provider_object":
        omicverse.single.geneset_aucell = lambda *_args, **_kwargs: None
    else:
        scgsea.fast_rank = lambda *_args, **_kwargs: None
    expected_after_rejection = scgsea.fast_rank

    with pytest.raises(RuntimeError, match="OmicVerse 1.7.5 AUCell"):
        with adapter.audited_fast_rank(omicverse):
            pass
    assert scgsea.fast_rank is expected_after_rejection


def test_compute_pathway_scores_records_audited_adapter_and_preserves_input(monkeypatch):
    omicverse, _, _ = _install_fake_provider(monkeypatch)
    monkeypatch.setattr(programs, "_require_omicverse", lambda: omicverse)
    matrix = sparse.csr_matrix([
        [2, 1, 0, 0, 0, 0],
        [0, 3, 1, 0, 0, 0],
        [1, 0, 4, 1, 0, 0],
        [2, 2, 2, 0, 0, 0],
    ], dtype=float)
    adata = AnnData(matrix.copy(), obs=pd.DataFrame({"kept": list("abcd")}, index=list("wxyz")),
                    var=pd.DataFrame(index=[f"G{i}" for i in range(6)]))
    before = (adata.X.data.copy(), adata.X.indices.copy(), adata.X.indptr.copy(),
              adata.obs.copy(), adata.obs_names.copy(), adata.var_names.copy())

    result = programs.compute_pathway_scores(
        adata, ["G0", "G2"], score_name="EMT", auc_threshold=0.5, seed=42
    )

    memory = result["metadata"]["memory_adapter"]
    assert memory["name"] == "omicverse_fast_rank_preallocated_v1"
    assert memory["call_count"] == 1
    assert memory["restored_original"] is True
    assert memory["seed"] == 42
    assert memory["n_genes"] == 6
    assert result["scores"].index.tolist() == list("wxyz")
    assert all(np.array_equal(left, right) for left, right in zip(
        before[:3], (adata.X.data, adata.X.indices, adata.X.indptr), strict=True
    ))
    pd.testing.assert_frame_equal(before[3], adata.obs)
    assert before[4].equals(adata.obs_names)
    assert before[5].equals(adata.var_names)
