import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise_analysis.methods.partition import compare_membership, compare_partitions
from revise_analysis.methods.partition import compute_partition
from revise_analysis.methods.regions import compute_window_diversity
from revise_analysis.methods.spatial import compute_moran, native_knn_graph
from revise_analysis.methods import programs


def test_compare_membership_uses_only_intersecting_identifiers_without_global_pairing():
    raw = pd.Series(["a", "a", "b"], index=["raw-only", "shared-2", "shared-1"])
    svc = pd.Series(["y", "x", "y"], index=["svc-only", "shared-1", "shared-2"])
    result = compare_membership(raw, svc)
    assert result.assignments.index.tolist() == ["shared-2", "shared-1"]
    assert result.summary.loc[0, "n_units"] == 2
    assert result.assignments.loc["shared-2", "raw_cluster"] == "a"
    assert result.assignments.loc["shared-1", "recon_cluster"] == "x"


def test_partition_comparison_retains_strict_paired_contract_for_regressions():
    with pytest.raises(ValueError, match="same unit IDs"):
        compare_partitions(pd.Series(["a"], index=["raw"]), pd.Series(["a"], index=["svc"]))


def test_membership_without_shared_identifiers_is_an_explicit_empty_comparison():
    result = compare_membership(pd.Series(["a"], index=["raw"]), pd.Series(["b"], index=["svc"]))
    assert result.assignments.empty
    assert result.summary.loc[0, "status"] == "no_common_units"


def test_native_knn_is_symmetric_binary_and_excludes_self_edges():
    graph = native_knn_graph(np.array([[0., 0.], [1., 0.], [2., 0.], [3., 0.]]), n_neighbors=1)
    assert (graph != graph.T).nnz == 0
    assert graph.diagonal().sum() == 0
    assert set(graph.data) == {1.0}


def test_moran_keeps_constant_gene_explicit_instead_of_zero(monkeypatch):
    def normalize_total(adata, **_kwargs):
        totals = np.asarray(adata.X.sum(axis=1)).ravel()
        totals[totals == 0] = 1
        adata.X = np.asarray(adata.X, dtype=float) / totals[:, None] * 1e4
    def log1p(adata):
        adata.X = np.log1p(adata.X)
    monkeypatch.setitem(sys.modules, "scanpy", SimpleNamespace(pp=SimpleNamespace(normalize_total=normalize_total, log1p=log1p)))
    adata = AnnData(np.array([[1., 5., 4.], [2., 5., 3.], [3., 5., 2.], [4., 5., 1.]]), obs=pd.DataFrame(index=list("abcd")), var=pd.DataFrame(index=["varying", "constant", "balancer"]))
    adata.obsm["spatial"] = np.c_[np.arange(4), np.zeros(4)]
    result = compute_moran(adata, n_neighbors=1, min_units=3)
    assert result.loc["constant", "status"] == "not_computable"
    assert result.loc["constant", "reason"] == "constant_expression"
    assert pd.isna(result.loc["constant", "moran"])
    assert result.loc["varying", "status"] == "computed"


def test_formal_methods_reject_negative_expression_before_transformation(monkeypatch):
    adata = AnnData(np.array([[-1., 1.], [1., 1.], [1., 1.]]), obs=pd.DataFrame(index=list("abc")), var=pd.DataFrame(index=["g1", "g2"]))
    adata.obsm["spatial"] = np.c_[np.arange(3), np.zeros(3)]
    with pytest.raises(ValueError, match="nonnegative"):
        compute_moran(adata)
    with pytest.raises(ValueError, match="nonnegative"):
        compute_partition(adata)
    with pytest.raises(ValueError, match="nonnegative"):
        programs.compute_pathway_scores(adata, ["g1"])


def test_partition_with_two_genes_exposes_one_component_embedding(monkeypatch):
    def normalize_total(*_args, **_kwargs):
        return None
    def log1p(*_args, **_kwargs):
        return None
    def pca(adata, *, n_comps, **_kwargs):
        assert n_comps == 1
        adata.obsm["X_pca"] = np.ones((adata.n_obs, 1))
    def neighbors(adata, **_kwargs):
        adata.obsp["connectivities"] = np.eye(adata.n_obs)
    def leiden(adata, *, key_added, **_kwargs):
        adata.obs[key_added] = ["0", "1", "1"]
    fake_scanpy = SimpleNamespace(pp=SimpleNamespace(normalize_total=normalize_total, log1p=log1p, neighbors=neighbors), tl=SimpleNamespace(pca=pca, leiden=leiden))
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    adata = AnnData(np.ones((3, 2)), obs=pd.DataFrame(index=list("abc")), var=pd.DataFrame(index=["g1", "g2"]))
    result = compute_partition(adata)
    assert result["labels"].index.tolist() == ["a", "b", "c"]
    assert result["embedding"].columns.tolist() == ["PC1"]


def test_single_side_window_diversity_is_deterministic_and_does_not_need_other_side_ids():
    coordinates = pd.DataFrame({"x": [0., 1., 2., 40.], "y": [0., 0., 0., 0.]}, index=["svc-a", "svc-b", "svc-c", "svc-d"])
    labels = pd.Series(["x", "x", "y", "z"], index=coordinates.index)
    first = compute_window_diversity(coordinates, labels, window_side_length=10., min_units=3, n_draws=30, random_state=42)
    second = compute_window_diversity(coordinates, labels, window_side_length=10., min_units=3, n_draws=30, random_state=42)
    pd.testing.assert_frame_equal(first, second)
    valid = first.set_index("window_id").loc["0_0"]
    invalid = first.set_index("window_id").loc["4_0"]
    assert valid["neff"] > 1
    assert bool(valid["valid_window"])
    assert not bool(invalid["valid_window"])
    assert pd.isna(invalid["neff"])


def test_pathway_scores_do_not_substitute_a_fallback_when_aucell_provider_is_missing(monkeypatch):
    adata = AnnData(np.ones((2, 1)), obs=pd.DataFrame(index=["a", "b"]), var=pd.DataFrame(index=["G1"]))
    def absent(_name):
        raise ModuleNotFoundError("No module named 'omicverse'", name="omicverse")
    monkeypatch.setattr(programs.importlib, "import_module", absent)
    with pytest.raises(ImportError, match="OmicVerse is required"):
        programs.compute_pathway_scores(adata, ["G1"])
