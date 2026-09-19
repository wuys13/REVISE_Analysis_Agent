import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise_analysis.methods import partition
from revise_analysis.methods.partition import PreparedPartition, prepare_partition_graph, select_matched_k_partition
from revise_analysis.methods.regions import (
    aggregate_parent_anatomy,
    assign_anatomy_candidates,
    assign_points_to_anatomy,
    assign_square_windows,
    compute_window_diversity,
    select_window_scale,
)
from revise_analysis.methods.spatial import compute_moran


def test_flat_support_curve_has_no_fake_recommendation_but_explicit_scale_still_computes():
    coordinates = pd.DataFrame(
        {"x": [0.0, 1.0, 2.0, 3.0], "y": [0.0, 0.0, 0.0, 0.0]},
        index=list("abcd"),
    )
    recommendation, sensitivity = select_window_scale(
        coordinates, candidate_window_sides=[10, 20, 40], min_parent_units=1, origin=(0, 0)
    )
    assert recommendation == {
        "status": "no_clear_recommendation",
        "window_side_length": None,
        "min_parent_units": 1,
        "reason": "flat_support_curve",
    }
    assert sensitivity.n_valid_windows.tolist() == [1, 1, 1]
    diversity = compute_window_diversity(
        coordinates,
        pd.Series(["x", "x", "y", "y"], index=coordinates.index),
        window_side_length=20,
        origin=(0, 0),
        min_units=4,
        n_draws=5,
    )
    assert diversity.valid_window.tolist() == [True]
    assert diversity.neff.iloc[0] == pytest.approx(2.0)


def test_cross_scale_anatomy_mapping_uses_points_and_retains_unknown_and_ties():
    raw_coordinates = pd.DataFrame(
        {"x": [0.0, 1.0, 20.0], "y": [0.0, 0.0, 0.0]}, index=["r1", "r2", "r3"]
    )
    anatomy_assignments = assign_square_windows(
        raw_coordinates, window_side_length=10, origin=(0, 0)
    )
    anatomy = assign_anatomy_candidates(
        anatomy_assignments,
        pd.Series(["Tumor", "Intestinal Epithelial", "Interface"], index=raw_coordinates.index),
    )
    # Literal upstream Interface is source information, not a window-level
    # Interface candidate; that third window is therefore Other.
    assert anatomy.set_index("window_id").level1_region.to_dict() == {
        "0_0": "Interface",
        "2_0": "Other",
    }

    svc_coordinates = pd.DataFrame(
        {"x": [0.5, 20.5, 100.0], "y": [0.0, 0.0, 0.0]}, index=["s1", "s2", "s3"]
    )
    mapped = assign_points_to_anatomy(
        svc_coordinates,
        anatomy,
        anatomy_window_side_length=10,
        origin=(0, 0),
    )
    assert mapped.anatomy_region.tolist() == ["Interface", "Other", "Unknown"]
    assert mapped.anatomy_covered.tolist() == [True, True, False]

    parent = assign_square_windows(svc_coordinates, window_side_length=50, origin=(0, 0))
    summary = aggregate_parent_anatomy(parent, mapped).set_index("window_id")
    assert summary.loc["0_0", "anatomy_dominant"] == "Tie"
    assert bool(summary.loc["0_0", "anatomy_tied"])
    assert summary.loc["0_0", "anatomy_tied_categories"] == "Interface|Other"
    assert summary.loc["2_0", "anatomy_unknown_fraction"] == 1.0
    assert summary.loc["2_0", "n_anatomy_covered"] == 0


def _fake_scanpy(call_log):
    def filter_cells(*_args, **_kwargs):
        return None

    def filter_genes(*_args, **_kwargs):
        return None

    def highly_variable_genes(adata, *, n_top_genes, flavor):
        call_log.append(("hvg", flavor, np.asarray(adata.X).copy()))
        adata.var["highly_variable"] = [index < n_top_genes for index in range(adata.n_vars)]

    def normalize_total(adata, **_kwargs):
        call_log.append(("normalize", None, np.asarray(adata.X).copy()))

    def log1p(adata):
        call_log.append(("log1p", None, np.asarray(adata.X).copy()))

    def pca(adata, *, n_comps, **_kwargs):
        adata.obsm["X_pca"] = np.ones((adata.n_obs, n_comps))

    def neighbors(adata, **_kwargs):
        adata.obsp["connectivities"] = np.eye(adata.n_obs)

    return SimpleNamespace(
        pp=SimpleNamespace(
            filter_cells=filter_cells,
            filter_genes=filter_genes,
            highly_variable_genes=highly_variable_genes,
            normalize_total=normalize_total,
            log1p=log1p,
            neighbors=neighbors,
        ),
        tl=SimpleNamespace(pca=pca),
    )


@pytest.mark.parametrize(
    ("state", "expected_calls", "expected_preprocessing"),
    [
        ("untransformed_nonnegative", ["normalize", "log1p", "hvg"], "normalize_total_1e4_log1p"),
        ("log1p_nonnegative", ["hvg"], "declared_log1p_no_transform"),
    ],
)
def test_partition_hvg_matches_declared_scale_without_double_transform(
    monkeypatch, state, expected_calls, expected_preprocessing
):
    calls = []
    monkeypatch.setitem(sys.modules, "scanpy", _fake_scanpy(calls))
    values = np.arange(16, dtype=float).reshape(4, 4) + 1
    adata = AnnData(values, obs=pd.DataFrame(index=list("abcd")), var=pd.DataFrame(index=list("wxyz")))
    prepared = prepare_partition_graph(adata, n_top_genes=2, transformation_state=state)
    assert prepared.summary["hvg_flavor"] == "seurat"
    assert prepared.summary["preprocessing"] == expected_preprocessing
    assert [call[0] for call in calls] == expected_calls
    hvg_call = next(call for call in calls if call[0] == "hvg")
    assert hvg_call[1] == "seurat"
    # HVG sees the full axis; feature subsetting happens only afterwards.
    assert hvg_call[2].shape == values.shape
    np.testing.assert_array_equal(adata.X, values)


def test_matched_k_selection_is_exact_first_and_does_not_mutate_prepared_baseline(monkeypatch):
    work = AnnData(np.ones((4, 2)), obs=pd.DataFrame(index=list("abcd")), var=pd.DataFrame(index=["g1", "g2"]))
    prepared = PreparedPartition(
        work=work,
        embedding=None,
        summary={"status": "prepared"},
        random_state=42,
        transformation_state="log1p_nonnegative",
    )
    baseline_obs = work.obs.copy()

    cluster_counts = {0.2: 2, 0.4: 3, 0.6: 3, 0.8: 4}

    def fake_partition(_prepared, *, resolution):
        k = cluster_counts[resolution]
        labels = pd.Series([str(index % k) for index in range(4)], index=work.obs_names)
        return {"labels": labels, "summary": {"n_clusters": k, "resolution": resolution}, "embedding": None}

    monkeypatch.setattr(partition, "partition_prepared_graph", fake_partition)
    result = select_matched_k_partition(
        prepared,
        candidate_resolutions=[0.8, 0.2, 0.6, 0.4],
        target_k=3,
        main_resolution=0.5,
    )
    # 0.4 and 0.6 are equally close; lower resolution breaks the tie.
    assert result["selection"] == {
        "status": "exact",
        "target_k": 3,
        "actual_k": 3,
        "k_delta": 0,
        "resolution": 0.4,
        "main_resolution": 0.5,
    }
    assert result["candidates"].selected.sum() == 1
    pd.testing.assert_frame_equal(work.obs, baseline_obs)


def test_matched_k_nearest_is_explicit_and_moran_log_input_is_not_retransformed(monkeypatch):
    work = AnnData(np.ones((4, 2)), obs=pd.DataFrame(index=list("abcd")), var=pd.DataFrame(index=["g1", "g2"]))
    prepared = PreparedPartition(work, None, {"status": "prepared"}, 42, "log1p_nonnegative")

    def fake_partition(_prepared, *, resolution):
        k = {0.2: 2, 0.8: 4}[resolution]
        return {"labels": pd.Series(["0"] * 4, index=work.obs_names), "summary": {"n_clusters": k}, "embedding": None}

    monkeypatch.setattr(partition, "partition_prepared_graph", fake_partition)
    result = select_matched_k_partition(
        prepared, candidate_resolutions=[0.2, 0.8], target_k=3, main_resolution=0.5
    )
    assert result["selection"]["status"] == "nearest"
    assert result["selection"]["resolution"] == 0.2
    assert result["selection"]["k_delta"] == -1

    def forbidden(*_args, **_kwargs):
        raise AssertionError("declared log1p values must not be transformed again")

    monkeypatch.setitem(
        sys.modules,
        "scanpy",
        SimpleNamespace(pp=SimpleNamespace(normalize_total=forbidden, log1p=forbidden)),
    )
    adata = AnnData(
        np.log1p(np.array([[1.0, 4.0], [2.0, 3.0], [3.0, 2.0], [4.0, 1.0]])),
        obs=pd.DataFrame(index=list("abcd")),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    adata.obsm["spatial"] = np.c_[np.arange(4), np.zeros(4)]
    moran = compute_moran(adata, n_neighbors=1, transformation_state="log1p_nonnegative")
    assert set(moran.transformation_state) == {"log1p_nonnegative"}
    assert set(moran.status) == {"computed"}


def test_unknown_expression_scale_is_a_method_prerequisite_failure():
    adata = AnnData(
        np.ones((4, 2)), obs=pd.DataFrame(index=list("abcd")), var=pd.DataFrame(index=["g1", "g2"])
    )
    adata.obsm["spatial"] = np.c_[np.arange(4), np.zeros(4)]
    with pytest.raises(ValueError, match="transformation_state"):
        prepare_partition_graph(adata, transformation_state="unknown")
    with pytest.raises(ValueError, match="transformation_state"):
        compute_moran(adata, transformation_state="unknown")
