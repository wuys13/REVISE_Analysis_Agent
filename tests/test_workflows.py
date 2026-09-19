from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from revise_analysis.io import Sample
from revise_analysis.analyses import reconstruction_impact as impact
from revise_analysis.analyses import spatial_autocorrelation as moran
from revise_analysis.analyses import pathway_activity as pathways


def _sample() -> Sample:
    raw = ad.AnnData(np.arange(24, dtype=float).reshape(6, 4))
    raw.obs_names = [f"u{i}" for i in range(6)]
    raw.var_names = [f"g{i}" for i in range(4)]
    raw.obs["Level1"] = ["Fibroblast", "T", "Mono_Macro", "T", "Fibroblast", "Other"]
    raw.obs["Level2"] = ["A", "B", "C", "B", "A", "D"]
    raw.obsm["spatial"] = np.c_[np.arange(6), np.arange(6)]
    svc = raw[:5].copy()
    svc.X = svc.X + 1
    svc.obs["SVC_cluster"] = [str(i % 2) for i in range(svc.n_obs)]
    return Sample("unit", raw, svc, {"expression": {side: {"identity": "test_expression", "scale": "untransformed_nonnegative"} for side in ("raw", "svc")}, "label_aliases": {"Mono_Macro": "Mono/Macro"}, "spatial": {"microns_per_coordinate": 1}}, Path("sample.yaml"))


def test_impact_keeps_sides_independent_and_membership_is_explicit(tmp_path, monkeypatch):
    sample = _sample()

    def fake_partition(adata, *, resolution, n_top_genes, random_state, **kwargs):
        labels = pd.Series([str(i % 2) for i in range(adata.n_obs)], index=adata.obs_names, name="partition")
        return {"labels": labels, "summary": {"n_units": adata.n_obs, "n_clusters": 2}}

    monkeypatch.setattr(impact, "compute_partitions", fake_partition)
    from revise_analysis.analyses import spatial_autocorrelation as spatial_module
    monkeypatch.setattr(spatial_module, "compute_moran", lambda adata, **_: pd.DataFrame({"gene_id": adata.var_names, "moran": [0.1] * adata.n_vars, "status": "ok", "reason": ""}))
    result = impact.run(sample, tmp_path, {"membership_same_units": True, "window_side_microns": 100})

    assert result["status"] == "partial"  # missing scopes are reported, not silently erased
    assert (tmp_path / "tables/raw_anatomy_context.csv").exists()
    assert (tmp_path / "tables/membership_All.csv").exists()
    assert (tmp_path / "tables/spatial_changed_map_All.csv").exists()
    assert (tmp_path / "tables/gain_All_common_valid_windows.csv").exists()
    assert (tmp_path / "tables/window_diversity_raw_level2_baseline_All.csv").exists()
    assert (tmp_path / "tables/raw_level2_baseline_vs_svc_All_common_valid_windows.csv").exists()
    assert (tmp_path / "figures/raw_anatomy_context.png").exists()
    assert (tmp_path / "figures/spatial_changed_map_All.png").exists()


def test_membership_scope_is_defined_by_raw_not_svc_broad_annotation(tmp_path, monkeypatch):
    sample = _sample()
    sample.raw.obs["Level1"] = "Fibroblast"
    sample.svc.obs["Level1"] = "Other"  # a changed broad identity must remain eligible for Raw-defined membership

    def fake_partition(adata, **_):
        return {"labels": pd.Series(["0"] * adata.n_obs, index=adata.obs_names), "summary": {"n_units": adata.n_obs, "n_clusters": 1}}

    monkeypatch.setattr(impact, "compute_partitions", fake_partition)
    from revise_analysis.analyses import spatial_autocorrelation as spatial_module
    monkeypatch.setattr(spatial_module, "compute_moran", lambda adata, **_: pd.DataFrame({"gene_id": adata.var_names, "moran": [0.1] * adata.n_vars, "status": "ok", "reason": ""}))
    impact.run(sample, tmp_path, {"membership_same_units": True})
    assert (tmp_path / "tables/membership_Fibroblast.csv").exists()


def test_impact_anatomy_uses_full_raw_level1_when_level2_is_absent(tmp_path, monkeypatch):
    sample = _sample()
    del sample.raw.obs["Level2"]
    del sample.svc.obs["Level2"]

    monkeypatch.setattr(impact, "compute_partitions", lambda adata, **_: {
        "labels": pd.Series(["0"] * adata.n_obs, index=adata.obs_names),
        "summary": {"n_units": adata.n_obs, "n_clusters": 1},
    })
    from revise_analysis.analyses import spatial_autocorrelation as spatial_module
    monkeypatch.setattr(spatial_module, "compute_moran", lambda adata, **_: pd.DataFrame({
        "gene_id": adata.var_names, "moran": [0.1] * adata.n_vars, "status": "ok", "reason": "",
    }))

    result = impact.run(sample, tmp_path, {"window_side_microns": 100})
    anatomy = pd.read_csv(tmp_path / "tables/raw_anatomy_context.csv")
    assert set(anatomy["broad_label"]) == set(sample.labels("raw"))
    assert not (tmp_path / "tables/window_diversity_raw_level2_baseline_All.csv").exists()
    assert any(item["component"] == "raw_level2_baseline" for item in result["unavailable"])


def test_impact_writes_state_and_positive_gain_masks_when_threshold_is_accepted(tmp_path, monkeypatch):
    sample = _sample()
    sample.raw.obs["Level1"] = "Parent"
    sample.svc.obs["Level1"] = "Parent"
    def fake_partition(adata, **_):
        labels = ["0"] * adata.n_obs if adata.n_obs == sample.raw.n_obs else [str(index % 2) for index in range(adata.n_obs)]
        return {"labels": pd.Series(labels, index=adata.obs_names),
                "summary": {"n_units": adata.n_obs, "n_clusters": len(set(labels))}}
    monkeypatch.setattr(impact, "compute_partitions", fake_partition)
    from revise_analysis.analyses import spatial_autocorrelation as spatial_module
    from revise_analysis.methods import regions
    monkeypatch.setattr(spatial_module, "compute_moran", lambda adata, **_: pd.DataFrame({
        "gene_id": adata.var_names, "moran": [0.1] * adata.n_vars, "status": "ok", "reason": "",
    }))
    monkeypatch.setattr(regions, "select_region_threshold", lambda values, **_: (
        {"status": "ok", "threshold": 0.01, "n_windows": len(values), "n_valid_bootstrap": 1},
        pd.DataFrame({"bootstrap_threshold": [0.01]}),
    ))

    impact.run(sample, tmp_path, {
        "scopes": ["Parent"], "window_side_microns": 100, "min_window_units": 2,
        "n_window_draws": 5,
    })
    state = pd.read_csv(tmp_path / "tables/state_Parent_region_windows.csv")
    gain = pd.read_csv(tmp_path / "tables/gain_Parent_region_windows.csv")
    extent = pd.read_csv(tmp_path / "tables/gain_Parent_region_extent.csv")
    assert state["in_region"].any()
    assert gain["in_region"].any()
    assert extent.loc[0, "window_side_microns"] == 100
    assert (tmp_path / "figures/gain_Parent_region_windows.png").exists()


def test_moran_and_pathway_are_independent_component_workflows(tmp_path, monkeypatch):
    sample = _sample()
    seen_neighbors = []
    def fake_moran(adata, **kwargs):
        seen_neighbors.append(kwargs["n_neighbors"])
        return pd.DataFrame({"gene_id": adata.var_names, "moran": [0.1] * adata.n_vars, "status": "ok", "reason": ""})
    monkeypatch.setattr(moran, "compute_moran", fake_moran)
    moran_result = moran.run(sample, tmp_path / "moran", {"moran_n_neighbors": 2})
    assert moran_result["status"] == "succeeded"
    assert (tmp_path / "moran/tables/moran_raw.csv").exists()
    assert seen_neighbors == [2, 2]

    def fake_scores(adata, *, genes, score_name, auc_threshold, seed):
        return {"scores": pd.Series(np.arange(adata.n_obs), index=adata.obs_names, name=score_name), "coverage": {"n_present": 1}, "metadata": {"provider": "test"}}

    monkeypatch.setattr(pathways, "compute_pathway_scores", fake_scores)
    pathway_result = pathways.run(sample, tmp_path / "programs", {"gene_sets": {"Test": ["g1"]}})
    assert pathway_result["status"] == "succeeded"
    assert (tmp_path / "programs/tables/pathway_svc_Test.csv").exists()


def test_uncomputable_moran_genes_are_not_a_complete_success(tmp_path, monkeypatch):
    monkeypatch.setattr(moran, "compute_moran", lambda adata, **_: pd.DataFrame({
        "gene_id": adata.var_names, "moran": float("nan"),
        "status": "not_computable", "reason": "constant_expression",
    }))
    result = moran.run(_sample(), tmp_path)
    assert result["status"] == "partial"
    assert all("constant_expression" in item["reason"] for item in result["unavailable"])


def test_pathway_no_overlap_is_explicit_without_provider(tmp_path):
    result = pathways.run(_sample(), tmp_path, {"gene_sets": {"Absent": ["not-a-gene"]}})
    assert result["status"] == "skipped"
    assert len(result["unavailable"]) == 2
    assert all("no overlap" in item["reason"] for item in result["unavailable"])


def test_native_sampling_uses_distinct_side_random_streams(tmp_path, monkeypatch):
    sample = _sample()
    sample.svc = sample.raw.copy()
    seen = []
    def record(adata, **_):
        seen.append(adata.obs_names.tolist())
        return pd.DataFrame({"gene_id": adata.var_names, "moran": 0.1, "status": "computed"})
    monkeypatch.setattr(moran, "compute_moran", record)
    moran.run(sample, tmp_path, {"sample_n_units": 3})
    assert seen[0] != seen[1]


def test_artifact_name_collisions_are_explicit(tmp_path):
    import pytest
    with pytest.raises(ValueError, match="collide"):
        impact.run(_sample(), tmp_path, {"scopes": ["A/B", "A_B"]})
    with pytest.raises(ValueError, match="collide"):
        pathways.run(_sample(), tmp_path, {"gene_sets": {"A/B": ["g1"], "A_B": ["g2"]}})


def test_unexpected_pathway_provider_valueerror_is_not_unavailable(tmp_path, monkeypatch):
    import pytest
    def broken_provider(*args, **kwargs):
        raise ValueError('provider execution bug')
    monkeypatch.setattr(pathways, 'compute_pathway_scores', broken_provider)
    with pytest.raises(ValueError, match='provider execution bug'):
        pathways.run(_sample(), tmp_path, {'gene_sets': {'Test': ['g1']}})
