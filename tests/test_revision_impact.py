from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from revise_analysis.analyses import reconstruction_impact as impact
from revise_analysis.analyses.impact_tables import changed_unit_summary, common_valid_delta, label_composition
from revise_analysis.io import Sample


def _sample(*, reconstruction=True, unknown_expression=False):
    n = 12
    raw = ad.AnnData(np.arange(n * 5, dtype=float).reshape(n, 5) + 1)
    raw.obs_names = [f"u{i}" for i in range(n)]
    raw.var_names = [f"g{i}" for i in range(5)]
    raw.obs["Level1"] = ["Tumor", "Intestinal Epithelial"] * 6
    raw.obs["Level2"] = ["L0", "L1", "L2"] * 4
    raw.obsm["spatial"] = np.c_[np.arange(n) * 10, np.zeros(n)]
    svc = raw.copy()
    svc.obs["Level1"] = "Fibroblast"
    if reconstruction:
        svc.obs["SVC_cluster"] = ["A", "B", "C"] * 4
    expression = ({side: {"identity": "unknown", "scale": "unknown"} for side in ("raw", "svc")}
                  if unknown_expression else
                  {side: {"identity": "test", "scale": "untransformed_nonnegative"} for side in ("raw", "svc")})
    config = {
        "expression": expression,
        "label_aliases": {"A": "aliased-broad-only"},
        "spatial": {"microns_per_coordinate": 1},
    }
    return Sample("unit", raw, svc, config, Path("sample.yaml"))


def test_reconstruction_labels_are_not_rewritten_by_broad_aliases(tmp_path):
    workflow = impact.ImpactWorkflow(_sample(), tmp_path)
    workflow.run_stage("input")
    assert set(workflow.state["svc_labels"]) == {"A", "B", "C"}
    saved = pd.read_csv(tmp_path / "tables/reconstruction_label_summary.csv")
    assert set(saved.reconstruction_label) == {"A", "B", "C"}
    assert set(saved.label_role) == {"reconstructed_state_label"}


def test_input_manifest_preserves_carrier_provenance(tmp_path):
    sample = _sample()
    sample.config["files"] = {"raw": "raw.h5ad", "svc": "SVC.h5ad"}
    sample.config["metadata"] = {"svc_carrier": "temporary spatial-label reconstruction carrier"}
    sample.config["provenance"] = {"svc": {"source_path": "/source/spatial.h5ad", "sha256": "abc"}}
    workflow = impact.ImpactWorkflow(sample, tmp_path)
    workflow.run_stage("input")
    manifest = json.loads((tmp_path / "tables/input_manifest.json").read_text())
    assert manifest["metadata"]["svc_carrier"].startswith("temporary")
    assert manifest["provenance"]["svc"]["sha256"] == "abc"
    assert manifest["resolved_files"]["svc"].endswith("SVC.h5ad")


@pytest.mark.parametrize("parameters", [
    {"scopes": "All"}, {"membership_same_units": "false"},
    {"n_window_draws": 2.5}, {"sample_n_units": 0},
    {"anatomy_window_side_microns": float("nan")},
    {"parent_window_side_microns": {"All": "40"}},
])
def test_effective_parameters_reject_ambiguous_values(parameters):
    with pytest.raises(ValueError):
        impact.effective_parameters(_sample(), parameters)


def test_anatomy_scale_does_not_change_parent_diversity(tmp_path):
    outputs = []
    for anatomy in (20, 80):
        workflow = impact.ImpactWorkflow(
            _sample(unknown_expression=True), tmp_path / str(anatomy),
            {"scopes": ["All"], "anatomy_window_side_microns": anatomy,
             "parent_window_side_microns": {"All": 40}, "min_window_units": 2,
             "n_window_draws": 10},
        )
        for stage in ("input", "baseline", "support", "diversity"):
            workflow.run_stage(stage)
        outputs.append(workflow.state["windows"]["svc", "All"])
    pd.testing.assert_frame_equal(outputs[0], outputs[1])


def test_threshold_failure_keeps_continuous_state_table(tmp_path):
    workflow = impact.ImpactWorkflow(
        _sample(unknown_expression=True), tmp_path,
        {"scopes": ["All"], "parent_window_side_microns": {"All": 40},
         "min_window_units": 2, "n_window_draws": 5},
    )
    for stage in ("input", "baseline", "support", "diversity", "regions"):
        workflow.run_stage(stage)
    table = pd.read_csv(tmp_path / "tables/state_All_region_windows.csv")
    assert "neff" in table
    assert table.neff.notna().any()
    assert table.in_region.isna().all()
    assert (tmp_path / "tables/state_All_threshold.json").exists()


def test_supported_outside_and_unsupported_are_distinct(tmp_path, monkeypatch):
    from revise_analysis.methods import regions
    monkeypatch.setattr(regions, "select_region_threshold", lambda values, **kwargs: (
        {"status": "ok", "threshold": 2.0, "n_windows": 2, "n_valid_bootstrap": 1},
        pd.DataFrame(),
    ))
    workflow = impact.ImpactWorkflow(_sample(), tmp_path, {"scopes": ["All"]})
    table = pd.DataFrame({
        "window_id": ["a", "b", "c"], "valid_window": [True, True, False],
        "n_units": [4, 4, 1], "neff": [3.0, 1.0, np.nan],
        "window_x": [0.0, 1.0, 2.0], "window_y": [0.0, 0.0, 0.0],
    })
    flagged = workflow._write_region(table, name="state_All", value_column="neff", side_microns=40)
    assert flagged.loc[0, "in_region"] == True  # noqa: E712
    assert flagged.loc[1, "in_region"] == False  # noqa: E712
    assert pd.isna(flagged.loc[2, "in_region"])
    assert flagged.loc[2, "region_status"] == "insufficient_support"


def test_unexpected_stage_error_retains_prior_artifacts_and_is_failed(tmp_path, monkeypatch):
    workflow = impact.ImpactWorkflow(_sample(), tmp_path, continue_on_error=True)
    monkeypatch.setattr(workflow, "stage_baseline", lambda: (_ for _ in ()).throw(RuntimeError("boom")))
    workflow.run_stage("input")
    workflow.run_stage("baseline")
    result = workflow.result()
    assert result["status"] == "failed"
    assert (tmp_path / "tables/input_overview.csv").exists()
    assert result["stage_errors"][0]["stage"] == "baseline"
    assert (tmp_path / "errors/baseline.txt").exists()


def test_figures_do_not_adopt_stale_tables(tmp_path):
    stale = tmp_path / "tables/gain_stale_common_valid_windows.csv"
    stale.parent.mkdir(parents=True)
    pd.DataFrame({"delta_neff_vs_raw_leiden": [1.0]}).to_csv(stale, index=False)
    workflow = impact.ImpactWorkflow(_sample(), tmp_path)
    workflow.run_stage("figures")
    assert "figures/gain_stale_common_valid_windows.png" not in workflow.outputs


def test_missing_reconstruction_label_is_expected_unavailable(tmp_path):
    workflow = impact.ImpactWorkflow(_sample(reconstruction=False), tmp_path, continue_on_error=True)
    workflow.run_stage("input")
    result = workflow.result()
    assert result["status"] == "partial"
    assert not result["stage_errors"]
    assert result["unavailable"][0]["component"] == "svc_reconstruction_labels"


def test_derived_tables_keep_denominators_and_signed_delta():
    raw = pd.DataFrame({"window_id": ["0_0", "1_0"], "valid_window": [True, True],
                        "neff": [3.0, 2.0], "n_units": [4, 4], "window_x": [1, 2], "window_y": [1, 1]})
    svc = raw.copy()
    svc["neff"] = [2.0, 4.0]
    delta = common_valid_delta(raw, svc, delta_column="delta")
    assert delta.delta.tolist() == [-1.0, 2.0]

    units = pd.DataFrame({"state_region": [True, True, False], "anatomy_region": ["Tumor"] * 3},
                         index=["a", "b", "c"])
    composition = label_composition(units, pd.Series(["x", "y", "x"], index=units.index),
                                    group_columns=["state_region", "anatomy_region"])
    assert set(composition.denominator_units) == {1, 2}

    changed = units.assign(unit_changed=[True, False, True])
    summary = changed_unit_summary(changed, group_columns=["state_region"])
    assert summary.n_compared.sum() == 3


def test_partition_execution_valueerror_is_not_scientific_unavailability(tmp_path, monkeypatch):
    sample = _sample()
    sample.config['expression'] = {'raw': {'identity': 'synthetic'}}
    workflow = impact.ImpactWorkflow(sample, tmp_path, {'scopes': ['All']}, continue_on_error=True)
    monkeypatch.setattr(impact, 'compute_partitions', lambda *a, **k: (_ for _ in ()).throw(ValueError('internal partition bug')))
    workflow.run_stage('baseline')
    result = workflow.result()
    assert result['status'] == 'failed'
    assert result['stage_errors'][0]['error_type'] == 'ValueError'
    assert not any('internal partition bug' in item['reason'] for item in result['unavailable'])


def test_raw_baseline_anatomy_uses_actual_partition_cohort(tmp_path):
    sample = _sample(unknown_expression=True)
    workflow = impact.ImpactWorkflow(sample, tmp_path, {'scopes': ['All']})
    workflow.run_stage('input')
    workflow.run_stage('baseline')
    workflow.state['partitions'] = {('raw', 'All'): pd.Series('0', index=sample.raw.obs_names[:3])}
    workflow.run_stage('support')
    workflow.run_stage('anatomy')
    native = workflow.state['parent_anatomy']['raw', 'All']
    baseline = workflow.state['parent_anatomy']['raw_baseline', 'All']
    assert native.n_units.sum() == sample.raw.n_obs
    assert baseline.n_units.sum() == 3
    assert set(baseline.population) == {'raw_leiden_cohort'}
