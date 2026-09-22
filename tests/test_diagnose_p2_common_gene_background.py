import importlib.util
from pathlib import Path

from anndata import AnnData
import numpy as np
import pandas as pd
import pytest

from revise_analysis.analyses._shared import deterministic_subset
from revise_analysis.io import Sample


def diagnostic():
    path = Path(__file__).parents[1] / "scripts/diagnose_p2_common_gene_background.py"
    spec = importlib.util.spec_from_file_location("diagnose_p2_common_gene_background", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_sample():
    ids = [f"u{i:02d}" for i in range(12)]
    labels = ["Fibroblast"] * 4 + ["Mono/Macro"] * 4 + ["T"] * 4
    raw = AnnData(
        np.arange(12 * 3, dtype=float).reshape(12, 3) + 1,
        obs=pd.DataFrame({"Level1": labels}, index=ids),
        var=pd.DataFrame(index=["b", "raw_only", "a"]),
    )
    svc = AnnData(
        np.arange(12 * 3, dtype=float).reshape(12, 3) + 2,
        obs=pd.DataFrame({"Level1": labels}, index=ids),
        var=pd.DataFrame(index=["svc_only", "a", "b"]),
    )
    raw.obsm["spatial"] = np.c_[np.arange(12), np.zeros(12)]
    svc.obsm["spatial"] = raw.obsm["spatial"].copy()
    config = {
        "expression": {
            "raw": {"identity": "measured_expression"},
            "svc": {"identity": "reconstructed_expression"},
        },
        "columns": {"broad": "Level1", "reconstruction": "SVC_cluster"},
        "spatial": {"key": "spatial"},
    }
    return Sample("P2CRC_Xenium", raw, svc, config, Path("sample.yaml"))


def cohort_manifest(module, sample, parameters):
    records = []
    for scope in module.SCOPES:
        for side in module.SIDES:
            eligible = module._scope(sample, side, scope)
            seed = parameters["random_state"] + (side == "svc")
            selected = deterministic_subset(eligible, parameters["sample_n_units"], seed)
            records.append(module.actual_cohort_record(
                side=side, scope=scope, eligible=eligible, selected=selected, seed=seed,
            ))
    pathway_scope = "|".join(module.SCOPES)
    for side in module.SIDES:
        eligible = getattr(sample, side)
        seed = parameters["random_state"] + (side == "svc")
        selected = deterministic_subset(eligible, parameters["sample_n_units"], seed)
        records.append(module.actual_cohort_record(
            side=side, scope=pathway_scope, eligible=eligible, selected=selected, seed=seed,
            component="pathway", method_random_state=parameters["random_state"],
        ))
    return {"hash_contract": module.UNIT_ID_HASH_CONTRACT, "records": records}


def test_reproduces_native_cohorts_and_rejects_ordered_id_drift():
    module = diagnostic()
    sample = make_sample()
    parameters = {"random_state": 42, "sample_n_units": 5, "moran_n_neighbors": 6}
    manifest = cohort_manifest(module, sample, parameters)

    cohorts, records = module.reproduce_cohorts(sample, parameters, manifest)

    assert cohorts["raw", "All"].n_obs == 5
    assert cohorts["svc", "All"].n_obs == 5
    all_records = {row["side"]: row for row in records if row["scope"] == "All"}
    assert all_records["raw"]["selection_random_state"] == 42
    assert all_records["svc"]["selection_random_state"] == 43
    for scope in module.SCOPES[1:]:
        assert cohorts["raw", scope].obs_names.tolist() == cohorts["svc", scope].obs_names.tolist()
    assert len(records) == 8

    manifest["records"][0]["selected_unit_ids_sha256"] = "0" * 64
    with pytest.raises(RuntimeError, match="cohort drift"):
        module.reproduce_cohorts(sample, parameters, manifest)


def test_common_axis_is_imposed_before_existing_method_and_nan_is_preserved():
    module = diagnostic()
    sample = make_sample()
    parameters = {"random_state": 42, "sample_n_units": 5, "moran_n_neighbors": 6}
    manifest = cohort_manifest(module, sample, parameters)
    cohorts, records = module.reproduce_cohorts(sample, parameters, manifest)
    calls = []

    def fake_moran(adata, **kwargs):
        calls.append((adata.var_names.tolist(), adata.obs_names.tolist(), kwargs))
        assert adata.var_names.tolist() == ["a", "b"]
        return pd.DataFrame({
            "gene_id": ["a", "b"],
            "moran": [0.25, np.nan],
            "status": ["computed", "not_computable"],
            "reason": ["", "constant_expression"],
            "n_units": adata.n_obs,
            "n_edges": 2,
            "transformation_state": "untransformed_nonnegative",
        }).set_index("gene_id", drop=False)

    tables = module.compute_common_axis_tables(
        sample=sample, cohorts=cohorts, records=records, parameters=parameters,
        genes=["a", "b"], moran_function=fake_moran,
    )

    assert len(calls) == 8
    assert all(axis == ["a", "b"] for axis, _, _ in calls)
    assert {call[2]["n_neighbors"] for call in calls} == {6}
    assert {call[2]["transformation_state"] for call in calls} == {
        "untransformed_nonnegative"
    }
    assert tables["raw", "All"].loc[1, "reason"] == "constant_expression"
    assert pd.isna(tables["raw", "All"].loc[1, "moran"])

    joint = module.joint_table(tables["raw", "All"], tables["svc", "All"], ["a", "b"], "All")
    assert joint["comparison_available"].tolist() == [True, False]
    assert pd.isna(joint.loc[1, "delta_moran"])
    assert joint.loc[1, "raw_moran"] != 0


def test_common_gene_axis_has_one_explicit_order_and_rejects_count_drift(monkeypatch):
    module = diagnostic()
    sample = make_sample()
    monkeypatch.setattr(module, "EXPECTED_COMMON_GENES", 2)
    assert module.common_gene_axis(sample) == ["a", "b"]
    monkeypatch.setattr(module, "EXPECTED_COMMON_GENES", 3)
    with pytest.raises(RuntimeError, match="common gene axis drift"):
        module.common_gene_axis(sample)


def test_pathway_reproduces_union_cohort_and_common_axis_before_provider():
    module = diagnostic()
    sample = make_sample()
    parameters = {
        "random_state": 42, "sample_n_units": 5, "moran_n_neighbors": 6,
        "pathway_auc_threshold": 0.01,
    }
    manifest = cohort_manifest(module, sample, parameters)
    cohorts, records = module.reproduce_pathway_cohorts(sample, parameters, manifest)
    calls = []

    def fake_pathway(adata, **kwargs):
        calls.append((adata.var_names.tolist(), adata.obs_names.tolist(), kwargs))
        assert adata.var_names.tolist() == ["a", "b"]
        requested = list(dict.fromkeys(kwargs["genes"]))
        return {
            "scores": pd.Series(np.linspace(0.1, 0.5, adata.n_obs), index=adata.obs_names),
            "coverage": {
                "requested_genes": requested, "n_requested": len(requested),
                "present_genes": ["a"], "n_present": 1,
                "n_missing": len(requested) - 1, "fraction_present": 1 / len(requested),
            },
            "metadata": {
                "provider": "test", "rank_cutoff": 2,
                "memory_adapter": {"call_count": 1, "restored_original": True},
            },
        }

    outputs = module.compute_common_axis_pathway(
        sample=sample, cohorts=cohorts, records=records, parameters=parameters,
        genes=["a", "b"], emt_genes=["a", "missing"],
        pathway_function=fake_pathway,
    )

    assert len(calls) == 2
    assert {call[2]["seed"] for call in calls} == {42}
    assert {call[2]["auc_threshold"] for call in calls} == {0.01}
    assert outputs["raw"]["scores"]["unit_id"].tolist() == cohorts["raw"].obs_names.tolist()
    assert outputs["svc"]["scores"]["unit_id"].tolist() == cohorts["svc"].obs_names.tolist()
    assert outputs["raw"]["metadata"]["coverage"]["n_present"] == 1
    assert outputs["raw"]["metadata"]["provider_metadata"]["memory_adapter"] == {
        "call_count": 1, "restored_original": True,
    }
    assert outputs["raw"]["metadata"]["detection_distribution"]["n_units"] == 5


def test_pathway_invalid_cutoff_is_explicit_unavailable_without_zero_scores():
    module = diagnostic()
    sample = make_sample()
    parameters = {
        "random_state": 42, "sample_n_units": 5,
        "pathway_auc_threshold": 0.01,
    }
    manifest = cohort_manifest(module, sample, parameters)
    cohorts, records = module.reproduce_pathway_cohorts(sample, parameters, manifest)

    def invalid_cutoff(*_args, **_kwargs):
        raise module.PrerequisiteUnavailable(
            "AUCell unavailable: detection-derived rank cutoff must be between 2 and n_genes-1"
        )

    outputs = module.compute_common_axis_pathway(
        sample=sample, cohorts=cohorts, records=records, parameters=parameters,
        genes=["a", "b"], emt_genes=["a"], pathway_function=invalid_cutoff,
    )

    for side in module.SIDES:
        assert outputs[side]["status"] == "unavailable"
        assert outputs[side]["scores"] is None
        assert "rank cutoff" in outputs[side]["reason"]
        assert len(outputs[side]["detections"]) == 5
        assert outputs[side]["metadata"]["coverage"]["n_present"] == 1
        assert outputs[side]["metadata"]["detection_derived_cutoff"]["valid"] is False


def test_delta_change_keeps_unavailable_values_as_nan():
    module = diagnostic()
    baseline = pd.DataFrame({
        "gene_id": ["a", "b"],
        "raw_moran": [0.2, np.nan],
        "svc_moran": [0.1, 0.3],
        "delta_moran": [-0.1, np.nan],
        "comparison_available": [True, False],
    })
    common = pd.DataFrame({
        "gene_id": ["a", "b"],
        "raw_moran": [0.4, np.nan],
        "svc_moran": [0.2, np.nan],
        "delta_moran": [-0.2, np.nan],
        "comparison_available": [True, False],
    })

    table = module.delta_change_table(baseline, common, ["a", "b"], "All")

    assert table.loc[0, "delta_moran_axis_change"] == pytest.approx(-0.1)
    assert pd.isna(table.loc[1, "baseline_delta_moran_native_axis"])
    assert pd.isna(table.loc[1, "common_axis_delta_moran"])
    assert pd.isna(table.loc[1, "delta_moran_axis_change"])


def _adapter_metadata(module):
    return {
        "provider": "omicverse",
        "provider_version": "1.7.5",
        "scorer": "single.geneset_aucell",
        "provider_module": "omicverse.single._scgsea",
        "score_name": module.EMT_NAME,
        "score_key": f"{module.EMT_NAME}_aucell",
        "auc_threshold": 0.01,
        "effective_rank_fraction": 0.1,
        "rank_cutoff": 2,
        "seed": 42,
        "memory_adapter": {
            "name": "omicverse_fast_rank_preallocated_v1",
            "active": True,
            "seed": 42,
            "n_genes": 3,
            "rank_dtype": "int32",
            "permutation_sha256": "a" * 64,
            "permutation_hash_contract": "sha256-int64-little-endian-c-order-v1",
            "omicverse_version": "1.7.5",
            "ctxcore_version": "0.2.2",
            "public_scorer": "omicverse.single.geneset_aucell",
            "source_hash_contract": "sha256-inspect-getsource-utf8-v1",
            "original_source_hashes": {"fast_rank": "b" * 64},
            "adapter_source_hash": "c" * 64,
            "call_count": 1,
            "full_cohort_single_provider_call": True,
            "restored_original": True,
        },
    }


def test_baseline_pathway_schema_and_provider_provenance_are_strict(tmp_path):
    module = diagnostic()
    sample = make_sample()
    parameters = {"random_state": 42, "sample_n_units": 5}
    manifest = cohort_manifest(module, sample, parameters)
    _, records = module.reproduce_pathway_cohorts(sample, parameters, manifest)
    record = next(row for row in records if row["side"] == "raw")
    selected = deterministic_subset(sample.raw, 5, 42)
    score_path = tmp_path / "scores.csv"
    pd.DataFrame({
        "unit_id": selected.obs_names,
        "score": np.linspace(0.1, 0.5, 5),
        "side": "raw",
        "program": module.EMT_NAME,
        "scope": "|".join(module.SCOPES),
    }).to_csv(score_path, index=False)
    table = module.read_baseline_pathway_scores(
        score_path, expected_ids=selected.obs_names.tolist(),
        expected_hash=record["selected_unit_ids_sha256"], side="raw",
    )
    assert table["unit_id"].tolist() == selected.obs_names.tolist()

    payload = {
        "coverage": {
            "requested_genes": ["a", "missing"],
            "present_genes": ["a"],
        },
        "metadata": _adapter_metadata(module),
        "scoring_scopes": list(module.SCOPES),
        "n_eligible_units": record["eligible_n_units"],
        "n_scored_units": record["selected_n_units"],
        "gene_axis": sample.raw.var_names.tolist(),
        "expression": sample.expression("raw"),
    }
    current = _adapter_metadata(module)
    current["rank_cutoff"] = 1
    current["memory_adapter"] = dict(current["memory_adapter"], n_genes=2,
                                      permutation_sha256="d" * 64)
    module.validate_baseline_pathway_metadata(
        payload, sample=sample, side="raw", record=record,
        emt_genes=["a", "missing"], current_provider_metadata=current,
    )
    payload["metadata"]["seed"] = 43
    with pytest.raises(RuntimeError, match="method metadata drift"):
        module.validate_baseline_pathway_metadata(
            payload, sample=sample, side="raw", record=record,
            emt_genes=["a", "missing"], current_provider_metadata=current,
        )

    broken = pd.read_csv(score_path).drop(columns="unit_id")
    broken.to_csv(score_path, index=False)
    with pytest.raises(ValueError, match="lacks columns"):
        module.read_baseline_pathway_scores(
            score_path, expected_ids=selected.obs_names.tolist(),
            expected_hash=record["selected_unit_ids_sha256"], side="raw",
        )


def test_input_and_registered_artifact_provenance_guards(tmp_path):
    module = diagnostic()
    output_root = tmp_path / "formal"
    relative = "P2CRC_Xenium/reconstruction_impact/tables/cohort_manifest.json"
    artifact = output_root / relative
    artifact.parent.mkdir(parents=True)
    artifact.write_text('{"schema_version": 1}', encoding="utf-8")
    recorded = module.file_identity(artifact)
    recorded["relative_to_output"] = relative
    index = {relative: recorded}

    module.require_registered_artifact(
        output_root=output_root, relative=relative, artifact_index=index,
        label="cohort_manifest",
    )
    artifact.write_text('{"schema_version": 2}', encoding="utf-8")
    with pytest.raises(RuntimeError, match="identity drift"):
        module.require_registered_artifact(
            output_root=output_root, relative=relative, artifact_index=index,
            label="cohort_manifest",
        )

    sample = tmp_path / "sample.yaml"
    sample.write_text("schema_version: 1\n", encoding="utf-8")
    identity = module.file_identity(sample)
    sample.write_text("schema_version: 2\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="identity drift"):
        module.require_file_identity(sample, identity, "sample_config")
    with pytest.raises(ValueError, match="schema_version"):
        module._registered_artifact_index({"schema_version": 2, "artifacts": []})


def test_final_reverification_rejects_post_preflight_artifact_drift(tmp_path, monkeypatch):
    module = diagnostic()
    analysis = tmp_path / "P2CRC_Xenium/reconstruction_impact"
    run = {
        "repo": tmp_path,
        "config_path": tmp_path / "project.yaml",
        "sample_yaml": tmp_path / "sample.yaml",
        "raw": tmp_path / "raw.h5ad",
        "svc": tmp_path / "SVC.h5ad",
        "geneset_path": tmp_path / "genes.gmt",
        "analysis_dir": analysis,
        "result_path": analysis / "result.json",
        "cohort_manifest_path": analysis / "tables/cohort_manifest.json",
    }
    input_paths = {
        "project_config": run["config_path"], "sample_config": run["sample_yaml"],
        "raw": run["raw"], "svc": run["svc"], "geneset": run["geneset_path"],
    }
    for path in input_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixed", encoding="utf-8")
    artifact_paths = {
        "formal_result": run["result_path"],
        "cohort_manifest": run["cohort_manifest_path"],
        **{
            f"baseline_moran_{scope}": analysis / f"tables/moran_shared_genes_{scope}.csv"
            for scope in module.SCOPES
        },
        **{
            f"baseline_pathway_scores_{side}": analysis / f"tables/pathway_{side}_{module.EMT_NAME}.csv"
            for side in module.SIDES
        },
        **{
            f"baseline_pathway_metadata_{side}": analysis / f"tables/pathway_{side}_{module.EMT_NAME}_metadata.json"
            for side in module.SIDES
        },
    }
    for path in artifact_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixed", encoding="utf-8")
    repository = {"commit": "current", "branch": "sl", "dirty": False, "dirty_paths": []}
    provenance = {
        "input_identities": {name: module.file_identity(path) for name, path in input_paths.items()},
        "artifact_identities": {name: module.file_identity(path) for name, path in artifact_paths.items()},
        "repository": repository,
        "receipt": {"repository": {"commit": "baseline"}},
        "code_identities": {"control": "fixed"},
    }
    monkeypatch.setattr(module, "git_identity", lambda _repo: repository)
    monkeypatch.setattr(module, "code_identities", lambda *_args: {"control": "fixed"})
    module.reverify_provenance(run, provenance)
    artifact_paths["baseline_moran_All"].write_text("drift", encoding="utf-8")
    with pytest.raises(RuntimeError, match="identity drift"):
        module.reverify_provenance(run, provenance)


def test_existing_output_is_refused_and_manifest_uses_final_paths(tmp_path):
    module = diagnostic()
    output = tmp_path / "moran-common-axis-001"
    output.mkdir()
    with pytest.raises(FileExistsError, match="already exists"):
        module.require_output_absent(output)

    staging = tmp_path / ".staging"
    staging.mkdir()
    (staging / "result.json").write_text("{}", encoding="utf-8")
    final = tmp_path / "final"
    manifest = module.build_artifact_manifest(staging, final)
    assert manifest["artifacts"][0]["path"] == str(final / "result.json")
    assert manifest["artifacts"][0]["relative_to_diagnostic"] == "result.json"
