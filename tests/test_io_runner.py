from pathlib import Path
from types import SimpleNamespace
import json

import numpy as np
import pandas as pd
import pytest
import yaml
from anndata import AnnData

from revise_analysis.io import load_sample, read_gene_sets
from revise_analysis.runner import run_analysis, resolve_parameters
from revise_analysis.batch import run_batch


@pytest.fixture
def sample_yaml(tmp_path):
    raw = AnnData(np.array([[1., 0], [2, 1], [0, 2]]), obs=pd.DataFrame({"Level1": ["Mono/Macro", "T", "T"]}, index=["a", "b", "c"]), var=pd.DataFrame(index=["G1", "G2"]))
    svc = AnnData(np.array([[0.5, 3], [1.5, 1]]), obs=pd.DataFrame({"Level1": ["Mono_Macro", "T"]}, index=["z", "b"]), var=pd.DataFrame(index=["G2", "G1"]))
    raw.write_h5ad(tmp_path / "raw.h5ad")
    svc.write_h5ad(tmp_path / "SVC.h5ad")
    config = {"schema_version": 1, "sample_id": "sample", "files": {"raw": "raw.h5ad", "svc": "SVC.h5ad"}, "expression": {"scale": "untransformed_nonnegative"}}
    path = tmp_path / "sample.yaml"
    path.write_text(yaml.safe_dump(config))
    return path


def test_loader_preserves_native_axes_and_values(sample_yaml):
    sample = load_sample(sample_yaml)
    assert sample.raw.obs_names.tolist() == ["a", "b", "c"]
    assert sample.svc.obs_names.tolist() == ["z", "b"]
    assert sample.svc.var_names.tolist() == ["G2", "G1"]
    np.testing.assert_array_equal(sample.svc.X, [[0.5, 3], [1.5, 1]])
    assert sample.labels("svc").iloc[0] == "Mono_Macro"
    assert sample.svc.obs.Level1.iloc[0] == "Mono_Macro"


def test_unknown_input_scale_is_not_guessed(sample_yaml):
    config = yaml.safe_load(sample_yaml.read_text())
    config["expression"] = {
        "scale": "unknown",
        "raw": {"identity": "known"},
        "svc": {"identity": "known"},
    }
    sample_yaml.write_text(yaml.safe_dump(config))
    sample = load_sample(sample_yaml)
    assert sample.expression("svc")["scale"] == "untransformed_nonnegative"
    assert "explicitly unknown" in sample.expression_unavailable("svc")


def test_failure_cannot_publish_stale_success(sample_yaml, monkeypatch, tmp_path):
    import revise_analysis.runner as runner
    def success(sample, path, parameters):
        (path / "table.csv").write_text("value\n7\n")
        return {"status": "succeeded", "parameters": parameters, "outputs": {"table": "table.csv"}}
    module = SimpleNamespace(run=success)
    original_import = runner.importlib.import_module
    monkeypatch.setattr(runner.importlib, "import_module", lambda name, *args, **kw: module if name.startswith("revise_analysis.analyses.") else original_import(name, *args, **kw))
    output = tmp_path / "output"
    first = run_analysis(sample_yaml, "spatial_autocorrelation", output)
    assert first["status"] == "succeeded"
    destination = output / "sample" / "spatial_autocorrelation"
    (destination / "manual-note.txt").write_text("preserve me")
    def failure(*args):
        raise RuntimeError("new calculation failed")
    module.run = failure
    result = run_analysis(sample_yaml, "spatial_autocorrelation", output)
    assert result["status"] == "failed"
    assert not (destination / "table.csv").exists()
    assert json.loads((destination / "result.json").read_text())["error"]["message"] == "new calculation failed"
    assert list((output / "sample" / ".previous").glob("*/manual-note.txt"))
    assert json.loads((output / "sample" / "index.json").read_text())["analyses"]["spatial_autocorrelation"] == "spatial_autocorrelation/result.json"


def test_unowned_result_directory_is_protected(sample_yaml, tmp_path):
    output = tmp_path / "output"
    destination = output / "sample" / "spatial_autocorrelation"
    destination.mkdir(parents=True)
    (destination / "manual.txt").write_text("mine")
    with pytest.raises(FileExistsError):
        run_analysis(sample_yaml, "spatial_autocorrelation", output)
    assert (destination / "manual.txt").read_text() == "mine"


def test_batch_continues_after_bad_sample(sample_yaml, tmp_path, monkeypatch):
    import revise_analysis.batch as batch
    calls = []
    def fake_run(sample, analysis, output, params):
        calls.append((sample, analysis))
        return {"sample_id": "sample", "status": "succeeded"}
    monkeypatch.setattr(batch, "run_analysis", fake_run)
    project = tmp_path / "project.yaml"
    project.write_text(yaml.safe_dump({"schema_version": 1, "samples": ["missing.yaml", sample_yaml.name], "analyses": {"spatial_autocorrelation": {}, "pathway_activity": {}}, "output_dir": "results"}))
    result = run_batch(project)
    assert [task["status"] for task in result["tasks"]] == ["failed", "failed", "succeeded", "succeeded"]
    assert len(calls) == 2


def test_resources_resolve_from_declaring_yaml(tmp_path):
    (tmp_path / "sets.gmt").write_text("A\tdescription\tG1\tG2\nB\tdescription\tG3\n")
    result = resolve_parameters({"geneset_path": "sets.gmt", "gene_set_names": ["A"]}, tmp_path)
    assert result == {"geneset_path": str(tmp_path / "sets.gmt"), "gene_set_names": ["A"]}
    assert read_gene_sets(result["geneset_path"], result["gene_set_names"]) == {"A": ["G1", "G2"]}


def test_traversal_sample_id_rejected(sample_yaml):
    config = yaml.safe_load(sample_yaml.read_text())
    config["sample_id"] = "../elsewhere"
    sample_yaml.write_text(yaml.safe_dump(config))
    with pytest.raises(ValueError, match="Unsafe"):
        load_sample(sample_yaml)


def test_invalid_declaration_rerun_replaces_old_success(sample_yaml, monkeypatch, tmp_path):
    import revise_analysis.runner as runner
    original = runner.importlib.import_module
    fake = SimpleNamespace(run=lambda *args: {"status": "succeeded", "outputs": {}, "parameters": {}})
    monkeypatch.setattr(runner.importlib, "import_module", lambda name, *a, **kw: fake if name.startswith("revise_analysis.analyses.") else original(name, *a, **kw))
    output = tmp_path / "output"
    run_analysis(sample_yaml, "spatial_autocorrelation", output)
    config = yaml.safe_load(sample_yaml.read_text())
    config["expression"]["svc"] = "not a mapping"
    sample_yaml.write_text(yaml.safe_dump(config))
    project = tmp_path / "project.yaml"
    project.write_text(yaml.safe_dump({"schema_version": 1, "samples": [sample_yaml.name], "analyses": {"spatial_autocorrelation": {}}, "output_dir": "output"}))
    assert run_batch(project)["tasks"][0]["status"] == "failed"
    current = json.loads((output / "sample/spatial_autocorrelation/result.json").read_text())
    assert current["status"] == "failed"


def test_corrupt_generated_index_is_rebuilt(sample_yaml, monkeypatch, tmp_path):
    import revise_analysis.runner as runner
    original = runner.importlib.import_module
    fake = SimpleNamespace(run=lambda *args: {"status": "succeeded", "outputs": {}, "parameters": {}})
    monkeypatch.setattr(runner.importlib, "import_module", lambda name, *a, **kw: fake if name.startswith("revise_analysis.analyses.") else original(name, *a, **kw))
    output = tmp_path / "output"
    run_analysis(sample_yaml, "spatial_autocorrelation", output)
    (output / "sample/index.json").write_text("broken")
    result = run_analysis(sample_yaml, "spatial_autocorrelation", output)
    assert result["status"] == "succeeded"
    assert json.loads((output / "sample/index.json").read_text())["sample_id"] == "sample"


def test_per_side_expression_contract_and_label_only_input(sample_yaml):
    import anndata as ad
    config = yaml.safe_load(sample_yaml.read_text())
    config['expression'] = {
        'raw': {'identity': 'measured_counts', 'scale': 'untransformed_nonnegative'},
        'svc': {'identity': 'unknown', 'scale': 'unknown'},
    }
    svc = ad.read_h5ad(sample_yaml.parent / 'SVC.h5ad')
    svc.X = None
    svc.obs['SVC_cluster'] = ['r0', 'r1']
    svc.write_h5ad(sample_yaml.parent / 'SVC.h5ad')
    sample_yaml.write_text(yaml.safe_dump(config))
    sample = load_sample(sample_yaml)
    assert sample.expression_unavailable('raw') is None
    assert 'unknown' in sample.expression_unavailable('svc')
    assert sample.labels('svc', sample.reconstruction_key).tolist() == ['r0', 'r1']
    assert sample.svc.X is None


def test_gene_set_resource_rejects_ambiguous_rows(tmp_path):
    resource = tmp_path / 'bad.gmt'
    resource.write_text('A\tdesc\tG1\nA\tdesc\tG2\n')
    with pytest.raises(ValueError, match='Duplicate'):
        read_gene_sets(resource)
    resource.write_text('A\tdesc\n')
    with pytest.raises(ValueError, match='Invalid GMT row'):
        read_gene_sets(resource)
