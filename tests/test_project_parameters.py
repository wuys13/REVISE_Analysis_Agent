from pathlib import Path

import pytest
import yaml

from revise_analysis.runner import resolve_analysis_parameters


def _write_yaml(path: Path, value: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(value))
    return path


def _sample(path: Path, parameters: dict | None = None) -> Path:
    config = {
        "schema_version": 1,
        "sample_id": "sample",
        "files": {"raw": "raw.h5ad", "svc": "SVC.h5ad"},
        "expression": {
            side: {"identity": "known"} for side in ("raw", "svc")
        },
    }
    if parameters is not None:
        config["analysis_parameters"] = {"reconstruction_impact": parameters}
    return _write_yaml(path, config)


def test_parameter_precedence_is_sample_then_project_then_override(tmp_path):
    sample = _sample(tmp_path / "samples" / "sample.yaml", {
        "sample_n_units": 10,
        "random_state": 1,
    })
    project = _write_yaml(tmp_path / "project" / "project.yaml", {
        "schema_version": 1,
        "samples": ["../samples/sample.yaml"],
        "analyses": {"reconstruction_impact": {
            "sample_n_units": 20,
            "moran_n_neighbors": 4,
        }},
    })
    result = resolve_analysis_parameters(
        sample,
        "reconstruction_impact",
        project_yaml=project,
        overrides={"sample_n_units": 30},
    )
    assert result == {
        "sample_n_units": 30,
        "random_state": 1,
        "moran_n_neighbors": 4,
    }


def test_relative_resource_path_retains_its_declaring_yaml_base(tmp_path):
    sample = _sample(tmp_path / "samples" / "sample.yaml", {
        "geneset_path": "sample_sets.gmt",
        "random_state": 1,
    })
    project = _write_yaml(tmp_path / "projects" / "project.yaml", {
        "schema_version": 1,
        "samples": ["../samples/sample.yaml"],
        "analyses": {"reconstruction_impact": {"sample_n_units": 20}},
    })
    inherited = resolve_analysis_parameters(sample, "reconstruction_impact", project_yaml=project)
    assert inherited["geneset_path"] == str((sample.parent / "sample_sets.gmt").resolve())

    project_config = yaml.safe_load(project.read_text())
    project_config["analyses"]["reconstruction_impact"]["geneset_path"] = "resources/project_sets.gmt"
    project.write_text(yaml.safe_dump(project_config))
    overridden = resolve_analysis_parameters(sample, "reconstruction_impact", project_yaml=project)
    assert overridden["geneset_path"] == str((project.parent / "resources/project_sets.gmt").resolve())


def test_project_must_declare_sample_and_parameter_layers_must_be_mappings(tmp_path):
    sample = _sample(tmp_path / "sample.yaml")
    other = _sample(tmp_path / "other.yaml")
    project = _write_yaml(tmp_path / "project.yaml", {
        "schema_version": 1,
        "samples": [other.name],
        "analyses": {"reconstruction_impact": {}},
    })
    with pytest.raises(ValueError, match="not declared by project"):
        resolve_analysis_parameters(sample, "reconstruction_impact", project_yaml=project)

    config = yaml.safe_load(sample.read_text())
    config["analysis_parameters"] = ["not", "a", "mapping"]
    sample.write_text(yaml.safe_dump(config))
    with pytest.raises(ValueError, match="analysis_parameters must be a mapping"):
        resolve_analysis_parameters(sample, "reconstruction_impact")

    with pytest.raises(ValueError, match="parameters must be a mapping"):
        resolve_analysis_parameters(other, "reconstruction_impact", overrides=["bad"])
