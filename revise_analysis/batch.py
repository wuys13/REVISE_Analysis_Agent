"""Small sequential sample × analysis execution, without reconstruction coupling."""
from pathlib import Path

from .io import read_yaml, safe_segment, resolve_path
from .runner import ANALYSES, resolve_analysis_parameters, run_analysis, write_json


def run_batch(project_yaml) -> dict:
    source = Path(project_yaml).expanduser().resolve()
    config = read_yaml(source)
    samples = config.get("samples")
    analyses = config.get("analyses")
    if not isinstance(samples, list) or not samples:
        raise ValueError("samples must be a nonempty list of sample YAML paths")
    if not isinstance(analyses, dict) or not analyses or set(analyses) - set(ANALYSES):
        raise ValueError(f"analyses must map supported names to parameter mappings: {ANALYSES}")
    for name, parameters in analyses.items():
        if parameters is not None and not isinstance(parameters, dict):
            raise ValueError(f"{name}: parameters must be a mapping")
    output = resolve_path(config.get("output_dir", "output"), source.parent)
    output.mkdir(parents=True, exist_ok=True)
    paths = [resolve_path(path, source.parent) for path in samples]
    # Duplicate sample identities would silently overwrite results, so reject them.
    seen = {}
    errors = {}
    for path in paths:
        try:
            sample_config = read_yaml(path)
            safe_segment(sample_config.get("sample_id"))
        except Exception as exc:
            errors[path] = f"{type(exc).__name__}: {exc}"
            continue
        identifier = sample_config["sample_id"]
        if identifier in seen:
            raise ValueError(f"Duplicate sample_id {identifier!r}: {seen[identifier]} and {path}")
        seen[identifier] = path
    tasks = []
    for sample in paths:
        for analysis, parameters in analyses.items():
            try:
                if sample in errors:
                    raise ValueError(errors[sample])
                effective = resolve_analysis_parameters(
                    sample, analysis, project_yaml=source,
                )
                result = run_analysis(sample, analysis, output, effective)
                tasks.append({"sample": str(sample), "sample_id": result["sample_id"], "analysis": analysis,
                              "status": result["status"], "result": f"{result['sample_id']}/{analysis}/result.json"})
            except Exception as exc:
                tasks.append({"sample": str(sample), "analysis": analysis, "status": "failed",
                              "error": f"{type(exc).__name__}: {exc}"})
    result = {"schema_version": 1, "status": "succeeded" if all(task["status"] == "succeeded" for task in tasks) else "incomplete", "tasks": tasks}
    write_json(output / "batch_result.json", result)
    return result
