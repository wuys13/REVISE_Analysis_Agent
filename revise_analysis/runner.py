"""Execute one analysis and publish only its current, saved results."""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import importlib
import json
import math
from pathlib import Path
import shutil
import tempfile
import traceback
from uuid import uuid4

from .io import load_sample, read_yaml, safe_segment, resolve_path

ANALYSES = ("reconstruction_impact", "spatial_autocorrelation", "pathway_activity")
STATUSES = {"succeeded", "partial", "skipped", "failed"}


def json_value(value):
    if isinstance(value, dict):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_value(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if hasattr(value, "item"):
        return json_value(value.item())
    return value


def write_json(path: Path, value: dict):
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    temporary.write_text(json.dumps(json_value(value), ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


@contextmanager
def _sample_lock(sample_dir: Path):
    with (sample_dir / ".run.lock").open("a") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"Another analysis is publishing to {sample_dir}") from exc
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def resolve_parameters(parameters: dict | None, base: Path) -> dict:
    """Resolve resource paths; each analysis checks only its own requirements."""
    effective = _resolve_parameter_layer(parameters, base)
    resource = effective.get("geneset_path")
    names = effective.get("gene_set_names")
    if resource is not None:
        if effective.get("gene_sets") is not None:
            raise ValueError("Use either gene_sets or geneset_path, not both")
    elif names is not None:
        raise ValueError("gene_set_names requires geneset_path")
    return effective


def _resolve_parameter_layer(parameters: dict | None, base: Path) -> dict:
    """Resolve paths in one precedence layer before it is overlaid."""
    if parameters is not None and not isinstance(parameters, dict):
        raise ValueError("analysis parameters must be a mapping")
    effective = dict(parameters or {})
    resource = effective.get("geneset_path")
    if resource is not None:
        effective["geneset_path"] = str(resolve_path(resource, base))
    return effective


def _analysis_parameters(config: dict, key: str, analysis: str) -> dict:
    configured = config.get(key, {})
    if configured is None:
        configured = {}
    if not isinstance(configured, dict):
        raise ValueError(f"{key} must be a mapping")
    parameters = configured.get(analysis, {})
    if parameters is None:
        parameters = {}
    if not isinstance(parameters, dict):
        raise ValueError(f"{key}.{analysis} must be a mapping")
    return parameters


def resolve_analysis_parameters(
    sample_yaml,
    analysis: str,
    *,
    project_yaml=None,
    overrides: dict | None = None,
) -> dict:
    """Merge explicit sample, project and call parameters with source-aware paths.

    Analysis implementations remain the owners of method defaults. Relative
    resource paths are made absolute before layers merge, so a surviving value
    always retains the base directory of the YAML that declared it. Explicit
    call overrides follow the direct-run convention and resolve from the sample
    YAML directory.
    """
    if analysis not in ANALYSES:
        raise ValueError(f"Unknown analysis {analysis!r}; choose from {ANALYSES}")
    sample_source = Path(sample_yaml).expanduser().resolve()
    sample_config = read_yaml(sample_source)
    safe_segment(sample_config.get("sample_id"))
    sample_parameters = _resolve_parameter_layer(
        _analysis_parameters(sample_config, "analysis_parameters", analysis),
        sample_source.parent,
    )

    project_parameters = {}
    if project_yaml is not None:
        project_source = Path(project_yaml).expanduser().resolve()
        project_config = read_yaml(project_source)
        samples = project_config.get("samples")
        if (not isinstance(samples, list) or not samples
                or any(not isinstance(value, (str, Path)) for value in samples)):
            raise ValueError("samples must be a nonempty list of sample YAML paths")
        project_samples = {resolve_path(value, project_source.parent) for value in samples}
        if sample_source not in project_samples:
            raise ValueError(f"Sample {sample_source} is not declared by project {project_source}")
        project_parameters = _resolve_parameter_layer(
            _analysis_parameters(project_config, "analyses", analysis),
            project_source.parent,
        )

    override_parameters = _resolve_parameter_layer(overrides, sample_source.parent)
    merged = {**sample_parameters, **project_parameters, **override_parameters}
    # All resource paths are absolute by now; this validates the merged choice
    # without changing which declaration source owns a surviving path.
    return resolve_parameters(merged, sample_source.parent)


def _check_destination(destination: Path, sample_id: str, analysis: str):
    if not destination.exists():
        return
    record = destination / "result.json"
    try:
        previous = json.loads(record.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise FileExistsError(f"Refusing to replace a directory not owned by this analysis: {destination}") from exc
    if previous.get("sample_id") != sample_id or previous.get("analysis") != analysis:
        raise FileExistsError(f"Existing result identity differs: {destination}")


def _publish(staging: Path, destination: Path):
    """Keep prior results separate so failures and user annotations remain honest."""
    previous = None
    if destination.exists():
        archive = destination.parent / ".previous"
        archive.mkdir(exist_ok=True)
        previous = archive / f"{destination.name}-{uuid4().hex[:12]}"
        destination.rename(previous)
    try:
        staging.rename(destination)
    except BaseException:
        if previous is not None:
            previous.rename(destination)
        raise


def run_analysis(sample_yaml, analysis: str, output_dir, parameters: dict | None = None) -> dict:
    """Run one registered analysis into output/<sample>/<analysis>.

    Invalid identities/configuration raise. Runtime errors publish a current failed
    record; batch callers can then proceed to the next task. No input is modified.
    """
    if analysis not in ANALYSES:
        raise ValueError(f"Unknown analysis {analysis!r}; choose from {ANALYSES}")
    source = Path(sample_yaml).expanduser().resolve()
    config = read_yaml(source)
    sample_id = safe_segment(config.get("sample_id"))
    sample_dir = Path(output_dir).expanduser().resolve() / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    destination = sample_dir / analysis
    with _sample_lock(sample_dir):
        _check_destination(destination, sample_id, analysis)
        staging = Path(tempfile.mkdtemp(prefix=f".{analysis}-", dir=sample_dir))
        result = {
            "schema_version": 1, "sample_id": sample_id, "analysis": analysis,
            "started_at": datetime.now(timezone.utc).isoformat(),
            "parameters": dict(parameters or {}), "outputs": {}, "unavailable": [],
        }
        try:
            try:
                effective = resolve_analysis_parameters(source, analysis, overrides=parameters)
                result["parameters"] = effective
                sample = load_sample(source)
                implementation = importlib.import_module(f"revise_analysis.analyses.{analysis}")
                if analysis == "reconstruction_impact":
                    computed = implementation.run(sample, staging, effective, continue_on_error=True)
                else:
                    computed = implementation.run(sample, staging, effective)
                if computed.get("status") not in STATUSES:
                    raise ValueError("Analysis must return an explicit supported status")
                result.update(computed)
                # Identity belongs to the runner, never to arbitrary analysis output.
                result.update(sample_id=sample_id, analysis=analysis, schema_version=1)
            except Exception as exc:
                # Remove only this fresh temporary run, never a published/user directory.
                for child in staging.iterdir():
                    shutil.rmtree(child) if child.is_dir() and not child.is_symlink() else child.unlink()
                (staging / "error.txt").write_text(traceback.format_exc(), encoding="utf-8")
                result.update(status="failed", outputs={"error": "error.txt"}, unavailable=[],
                              error={"type": type(exc).__name__, "message": str(exc)})
            result["finished_at"] = datetime.now(timezone.utc).isoformat()
            write_json(staging / "result.json", result)
            try:
                from .reporting import render_report
                render_report(staging)
                result["outputs"]["report"] = "report.html"
            except Exception as exc:
                result.setdefault("unavailable", []).append({"component": "report", "reason": f"{type(exc).__name__}: {exc}"})
                if result["status"] == "succeeded":
                    result["status"] = "partial"
            write_json(staging / "result.json", result)
            index_path = sample_dir / "index.json"
            # Index is generated navigation: rebuild from actual saved records.
            index = {"sample_id": sample_id, "analyses": {}}
            for record in sample_dir.glob("*/result.json"):
                if record.parent.name.startswith(".") or record.parent == destination:
                    continue
                index["analyses"][record.parent.name] = f"{record.parent.name}/result.json"
            _publish(staging, destination)
            index["analyses"][analysis] = f"{analysis}/result.json"
            write_json(index_path, index)
            return json_value(result)
        finally:
            if staging.exists():
                shutil.rmtree(staging)
