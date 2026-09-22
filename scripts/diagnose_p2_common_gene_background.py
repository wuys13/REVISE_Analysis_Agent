#!/usr/bin/env python3
"""Run the P2 common-gene-background diagnostic after the formal Impact run.

This is a run-specific, read-only consumer of the three formal sample inputs.
It reproduces every native Moran and pathway cohort from the formal cohort
manifest, restricts both sides to one ordered common-gene axis, and then calls
the existing Moran and AUCell methods.  It never modifies the Impact result or
the input H5AD files.
"""
from __future__ import annotations

# ruff: noqa: E402

import argparse
from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import platform
import resource
import socket
import subprocess
import sys
import tempfile
import traceback
from typing import Any, Callable

# Direct ``python scripts/...py`` execution must import this exact checkout,
# even when the package has not been installed into the selected interpreter.
SCRIPT_REPO = Path(__file__).resolve().parents[1]
if str(SCRIPT_REPO) not in sys.path:
    sys.path.insert(0, str(SCRIPT_REPO))

import numpy as np
import pandas as pd
from scipy import sparse

from revise_analysis.analyses._shared import (
    UNIT_ID_HASH_CONTRACT,
    deterministic_subset,
    unit_id_digest,
)
from revise_analysis.io import load_sample, read_yaml, resolve_path
from revise_analysis.methods.errors import PrerequisiteUnavailable
from revise_analysis.methods.spatial import compute_moran
from revise_analysis.analyses.pathway_activity import (
    compute_pathway_scores,
    resolve_gene_sets,
)


TASK_ID = "p2-xenium-impact-20260922"
ITERATION = "iter-001"
DIAGNOSTIC_ID = "common-gene-background-001"
ANALYSIS = "reconstruction_impact"
SAMPLE_ID = "P2CRC_Xenium"
SCOPES = ("All", "Fibroblast", "Mono_Macro", "T")
SIDES = ("raw", "svc")
EXPECTED_COMMON_GENES = 343
EXPECTED_SAMPLE_N_UNITS = 30_000
EXPECTED_RANDOM_STATE = 42
EXPECTED_MORAN_NEIGHBORS = 6
EXPECTED_PATHWAY_AUC_THRESHOLD = 0.01
EMT_NAME = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
GENE_AXIS_HASH_CONTRACT = "sha256-json-utf8-ordered-v1"
DEFAULT_CONFIG = Path("configs/p2_xenium_impact_20260922_iter_001.yaml")
CONTROL_CODE_PATHS = (
    "revise_analysis/methods/spatial.py",
    "revise_analysis/analyses/_shared.py",
    "revise_analysis/analyses/impact_molecular.py",
    "revise_analysis/analyses/pathway_activity.py",
    "revise_analysis/methods/programs.py",
    "revise_analysis/methods/omicverse_aucell_adapter.py",
    "revise_analysis/io.py",
)
PACKAGE_NAMES = (
    "revise-analysis-agent", "anndata", "numpy", "pandas", "scipy",
    "scikit-learn", "scanpy", "omicverse", "ctxcore", "gseapy",
)


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    path = path.resolve()
    stat = path.stat()
    return {"path": str(path), "size_bytes": stat.st_size, "sha256": sha256(path)}


def ordered_digest(values: Any) -> str:
    payload = json.dumps([str(value) for value in values], ensure_ascii=False,
                         separators=(",", ":")).encode("utf-8")
    return sha256_bytes(payload)


def read_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False),
                    encoding="utf-8")


def environment_identity() -> dict[str, Any]:
    packages = {}
    for name in PACKAGE_NAMES:
        try:
            packages[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            packages[name] = None
    return {
        "hostname": socket.gethostname(),
        "python_executable": sys.executable,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "packages": packages,
    }


def git_identity(repo: Path) -> dict[str, Any]:
    def run(*args: str) -> str:
        return subprocess.run(("git", "-C", str(repo), *args), check=True,
                              text=True, capture_output=True).stdout.strip()

    dirty_paths = run("status", "--porcelain").splitlines()
    return {
        "commit": run("rev-parse", "HEAD"),
        "branch": run("branch", "--show-current"),
        "dirty": bool(dirty_paths),
        "dirty_paths": dirty_paths,
    }


def require_clean_worktree(identity: dict[str, Any]) -> None:
    if identity.get("dirty"):
        raise RuntimeError("diagnostic execution requires a clean committed worktree")


def require_output_absent(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"diagnostic output already exists: {path}")


def require_file_identity(path: Path, recorded: dict[str, Any], label: str) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"{label} is absent: {path}")
    actual = file_identity(path)
    recorded_sha = recorded.get("sha256")
    if (not isinstance(recorded.get("path"), str)
            or Path(recorded["path"]).resolve() != path.resolve()):
        raise RuntimeError(f"{label} recorded path drift")
    if (not isinstance(recorded.get("size_bytes"), int)
            or not isinstance(recorded_sha, str) or len(recorded_sha) != 64
            or any(char not in "0123456789abcdef" for char in recorded_sha)):
        raise ValueError(f"{label} has an invalid recorded file identity")
    for key in ("size_bytes", "sha256"):
        if recorded.get(key) != actual[key]:
            raise RuntimeError(
                f"{label} identity drift: recorded {key}={recorded.get(key)!r}, "
                f"actual {key}={actual[key]!r}"
            )
    return actual


def code_identities(repo: Path, baseline_commit: str) -> dict[str, Any]:
    values = {}
    for relative in CONTROL_CODE_PATHS:
        current_path = repo / relative
        current = file_identity(current_path)
        completed = subprocess.run(
            ("git", "-C", str(repo), "show", f"{baseline_commit}:{relative}"),
            check=True, capture_output=True,
        )
        baseline_sha = sha256_bytes(completed.stdout)
        if current["sha256"] != baseline_sha:
            raise RuntimeError(f"scientific control code drift since formal Impact: {relative}")
        values[relative] = {
            "current_sha256": current["sha256"],
            "formal_impact_commit": baseline_commit,
            "formal_impact_blob_sha256": baseline_sha,
            "unchanged": True,
        }
    script = file_identity(Path(__file__))
    values[str(Path(__file__).resolve().relative_to(repo.resolve()))] = {
        "current_sha256": script["sha256"],
        "formal_impact_commit": None,
        "formal_impact_blob_sha256": None,
        "unchanged": None,
    }
    return values


def _registered_artifact_index(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    if manifest.get("schema_version") != 1:
        raise ValueError("formal artifact manifest schema_version must be 1")
    if manifest.get("missing_registered_artifacts"):
        raise RuntimeError("formal Impact artifact manifest reports missing registered artifacts")
    rows = manifest.get("artifacts")
    if not isinstance(rows, list):
        raise ValueError("formal artifact manifest has no artifacts list")
    index = {}
    for row in rows:
        relative = row.get("relative_to_output")
        if not isinstance(relative, str) or relative in index:
            raise ValueError("formal artifact manifest contains invalid or duplicate paths")
        index[relative] = row
    return index


def require_registered_artifact(
    *, output_root: Path, relative: str, artifact_index: dict[str, dict[str, Any]], label: str,
) -> dict[str, Any]:
    try:
        recorded = artifact_index[relative]
    except KeyError as exc:
        raise RuntimeError(f"{label} is not registered by the formal artifact manifest: {relative}") from exc
    path = (output_root / relative).resolve()
    path.relative_to(output_root.resolve())
    return require_file_identity(path, recorded, label)


def _stage_status(result: dict[str, Any], stage: str) -> str | None:
    matches = [row.get("status") for row in result.get("stages", []) if row.get("stage") == stage]
    if len(matches) != 1:
        raise ValueError(f"formal result must contain exactly one {stage!r} stage record")
    return matches[0]


def validate_parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    expected = {
        "scopes": list(SCOPES),
        "sample_n_units": EXPECTED_SAMPLE_N_UNITS,
        "random_state": EXPECTED_RANDOM_STATE,
        "moran_n_neighbors": EXPECTED_MORAN_NEIGHBORS,
        "pathway_auc_threshold": EXPECTED_PATHWAY_AUC_THRESHOLD,
        "gene_set_names": [EMT_NAME],
    }
    changed = [key for key, value in expected.items() if parameters.get(key) != value]
    if changed:
        raise ValueError(f"formal Moran parameters drifted: {changed}")
    return expected


def resolve_run(
    repo: Path, config_path: Path, *, baseline_run: Path | None = None,
    output: Path | None = None,
) -> dict[str, Any]:
    repo = repo.expanduser().resolve()
    if repo != SCRIPT_REPO:
        raise ValueError(f"--repo must identify the checkout containing this script: {SCRIPT_REPO}")
    config_path = config_path.expanduser().resolve()
    config = read_yaml(config_path)
    analyses = config.get("analyses")
    if not isinstance(analyses, dict) or list(analyses) != [ANALYSIS]:
        raise ValueError(f"project must contain only {ANALYSIS!r}")
    parameters = analyses[ANALYSIS]
    validate_parameters(parameters)
    samples = config.get("samples")
    if not isinstance(samples, list) or len(samples) != 1:
        raise ValueError("project must contain exactly one sample")
    output_root = resolve_path(config.get("output_dir"), config_path.parent)
    expected_output = (repo / "output" / TASK_ID / ITERATION).resolve()
    if output_root != expected_output:
        raise ValueError(f"formal output root must be {expected_output}")
    if baseline_run is not None and baseline_run.expanduser().resolve() != output_root:
        raise ValueError(f"--baseline-run must resolve to {output_root}")
    sample_yaml = resolve_path(samples[0], config_path.parent)
    sample_config = read_yaml(sample_yaml)
    if sample_config.get("sample_id") != SAMPLE_ID:
        raise ValueError(f"sample_id must be {SAMPLE_ID}")
    raw = resolve_path(sample_config["files"]["raw"], sample_yaml.parent)
    svc = resolve_path(sample_config["files"]["svc"], sample_yaml.parent)
    geneset_path = resolve_path(parameters["geneset_path"], config_path.parent)
    parameters = dict(parameters, geneset_path=str(geneset_path))
    analysis_dir = output_root / SAMPLE_ID / ANALYSIS
    diagnostic_dir = output_root / "diagnostics" / DIAGNOSTIC_ID
    if output is not None and output.expanduser().resolve() != diagnostic_dir:
        raise ValueError(f"--output must resolve to {diagnostic_dir}")
    return {
        "repo": repo.resolve(),
        "config_path": config_path,
        "config": config,
        "parameters": parameters,
        "output_root": output_root,
        "sample_yaml": sample_yaml,
        "raw": raw,
        "svc": svc,
        "geneset_path": geneset_path,
        "analysis_dir": analysis_dir,
        "result_path": analysis_dir / "result.json",
        "cohort_manifest_path": analysis_dir / "tables/cohort_manifest.json",
        "receipt_path": output_root / "_run/run_receipt.json",
        "artifact_manifest_path": output_root / "_run/artifact_hashes.json",
        "diagnostic_dir": diagnostic_dir,
    }


def verify_formal_provenance(run: dict[str, Any]) -> dict[str, Any]:
    receipt = read_json(run["receipt_path"])
    if receipt.get("schema_version") != 1:
        raise ValueError("formal Impact receipt schema_version must be 1")
    if receipt.get("status") != "execution_complete":
        raise RuntimeError("formal Impact receipt is not execution_complete")
    if receipt.get("task_id") != TASK_ID or receipt.get("iteration") != ITERATION:
        raise RuntimeError("formal Impact receipt belongs to another run")
    result = read_json(run["result_path"])
    if result.get("status") not in {"succeeded", "partial"}:
        raise RuntimeError("formal Impact result is not a completed result")
    validate_parameters(result.get("parameters", {}))
    if result.get("stage_errors"):
        raise RuntimeError("formal Impact result contains stage errors")
    if _stage_status(result, "input") != "completed":
        raise RuntimeError("formal Impact input stage is not completed")
    if _stage_status(result, "molecular") not in {"completed", "partial"}:
        raise RuntimeError("formal Impact molecular stage is not completed or partial")

    artifact_manifest = read_json(run["artifact_manifest_path"])
    artifact_index = _registered_artifact_index(artifact_manifest)
    sample_prefix = f"{SAMPLE_ID}/{ANALYSIS}"
    required_relatives = {
        "formal_result": f"{sample_prefix}/result.json",
        "cohort_manifest": f"{sample_prefix}/tables/cohort_manifest.json",
    }
    for scope in SCOPES:
        required_relatives[f"baseline_moran_{scope}"] = (
            f"{sample_prefix}/tables/moran_shared_genes_{scope}.csv"
        )
    for side in SIDES:
        stem = f"pathway_{side}_{EMT_NAME}"
        required_relatives[f"baseline_pathway_scores_{side}"] = (
            f"{sample_prefix}/tables/{stem}.csv"
        )
        required_relatives[f"baseline_pathway_metadata_{side}"] = (
            f"{sample_prefix}/tables/{stem}_metadata.json"
        )
    registered_outputs = result.get("outputs")
    if not isinstance(registered_outputs, dict):
        raise ValueError("formal Impact result has no outputs mapping")
    registered_output_paths = set(registered_outputs.values())
    if not all(isinstance(path, str) for path in registered_output_paths):
        raise ValueError("formal Impact result contains non-string output paths")
    missing_outputs = [
        relative.removeprefix(f"{sample_prefix}/")
        for relative in required_relatives.values()
        if relative.startswith(f"{sample_prefix}/")
        and relative.removeprefix(f"{sample_prefix}/") != "result.json"
        and relative.removeprefix(f"{sample_prefix}/") not in registered_output_paths
    ]
    if missing_outputs:
        raise RuntimeError(f"formal Impact result does not register required outputs: {missing_outputs}")
    artifact_identities = {
        label: require_registered_artifact(
            output_root=run["output_root"], relative=relative,
            artifact_index=artifact_index, label=label,
        )
        for label, relative in required_relatives.items()
    }

    input_identity = receipt.get("input_identity")
    if not isinstance(input_identity, dict):
        raise ValueError("formal receipt has no input_identity mapping")
    current_inputs = {
        "project_config": run["config_path"],
        "sample_config": run["sample_yaml"],
        "raw": run["raw"],
        "svc": run["svc"],
        "geneset": run["geneset_path"],
    }
    verified_inputs = {}
    for name, path in current_inputs.items():
        if name not in input_identity:
            raise RuntimeError(f"formal receipt does not identify {name}")
        verified_inputs[name] = require_file_identity(path, input_identity[name], name)

    repository = receipt.get("repository")
    if not isinstance(repository, dict) or not repository.get("commit"):
        raise ValueError("formal receipt has no repository commit")
    if repository.get("dirty"):
        raise RuntimeError("formal Impact receipt records a dirty worktree")
    code = code_identities(run["repo"], repository["commit"])
    return {
        "receipt": receipt,
        "result": result,
        "artifact_manifest": artifact_manifest,
        "artifact_identities": artifact_identities,
        "input_identities": verified_inputs,
        "code_identities": code,
    }


def reverify_provenance(run: dict[str, Any], provenance: dict[str, Any]) -> None:
    for name, path in {
        "project_config": run["config_path"],
        "sample_config": run["sample_yaml"],
        "raw": run["raw"],
        "svc": run["svc"],
        "geneset": run["geneset_path"],
    }.items():
        require_file_identity(path, provenance["input_identities"][name], name)
    label_to_path = {
        "formal_result": run["result_path"],
        "cohort_manifest": run["cohort_manifest_path"],
        **{
            f"baseline_moran_{scope}": run["analysis_dir"] / f"tables/moran_shared_genes_{scope}.csv"
            for scope in SCOPES
        },
        **{
            f"baseline_pathway_scores_{side}": run["analysis_dir"] / f"tables/pathway_{side}_{EMT_NAME}.csv"
            for side in SIDES
        },
        **{
            f"baseline_pathway_metadata_{side}": run["analysis_dir"] / f"tables/pathway_{side}_{EMT_NAME}_metadata.json"
            for side in SIDES
        },
    }
    for label, path in label_to_path.items():
        require_file_identity(path, provenance["artifact_identities"][label], label)
    repository = git_identity(run["repo"])
    require_clean_worktree(repository)
    if repository["commit"] != provenance["repository"]["commit"]:
        raise RuntimeError("diagnostic repository commit drifted during execution")
    if code_identities(run["repo"], provenance["receipt"]["repository"]["commit"]) != provenance["code_identities"]:
        raise RuntimeError("diagnostic/control code drifted during execution")


def _scope(sample: Any, side: str, scope: str):
    data = getattr(sample, side)
    if scope == "All":
        return data
    if sample.broad_key not in data.obs:
        raise RuntimeError(f"{side}: broad-label column is unavailable for {scope}")
    labels = sample.labels(side)
    selected = labels.index[labels == scope]
    if selected.empty:
        raise RuntimeError(f"{side}: no units for scope {scope}")
    return data[selected]


def _baseline_cohort_record(
    manifest: dict[str, Any], side: str, scope: str, *, component: str = "moran",
) -> dict[str, Any]:
    matches = [
        row for row in manifest.get("records", [])
        if row.get("stage") == "molecular" and row.get("component") == component
        and row.get("side") == side and row.get("scope") == scope
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"formal cohort manifest needs one molecular {component} record for {side}/{scope}"
        )
    return matches[0]


def actual_cohort_record(
    *, side: str, scope: str, eligible: Any, selected: Any, seed: int,
    component: str = "moran", method_random_state: int | None = None,
) -> dict[str, Any]:
    method = ("all_eligible" if selected.n_obs == eligible.n_obs
              else "deterministic_without_replacement")
    return {
        "stage": "molecular",
        "component": component,
        "side": side,
        "scope": scope,
        "eligible_n_units": int(eligible.n_obs),
        "selected_n_units": int(selected.n_obs),
        "selection_method": method,
        "selection_random_state": seed,
        "selection_seed_applied": method == "deterministic_without_replacement",
        "method_random_state": method_random_state,
        "eligible_unit_ids_sha256": unit_id_digest(eligible.obs_names),
        "selected_unit_ids_sha256": unit_id_digest(selected.obs_names),
    }


def require_same_cohort(actual: dict[str, Any], expected: dict[str, Any]) -> None:
    keys = (
        "stage", "component", "side", "scope", "eligible_n_units", "selected_n_units",
        "selection_method", "selection_random_state", "selection_seed_applied",
        "method_random_state", "eligible_unit_ids_sha256", "selected_unit_ids_sha256",
    )
    changed = [key for key in keys if actual.get(key) != expected.get(key)]
    if changed:
        component, side, scope = actual["component"], actual["side"], actual["scope"]
        raise RuntimeError(f"formal {component} cohort drift for {side}/{scope}: {changed}")


def reproduce_cohorts(sample: Any, parameters: dict[str, Any], manifest: dict[str, Any]):
    if manifest.get("hash_contract") != UNIT_ID_HASH_CONTRACT:
        raise RuntimeError("formal cohort manifest uses an unexpected unit-ID hash contract")
    cohorts, records = {}, []
    for scope in SCOPES:
        for side in SIDES:
            eligible = _scope(sample, side, scope)
            seed = int(parameters["random_state"] + (side == "svc"))
            selected = deterministic_subset(eligible, parameters["sample_n_units"], seed)
            actual = actual_cohort_record(
                side=side, scope=scope, eligible=eligible, selected=selected, seed=seed,
            )
            require_same_cohort(actual, _baseline_cohort_record(manifest, side, scope))
            cohorts[side, scope] = selected
            records.append(actual)
    for scope in SCOPES[1:]:
        raw = next(row for row in records if row["side"] == "raw" and row["scope"] == scope)
        svc = next(row for row in records if row["side"] == "svc" and row["scope"] == scope)
        if raw["selected_unit_ids_sha256"] != svc["selected_unit_ids_sha256"]:
            raise RuntimeError(f"expected complete native parent cohorts to share ordered IDs for {scope}")
    return cohorts, records


def reproduce_pathway_cohorts(
    sample: Any, parameters: dict[str, Any], manifest: dict[str, Any],
):
    if manifest.get("hash_contract") != UNIT_ID_HASH_CONTRACT:
        raise RuntimeError("formal cohort manifest uses an unexpected unit-ID hash contract")
    scope = "|".join(SCOPES)
    cohorts, records = {}, []
    for side in SIDES:
        full = getattr(sample, side)
        labels = sample.labels(side) if sample.broad_key in full.obs else None
        ids = full.obs_names if "All" in SCOPES else labels.index[labels.isin(SCOPES)]
        eligible = full[ids]
        seed = int(parameters["random_state"] + (side == "svc"))
        selected = deterministic_subset(eligible, parameters["sample_n_units"], seed)
        actual = actual_cohort_record(
            side=side, scope=scope, eligible=eligible, selected=selected, seed=seed,
            component="pathway", method_random_state=int(parameters["random_state"]),
        )
        require_same_cohort(
            actual, _baseline_cohort_record(
                manifest, side, scope, component="pathway",
            ),
        )
        cohorts[side] = selected
        records.append(actual)
    return cohorts, records


def common_gene_axis(sample: Any) -> list[str]:
    raw_genes = {str(value) for value in sample.raw.var_names}
    svc_genes = {str(value) for value in sample.svc.var_names}
    genes = sorted(raw_genes & svc_genes)
    if len(genes) != EXPECTED_COMMON_GENES:
        raise RuntimeError(
            f"common gene axis drift: expected {EXPECTED_COMMON_GENES}, observed {len(genes)}"
        )
    return genes


def _validate_moran_table(table: pd.DataFrame, genes: list[str], side: str, scope: str) -> pd.DataFrame:
    required = {"gene_id", "moran", "status", "reason"}
    if not required.issubset(table.columns):
        raise ValueError(f"Moran output lacks columns: {sorted(required - set(table.columns))}")
    if table["gene_id"].duplicated().any():
        raise ValueError("Moran output has duplicate gene IDs")
    indexed = table.set_index(table["gene_id"].astype(str), drop=False)
    if set(indexed.index) != set(genes):
        raise RuntimeError(f"Moran output gene axis drift for {side}/{scope}")
    indexed = indexed.loc[genes].copy()
    values = pd.to_numeric(indexed["moran"], errors="coerce")
    computed = indexed["status"].eq("computed")
    if (~computed & values.notna()).any():
        raise RuntimeError("uncomputed Moran rows must retain NaN values")
    if (computed & ~np.isfinite(values)).any():
        raise RuntimeError("computed Moran rows must be finite")
    indexed["side"] = side
    indexed["scope"] = scope
    return indexed.reset_index(drop=True)


def compute_common_axis_tables(
    *, sample: Any, cohorts: dict[tuple[str, str], Any], records: list[dict[str, Any]],
    parameters: dict[str, Any], genes: list[str],
    moran_function: Callable[..., pd.DataFrame] = compute_moran,
) -> dict[tuple[str, str], pd.DataFrame]:
    record_index = {(row["side"], row["scope"]): row for row in records}
    axis_sha = ordered_digest(genes)
    tables = {}
    for scope in SCOPES:
        for side in SIDES:
            cohort = cohorts[side, scope]
            # This slice is deliberately the last expression operation before
            # the existing method performs normalize_total/log1p on its copy.
            work = cohort[:, genes].copy()
            if list(map(str, work.var_names)) != genes:
                raise RuntimeError(f"failed to impose common gene order for {side}/{scope}")
            table = moran_function(
                work,
                spatial_key=sample.spatial_key,
                n_neighbors=parameters["moran_n_neighbors"],
                transformation_state=sample.expression(side)["scale"],
            ).copy()
            table = _validate_moran_table(table, genes, side, scope)
            record = record_index[side, scope]
            table["n_eligible_units"] = record["eligible_n_units"]
            table["n_selected_units"] = record["selected_n_units"]
            table["selected_unit_ids_sha256"] = record["selected_unit_ids_sha256"]
            table["selection_random_state"] = record["selection_random_state"]
            table["selection_seed_applied"] = record["selection_seed_applied"]
            table["gene_axis_n"] = len(genes)
            table["gene_axis_sha256"] = axis_sha
            tables[side, scope] = table
    return tables


def detection_table(adata: Any, *, side: str, gene_axis_sha256: str) -> pd.DataFrame:
    detected = (
        np.asarray((adata.X != 0).sum(axis=1)).ravel()
        if sparse.issparse(adata.X) else np.count_nonzero(np.asarray(adata.X), axis=1)
    )
    return pd.DataFrame({
        "unit_id": adata.obs_names.astype(str),
        "side": side,
        "n_detected_genes": detected.astype(int),
        "detected_fraction": detected.astype(float) / adata.n_vars,
        "gene_axis_n": int(adata.n_vars),
        "gene_axis_sha256": gene_axis_sha256,
    })


def detection_summary(table: pd.DataFrame) -> dict[str, Any]:
    values = table["n_detected_genes"].to_numpy(dtype=float)
    quantiles = np.quantile(values, [0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1])
    return {
        "n_units": int(len(values)),
        "min": int(values.min()),
        "max": int(values.max()),
        "mean": float(values.mean()),
        "quantiles": {
            key: float(value) for key, value in zip(
                ("q0", "q01", "q05", "q25", "q50", "q75", "q95", "q99", "q100"),
                quantiles, strict=True,
            )
        },
        "value_counts": {
            str(int(key)): int(value)
            for key, value in table["n_detected_genes"].value_counts().sort_index().items()
        },
    }


def compute_common_axis_pathway(
    *, sample: Any, cohorts: dict[str, Any], records: list[dict[str, Any]],
    parameters: dict[str, Any], genes: list[str], emt_genes: list[str],
    pathway_function: Callable[..., dict[str, Any]] = compute_pathway_scores,
) -> dict[str, dict[str, Any]]:
    axis_sha = ordered_digest(genes)
    record_index = {row["side"]: row for row in records}
    outputs = {}
    for side in SIDES:
        cohort = cohorts[side]
        work = cohort[:, genes].copy()
        if list(map(str, work.var_names)) != genes:
            raise RuntimeError(f"failed to impose common gene order for pathway/{side}")
        detections = detection_table(work, side=side, gene_axis_sha256=axis_sha)
        requested = list(dict.fromkeys(map(str, emt_genes)))
        present = [gene for gene in requested if gene in set(genes)]
        expected_coverage = {
            "requested_genes": requested,
            "n_requested": len(requested),
            "present_genes": present,
            "n_present": len(present),
            "n_missing": len(requested) - len(present),
            "fraction_present": len(present) / len(requested),
        }
        detected = detections["n_detected_genes"]
        effective_fraction = float(
            detected.quantile(parameters["pathway_auc_threshold"]) / len(genes)
        )
        rank_cutoff = int(round(effective_fraction * len(genes)))
        base_metadata = {
            "side": side,
            "program": EMT_NAME,
            "status": "computed",
            "cohort": record_index[side],
            "gene_axis_n": len(genes),
            "gene_axis_sha256": axis_sha,
            "coverage": expected_coverage,
            "detection_distribution": detection_summary(detections),
            "detection_derived_cutoff": {
                "auc_threshold": parameters["pathway_auc_threshold"],
                "effective_rank_fraction": effective_fraction,
                "rank_cutoff": rank_cutoff,
                "valid": 1 < rank_cutoff < len(genes),
                "valid_range": [2, len(genes) - 1],
            },
        }
        try:
            result = pathway_function(
                work,
                genes=emt_genes,
                score_name=EMT_NAME,
                auc_threshold=parameters["pathway_auc_threshold"],
                seed=parameters["random_state"],
            )
        except (ImportError, PrerequisiteUnavailable) as exc:
            outputs[side] = {
                "status": "unavailable",
                "reason": str(exc),
                "scores": None,
                "detections": detections,
                "metadata": {**base_metadata, "status": "unavailable", "reason": str(exc)},
            }
            continue
        scores = result.get("scores")
        if not isinstance(scores, pd.Series):
            raise RuntimeError("AUCell result has no score Series")
        if scores.index.astype(str).tolist() != work.obs_names.astype(str).tolist():
            raise RuntimeError(f"AUCell score cohort/order drift for {side}")
        if not np.isfinite(scores.to_numpy(dtype=float)).all():
            raise RuntimeError(f"AUCell returned nonfinite scores for {side}")
        coverage = result.get("coverage")
        metadata = result.get("metadata")
        if not isinstance(coverage, dict) or not isinstance(metadata, dict):
            raise RuntimeError("AUCell result lacks coverage/provider metadata")
        if coverage != expected_coverage:
            raise RuntimeError(f"AUCell coverage drift for {side}")
        if metadata.get("rank_cutoff") != rank_cutoff:
            raise RuntimeError(f"AUCell detection-derived cutoff drift for {side}")
        if coverage.get("requested_genes") != requested:
            raise RuntimeError(f"AUCell requested-gene identity drift for {side}")
        score_table = detections.copy()
        score_table["score"] = scores.to_numpy(dtype=float)
        score_table["program"] = EMT_NAME
        outputs[side] = {
            "status": "computed",
            "reason": None,
            "scores": score_table,
            "detections": detections,
            "metadata": {
                **base_metadata,
                "coverage": coverage,
                "provider_metadata": metadata,
            },
        }
    return outputs


def joint_table(raw: pd.DataFrame, svc: pd.DataFrame, genes: list[str], scope: str) -> pd.DataFrame:
    left = raw.set_index("gene_id").loc[genes]
    right = svc.set_index("gene_id").loc[genes]
    table = pd.DataFrame({
        "gene_id": genes,
        "raw_moran": left["moran"].to_numpy(),
        "svc_moran": right["moran"].to_numpy(),
        "raw_status": left["status"].to_numpy(),
        "svc_status": right["status"].to_numpy(),
        "raw_reason": left["reason"].to_numpy(),
        "svc_reason": right["reason"].to_numpy(),
    })
    table["delta_moran"] = table["svc_moran"] - table["raw_moran"]
    table["comparison_available"] = (
        table["raw_status"].eq("computed") & table["svc_status"].eq("computed")
        & np.isfinite(table["raw_moran"]) & np.isfinite(table["svc_moran"])
    )
    table["scope"] = scope
    return table


def read_baseline_table(path: Path, genes: list[str], scope: str) -> pd.DataFrame:
    table = pd.read_csv(path)
    required = {"gene_id", "raw_moran", "svc_moran", "delta_moran"}
    if not required.issubset(table.columns):
        raise ValueError(f"baseline Moran table lacks columns: {sorted(required - set(table.columns))}")
    if table["gene_id"].astype(str).duplicated().any():
        raise ValueError(f"baseline Moran table has duplicate genes for {scope}")
    table["gene_id"] = table["gene_id"].astype(str)
    if set(table["gene_id"]) != set(genes):
        raise RuntimeError(f"baseline Moran shared-gene identity drift for {scope}")
    return table.set_index("gene_id").loc[genes].reset_index()


def delta_change_table(
    baseline: pd.DataFrame, common_axis: pd.DataFrame, genes: list[str], scope: str,
) -> pd.DataFrame:
    old = baseline.set_index("gene_id").loc[genes]
    new = common_axis.set_index("gene_id").loc[genes]
    table = pd.DataFrame({
        "gene_id": genes,
        "baseline_raw_moran_native_axis": old["raw_moran"].to_numpy(),
        "baseline_svc_moran_native_axis": old["svc_moran"].to_numpy(),
        "baseline_delta_moran_native_axis": old["delta_moran"].to_numpy(),
        "common_axis_raw_moran": new["raw_moran"].to_numpy(),
        "common_axis_svc_moran": new["svc_moran"].to_numpy(),
        "common_axis_delta_moran": new["delta_moran"].to_numpy(),
        "baseline_comparison_available": old.get(
            "comparison_available", pd.Series(False, index=old.index)
        ).astype(bool).to_numpy(),
        "common_axis_comparison_available": new["comparison_available"].to_numpy(),
    })
    table["raw_moran_axis_change"] = (
        table["common_axis_raw_moran"] - table["baseline_raw_moran_native_axis"]
    )
    table["svc_moran_axis_change"] = (
        table["common_axis_svc_moran"] - table["baseline_svc_moran_native_axis"]
    )
    table["delta_moran_axis_change"] = (
        table["common_axis_delta_moran"] - table["baseline_delta_moran_native_axis"]
    )
    table["scope"] = scope
    return table


def read_baseline_pathway_scores(
    path: Path, *, expected_ids: list[str], expected_hash: str, side: str,
) -> pd.DataFrame:
    table = pd.read_csv(path)
    required = {"unit_id", "score", "side", "program", "scope"}
    if not required.issubset(table.columns):
        raise ValueError(
            f"baseline pathway table lacks columns for {side}: {sorted(required - set(table.columns))}"
        )
    if table["unit_id"].astype(str).duplicated().any():
        raise ValueError(f"baseline pathway table has duplicate unit IDs for {side}")
    ids = table["unit_id"].astype(str).tolist()
    if ids != expected_ids or unit_id_digest(ids) != expected_hash:
        raise RuntimeError(f"baseline pathway score cohort/order drift for {side}")
    if set(table["side"].astype(str)) != {side}:
        raise RuntimeError(f"baseline pathway side drift for {side}")
    if set(table["program"].astype(str)) != {EMT_NAME}:
        raise RuntimeError(f"baseline pathway program drift for {side}")
    if set(table["scope"].astype(str)) != {"|".join(SCOPES)}:
        raise RuntimeError(f"baseline pathway scope drift for {side}")
    if not np.isfinite(pd.to_numeric(table["score"], errors="coerce")).all():
        raise RuntimeError(f"baseline pathway scores are nonfinite for {side}")
    return pd.DataFrame({
        "unit_id": ids,
        "baseline_score_native_axis": pd.to_numeric(table["score"]).to_numpy(dtype=float),
    })


def validate_baseline_pathway_metadata(
    payload: dict[str, Any], *, sample: Any, side: str, record: dict[str, Any],
    emt_genes: list[str], current_provider_metadata: dict[str, Any] | None,
) -> None:
    required = {
        "coverage", "metadata", "scoring_scopes", "n_eligible_units",
        "n_scored_units", "gene_axis", "expression",
    }
    if not required.issubset(payload):
        raise ValueError(
            f"baseline pathway metadata lacks fields for {side}: {sorted(required - set(payload))}"
        )
    expected_requested = list(dict.fromkeys(map(str, emt_genes)))
    coverage = payload["coverage"]
    provider = payload["metadata"]
    if not isinstance(coverage, dict) or not isinstance(provider, dict):
        raise ValueError(f"baseline pathway metadata has invalid coverage/provider payload for {side}")
    if coverage.get("requested_genes") != expected_requested:
        raise RuntimeError(f"baseline EMT requested-gene identity drift for {side}")
    native_genes = list(map(str, getattr(sample, side).var_names))
    expected_present = [gene for gene in expected_requested if gene in set(native_genes)]
    if coverage.get("present_genes") != expected_present:
        raise RuntimeError(f"baseline EMT present-gene identity drift for {side}")
    if payload["scoring_scopes"] != list(SCOPES):
        raise RuntimeError(f"baseline pathway scoring scopes drift for {side}")
    if payload["n_eligible_units"] != record["eligible_n_units"]:
        raise RuntimeError(f"baseline pathway eligible cohort size drift for {side}")
    if payload["n_scored_units"] != record["selected_n_units"]:
        raise RuntimeError(f"baseline pathway selected cohort size drift for {side}")
    if list(map(str, payload["gene_axis"])) != native_genes:
        raise RuntimeError(f"baseline pathway native gene axis drift for {side}")
    if payload["expression"] != sample.expression(side):
        raise RuntimeError(f"baseline pathway expression declaration drift for {side}")
    expected_provider = {
        "score_name": EMT_NAME,
        "score_key": f"{EMT_NAME}_aucell",
        "auc_threshold": EXPECTED_PATHWAY_AUC_THRESHOLD,
        "seed": EXPECTED_RANDOM_STATE,
    }
    changed = [key for key, value in expected_provider.items() if provider.get(key) != value]
    if changed:
        raise RuntimeError(f"baseline pathway method metadata drift for {side}: {changed}")
    rank_cutoff = provider.get("rank_cutoff")
    if not isinstance(rank_cutoff, int) or not 1 < rank_cutoff < len(native_genes):
        raise RuntimeError(f"baseline pathway rank cutoff is invalid for {side}")
    for key in ("provider", "provider_version", "scorer", "provider_module"):
        if not isinstance(provider.get(key), str) or not provider[key]:
            raise ValueError(f"baseline pathway provider metadata lacks {key} for {side}")
    memory = provider.get("memory_adapter")
    if not isinstance(memory, dict) or memory.get("call_count") != 1 or memory.get("restored_original") is not True:
        raise RuntimeError(f"baseline pathway adapter audit is invalid for {side}")
    if current_provider_metadata is None:
        return
    invariant_provider_keys = (
        "provider", "provider_version", "scorer", "provider_module", "score_name",
        "score_key", "auc_threshold", "seed",
    )
    provider_changes = [
        key for key in invariant_provider_keys
        if provider.get(key) != current_provider_metadata.get(key)
    ]
    current_memory = current_provider_metadata.get("memory_adapter")
    if not isinstance(current_memory, dict):
        raise ValueError(f"common-axis pathway metadata lacks adapter audit for {side}")
    invariant_adapter_keys = (
        "name", "active", "seed", "rank_dtype", "permutation_hash_contract",
        "omicverse_version", "ctxcore_version", "public_scorer",
        "source_hash_contract", "original_source_hashes", "adapter_source_hash",
        "call_count", "full_cohort_single_provider_call", "restored_original",
    )
    adapter_changes = [
        key for key in invariant_adapter_keys if memory.get(key) != current_memory.get(key)
    ]
    if provider_changes or adapter_changes:
        raise RuntimeError(
            f"pathway provider/adapter provenance drift for {side}: "
            f"provider={provider_changes}, adapter={adapter_changes}"
        )


def pathway_score_change_table(
    baseline: pd.DataFrame, common_axis: pd.DataFrame, *, side: str,
) -> pd.DataFrame:
    if baseline["unit_id"].tolist() != common_axis["unit_id"].tolist():
        raise RuntimeError(f"pathway score order drift while comparing {side}")
    table = common_axis.copy()
    table.insert(
        table.columns.get_loc("score"),
        "baseline_score_native_axis",
        baseline["baseline_score_native_axis"].to_numpy(dtype=float),
    )
    table = table.rename(columns={"score": "common_axis_score"})
    table["score_axis_change"] = (
        table["common_axis_score"] - table["baseline_score_native_axis"]
    )
    return table


def _status_summary(table: pd.DataFrame) -> dict[str, Any]:
    counts = {str(key): int(value) for key, value in table["status"].value_counts(dropna=False).items()}
    reasons = {
        str(key): int(value)
        for key, value in table.loc[table["status"] != "computed", "reason"].value_counts(dropna=False).items()
    }
    return {"status_counts": counts, "exclusion_reason_counts": reasons}


def _write_outputs(
    staging: Path, *, run: dict[str, Any], provenance: dict[str, Any], sample: Any,
    genes: list[str], records: list[dict[str, Any]],
    side_tables: dict[tuple[str, str], pd.DataFrame],
    pathway_outputs: dict[str, dict[str, Any]], emt_genes: list[str],
) -> dict[str, Any]:
    outputs = []

    def save_csv(relative: str, table: pd.DataFrame) -> None:
        path = staging / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(path, index=False)
        outputs.append(relative)

    def save_json(relative: str, value: Any) -> None:
        write_json(staging / relative, value)
        outputs.append(relative)

    axis_sha = ordered_digest(genes)
    save_csv("tables/gene_axis.csv", pd.DataFrame({
        "axis_index": np.arange(len(genes), dtype=int), "gene_id": genes,
    }))
    save_json("tables/gene_axis.json", {
        "schema_version": 1,
        "ordering": "ascending Unicode code point order",
        "hash_contract": GENE_AXIS_HASH_CONTRACT,
        "n_genes": len(genes),
        "ordered_gene_ids_sha256": axis_sha,
        "genes": genes,
        "raw_native_n_genes": int(sample.raw.n_vars),
        "svc_native_n_genes": int(sample.svc.n_vars),
    })
    save_json("tables/cohort_manifest.json", {
        "schema_version": 1,
        "sample_id": SAMPLE_ID,
        "hash_contract": UNIT_ID_HASH_CONTRACT,
        "baseline_manifest_sha256": provenance["artifact_identities"]["cohort_manifest"]["sha256"],
        "records": records,
        "parent_cross_side_identity": {
            scope: {
                "ordered_selected_ids_identical": True,
                "pairing_used_by_moran": False,
            }
            for scope in SCOPES[1:]
        },
    })

    summaries = {}
    unavailable = []
    exclusions = []
    for scope in SCOPES:
        for side in SIDES:
            table = side_tables[side, scope]
            save_csv(f"tables/moran_common_axis_{side}_{scope}.csv", table)
            summaries[f"{side}:{scope}"] = _status_summary(table)
            excluded = table.loc[table["status"] != "computed",
                                 ["gene_id", "side", "scope", "status", "reason"]]
            exclusions.append(excluded)
        joint = joint_table(side_tables["raw", scope], side_tables["svc", scope], genes, scope)
        save_csv(f"tables/moran_common_axis_joint_{scope}.csv", joint)
        baseline_path = run["analysis_dir"] / f"tables/moran_shared_genes_{scope}.csv"
        baseline = read_baseline_table(baseline_path, genes, scope)
        changes = delta_change_table(baseline, joint, genes, scope)
        save_csv(f"tables/moran_common_axis_delta_change_{scope}.csv", changes)
        summaries[f"joint:{scope}"] = {
            "n_genes": len(joint),
            "n_jointly_computable": int(joint["comparison_available"].sum()),
            "n_negative_delta": int((joint.loc[joint["comparison_available"], "delta_moran"] < 0).sum()),
            "n_positive_delta": int((joint.loc[joint["comparison_available"], "delta_moran"] > 0).sum()),
            "n_zero_delta": int((joint.loc[joint["comparison_available"], "delta_moran"] == 0).sum()),
        }
    exclusion_table = pd.concat(exclusions, ignore_index=True) if exclusions else pd.DataFrame(
        columns=["gene_id", "side", "scope", "status", "reason"]
    )
    save_csv("tables/moran_common_axis_exclusions.csv", exclusion_table)

    pathway_record_index = {
        row["side"]: row for row in records if row["component"] == "pathway"
    }
    for side in SIDES:
        pathway = pathway_outputs[side]
        prefix = f"tables/pathway_common_axis_{side}_{EMT_NAME}"
        save_csv(f"{prefix}_detection.csv", pathway["detections"])
        baseline_metadata_path = (
            run["analysis_dir"] / f"tables/pathway_{side}_{EMT_NAME}_metadata.json"
        )
        baseline_metadata = read_json(baseline_metadata_path)
        record = pathway_record_index[side]
        current_provider_metadata = (
            pathway["metadata"].get("provider_metadata")
            if pathway["status"] == "computed" else None
        )
        validate_baseline_pathway_metadata(
            baseline_metadata, sample=sample, side=side, record=record,
            emt_genes=emt_genes, current_provider_metadata=current_provider_metadata,
        )
        metadata_payload = {
            **pathway["metadata"],
            "baseline_native_axis_metadata": baseline_metadata,
        }
        save_json(f"{prefix}_metadata.json", metadata_payload)
        if pathway["status"] == "unavailable":
            unavailable.append({
                "component": f"pathway:{side}:{EMT_NAME}",
                "reason": pathway["reason"],
            })
            summaries[f"pathway:{side}"] = {
                "status": "unavailable", "reason": pathway["reason"],
            }
            continue
        scores = pathway["scores"]
        save_csv(f"{prefix}_scores.csv", scores)
        baseline_scores_path = run["analysis_dir"] / f"tables/pathway_{side}_{EMT_NAME}.csv"
        baseline_scores = read_baseline_pathway_scores(
            baseline_scores_path,
            expected_ids=scores["unit_id"].tolist(),
            expected_hash=record["selected_unit_ids_sha256"],
            side=side,
        )
        changes = pathway_score_change_table(baseline_scores, scores, side=side)
        save_csv(f"{prefix}_score_change.csv", changes)
        summaries[f"pathway:{side}"] = {
            "status": "computed",
            "n_scored_units": int(len(scores)),
            "score_min": float(scores["score"].min()),
            "score_max": float(scores["score"].max()),
            "score_mean": float(scores["score"].mean()),
            "coverage": pathway["metadata"]["coverage"],
            "rank_cutoff": pathway["metadata"]["provider_metadata"].get("rank_cutoff"),
        }

    result = {
        "schema_version": 1,
        "task_id": TASK_ID,
        "iteration": ITERATION,
        "diagnostic_id": DIAGNOSTIC_ID,
        "status": "diagnostic_complete",
        "acceptance_level": "technical_diagnostic_complete_scientific_review_pending",
        "scientific_acceptance": False,
        "question": "How do Moran and EMT AUCell results change when each existing native cohort is analyzed after restricting both sides to the same genes?",
        "interpretation_boundary": (
            "This isolates the common-gene background while preserving each native unit cohort and all method parameters. "
            "A Moran or EMT-score increase or decrease is a diagnostic change, not proof of biological improvement or deterioration."
        ),
        "fixed_parameters": {
            "scopes": list(SCOPES),
            "sample_n_units": EXPECTED_SAMPLE_N_UNITS,
            "raw_selection_seed": EXPECTED_RANDOM_STATE,
            "svc_selection_seed": EXPECTED_RANDOM_STATE + 1,
            "moran_n_neighbors": EXPECTED_MORAN_NEIGHBORS,
            "normalization": "normalize_total(target_sum=1e4) then log1p inside methods.spatial.compute_moran",
            "pathway_program": EMT_NAME,
            "pathway_auc_threshold": EXPECTED_PATHWAY_AUC_THRESHOLD,
            "pathway_method_seed": EXPECTED_RANDOM_STATE,
            "unit_pairing": False,
            "gene_axis_n": len(genes),
            "gene_axis_sha256": axis_sha,
        },
        "summaries": summaries,
        "unavailable": unavailable,
        "outputs": sorted(outputs),
    }
    save_json("result.json", result)
    return result


def build_artifact_manifest(staging: Path, final_dir: Path) -> dict[str, Any]:
    rows = []
    for path in sorted(staging.rglob("*")):
        if not path.is_file() or path.name in {"run_receipt.json", "artifact_hashes.json"}:
            continue
        identity = file_identity(path)
        relative = path.relative_to(staging)
        identity["path"] = str(final_dir / relative)
        identity["relative_to_diagnostic"] = str(relative)
        rows.append(identity)
    return {
        "schema_version": 1,
        "generated_at": now(),
        "artifacts": rows,
        "missing_artifacts": [],
    }


def preflight(
    repo: Path, config_path: Path, *, baseline_run: Path | None = None,
    output: Path | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    run = resolve_run(repo, config_path, baseline_run=baseline_run, output=output)
    require_output_absent(run["diagnostic_dir"])
    repository = git_identity(repo)
    require_clean_worktree(repository)
    provenance = verify_formal_provenance(run)
    manifest = read_json(run["cohort_manifest_path"])
    if manifest.get("schema_version") != 1 or manifest.get("sample_id") != SAMPLE_ID:
        raise ValueError("formal cohort manifest schema/sample identity is invalid")
    expected = {(side, scope) for side in SIDES for scope in SCOPES}
    found = {
        (row.get("side"), row.get("scope"))
        for row in manifest.get("records", [])
        if row.get("stage") == "molecular" and row.get("component") == "moran"
    }
    if found != expected:
        raise RuntimeError("formal cohort manifest lacks the exact Raw/SVC Moran scope matrix")
    pathway_scope = "|".join(SCOPES)
    pathway_found = {
        (row.get("side"), row.get("scope"))
        for row in manifest.get("records", [])
        if row.get("stage") == "molecular" and row.get("component") == "pathway"
    }
    if pathway_found != {(side, pathway_scope) for side in SIDES}:
        raise RuntimeError("formal cohort manifest lacks the exact Raw/SVC pathway union cohorts")
    provenance["repository"] = repository
    provenance["formal_cohort_manifest"] = manifest
    return run, provenance


def execute(
    repo: Path, config_path: Path, *, baseline_run: Path | None = None,
    output: Path | None = None,
) -> dict[str, Any]:
    started_at = now()
    run, provenance = preflight(
        repo, config_path, baseline_run=baseline_run, output=output,
    )
    diagnostics_root = run["diagnostic_dir"].parent
    diagnostics_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{DIAGNOSTIC_ID}-", dir=diagnostics_root))
    receipt = {
        "schema_version": 1,
        "task_id": TASK_ID,
        "iteration": ITERATION,
        "diagnostic_id": DIAGNOSTIC_ID,
        "owner": "REVISE_Analysis_Agent",
        "upstream": "formal P2CRC_Xenium/random Impact iter-001",
        "downstream": "independent scientific interpretation",
        "status": "running",
        "acceptance_level": "scientific_review_pending",
        "scientific_acceptance": False,
        "started_at": started_at,
        "observed_at": now(),
        "execution_location": socket.gethostname(),
        "command": [sys.executable, *sys.argv],
        "repository": provenance["repository"],
        "environment": environment_identity(),
        "input_identity": provenance["input_identities"],
        "baseline_artifact_identity": provenance["artifact_identities"],
        "code_identity": provenance["code_identities"],
        "output_dir": str(run["diagnostic_dir"]),
    }
    write_json(staging / "run_receipt.json", receipt)
    try:
        sample = load_sample(run["sample_yaml"])
        for side in SIDES:
            reason = sample.expression_unavailable(side)
            if reason:
                raise RuntimeError(reason)
        cohorts, moran_records = reproduce_cohorts(
            sample, run["parameters"], provenance["formal_cohort_manifest"]
        )
        pathway_cohorts, pathway_records = reproduce_pathway_cohorts(
            sample, run["parameters"], provenance["formal_cohort_manifest"]
        )
        records = [*moran_records, *pathway_records]
        genes = common_gene_axis(sample)
        side_tables = compute_common_axis_tables(
            sample=sample, cohorts=cohorts, records=moran_records,
            parameters=run["parameters"], genes=genes,
        )
        gene_sets, gene_set_error = resolve_gene_sets(run["parameters"])
        if gene_set_error:
            raise RuntimeError(gene_set_error)
        if list(gene_sets) != [EMT_NAME]:
            raise RuntimeError(f"expected exactly the EMT gene set {EMT_NAME}")
        emt_genes = gene_sets[EMT_NAME]
        pathway_outputs = compute_common_axis_pathway(
            sample=sample, cohorts=pathway_cohorts, records=pathway_records,
            parameters=run["parameters"], genes=genes, emt_genes=emt_genes,
        )
        result = _write_outputs(
            staging, run=run, provenance=provenance, sample=sample,
            genes=genes, records=records, side_tables=side_tables,
            pathway_outputs=pathway_outputs, emt_genes=emt_genes,
        )
        artifacts = build_artifact_manifest(staging, run["diagnostic_dir"])
        write_json(staging / "artifact_hashes.json", artifacts)
        result_identity = file_identity(staging / "result.json")
        result_identity["path"] = str(run["diagnostic_dir"] / "result.json")
        artifact_manifest_identity = file_identity(staging / "artifact_hashes.json")
        artifact_manifest_identity["path"] = str(
            run["diagnostic_dir"] / "artifact_hashes.json"
        )
        receipt.update(
            status="diagnostic_complete",
            result_status=result["status"],
            result="result.json",
            result_identity=result_identity,
            artifact_manifest="artifact_hashes.json",
            artifact_manifest_identity=artifact_manifest_identity,
            finished_at=now(), observed_at=now(),
            peak_rss_kib=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        )
        write_json(staging / "run_receipt.json", receipt)
        reverify_provenance(run, provenance)
        if run["diagnostic_dir"].exists():
            raise FileExistsError(f"diagnostic output appeared during execution: {run['diagnostic_dir']}")
        os.rename(staging, run["diagnostic_dir"])
        return receipt
    except BaseException as exc:
        receipt.update(
            status="diagnostic_failed",
            finished_at=now(), observed_at=now(),
            peak_rss_kib=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
            error={
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
            failed_staging_dir=str(staging),
        )
        write_json(staging / "run_receipt.json", receipt)
        raise


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("preflight", "run"))
    parser.add_argument("--repo")
    parser.add_argument("--config")
    parser.add_argument("--baseline-run")
    parser.add_argument("--output")
    args = parser.parse_args(argv)
    repo = Path(args.repo).expanduser().resolve() if args.repo else Path(__file__).resolve().parents[1]
    config_path = Path(args.config) if args.config else repo / DEFAULT_CONFIG
    baseline_run = Path(args.baseline_run) if args.baseline_run else None
    output = Path(args.output) if args.output else None
    if args.command == "preflight":
        run, provenance = preflight(
            repo, config_path, baseline_run=baseline_run, output=output,
        )
        value = {
            "status": "ready",
            "task_id": TASK_ID,
            "iteration": ITERATION,
            "diagnostic_id": DIAGNOSTIC_ID,
            "formal_result": str(run["result_path"]),
            "formal_cohort_manifest": str(run["cohort_manifest_path"]),
            "diagnostic_dir": str(run["diagnostic_dir"]),
            "verified_inputs": provenance["input_identities"],
            "verified_baseline_artifacts": provenance["artifact_identities"],
            "code_identity": provenance["code_identities"],
        }
        print(json.dumps(value, ensure_ascii=False, indent=2))
        return 0
    receipt = execute(repo, config_path, baseline_run=baseline_run, output=output)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0 if receipt["status"] == "diagnostic_complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
