#!/usr/bin/env python3
"""Run one guarded P2 Impact delivery and preserve an auditable receipt."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import platform
import resource
import re
import socket
import subprocess
import sys
import threading
import traceback

from revise_analysis.batch import run_batch
from revise_analysis.io import read_yaml, resolve_path
from revise_analysis.runner import write_json


TASK_ID = "p2-xenium-impact-20260922"
DEFAULT_ITERATION = "iter-001"
EXPECTED_SCOPES = ["All", "Fibroblast", "Mono_Macro", "T"]
EXPECTED_ANALYSIS = "reconstruction_impact"
EXPECTED_CONFIG = Path("configs/p2_xenium_impact_20260922_iter_001.yaml")
BASELINE_CONFIG = Path("configs/p2_full_project.yaml")
PACKAGE_NAMES = (
    "revise-analysis-agent", "anndata", "numpy", "pandas", "scipy", "scanpy",
    "omicverse", "gseapy", "torch", "nbconvert",
)


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Path) -> dict:
    stat = path.stat()
    return {"path": str(path), "size_bytes": stat.st_size, "sha256": sha256(path)}


def git_identity(repo: Path) -> dict:
    def run(*args):
        return subprocess.run(("git", "-C", str(repo), *args), check=True,
                              text=True, capture_output=True).stdout.strip()
    dirty = run("status", "--porcelain").splitlines()
    return {
        "commit": run("rev-parse", "HEAD"),
        "branch": run("branch", "--show-current"),
        "dirty": bool(dirty),
        "dirty_paths": dirty,
    }


def environment_identity() -> dict:
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


def cgroup_identity() -> dict:
    values = {}
    for name in ("cpu.max", "memory.max", "memory.current"):
        path = Path("/sys/fs/cgroup") / name
        if path.exists():
            values[name.replace(".", "_")] = path.read_text(encoding="utf-8").strip()
    return values


def normalized_parameters(parameters: dict) -> dict:
    value = dict(parameters)
    for key, default in (("moran_n_neighbors", 6), ("raw_k_control", False),
                         ("svc_leiden_baseline", False)):
        value.setdefault(key, default)
    return value


def require_frozen_parameters(parameters: dict, repo: Path) -> None:
    baseline = read_yaml(repo / BASELINE_CONFIG)["analyses"][EXPECTED_ANALYSIS]
    expected = normalized_parameters(baseline)
    actual = normalized_parameters(parameters)
    if actual != expected:
        changed = sorted(key for key in set(actual) | set(expected)
                         if actual.get(key) != expected.get(key))
        raise ValueError(f"Impact parameters differ from the frozen P2 full baseline: {changed}")


def require_clean_worktree(identity: dict) -> None:
    if identity["dirty"]:
        raise RuntimeError("formal execution requires a clean committed worktree")


def execution_status(analysis_status: str | None, missing_artifacts: list[str]) -> str:
    if missing_artifacts:
        return "artifact_validation_failed"
    return ("execution_complete" if analysis_status in {"succeeded", "partial"}
            else "execution_failed")


def iteration_segment(value: str) -> str:
    if re.fullmatch(r"iter-[0-9]{3}", value) is None:
        raise ValueError("iteration must match iter-NNN")
    return value


def expected_run_paths(repo: Path, iteration: str, producer_iteration: str) -> tuple[Path, Path]:
    output = (repo / "output" / TASK_ID / iteration_segment(iteration)).resolve()
    sample = (repo.parent / "REVISE" / "results" / TASK_ID
              / iteration_segment(producer_iteration) / "P2CRC_Xenium" / "random"
              / "sample.yaml").resolve()
    return output, sample


def validate_project(config_path: Path, repo: Path, *, require_inputs: bool,
                     iteration=DEFAULT_ITERATION, producer_iteration=DEFAULT_ITERATION) -> dict:
    config_path = config_path.expanduser().resolve()
    config = read_yaml(config_path)
    analyses = config.get("analyses")
    if not isinstance(analyses, dict) or list(analyses) != [EXPECTED_ANALYSIS]:
        raise ValueError(f"analyses must contain only {EXPECTED_ANALYSIS!r}")
    parameters = analyses[EXPECTED_ANALYSIS]
    if parameters.get("scopes") != EXPECTED_SCOPES:
        raise ValueError(f"Impact scopes must be {EXPECTED_SCOPES}")
    require_frozen_parameters(parameters, repo)
    samples = config.get("samples")
    if not isinstance(samples, list) or len(samples) != 1:
        raise ValueError("project must declare exactly one sample")
    output_dir = resolve_path(config.get("output_dir"), config_path.parent)
    expected_output, expected_sample = expected_run_paths(repo, iteration, producer_iteration)
    if output_dir != expected_output:
        raise ValueError(f"output_dir must resolve to {expected_output}")
    sample_yaml = resolve_path(samples[0], config_path.parent)
    if sample_yaml != expected_sample:
        raise ValueError(f"sample must resolve to {expected_sample}")
    geneset_path = resolve_path(parameters.get("geneset_path"), config_path.parent)
    required = [config_path, geneset_path]
    if require_inputs:
        required.append(sample_yaml)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"required files are absent: {missing}")
    raw = svc = None
    sample_id = "P2CRC_Xenium"
    if sample_yaml.is_file():
        sample = read_yaml(sample_yaml)
        sample_id = sample.get("sample_id")
        if sample_id != "P2CRC_Xenium":
            raise ValueError("formal sample_id must be P2CRC_Xenium")
        files = sample.get("files", {})
        raw = resolve_path(files.get("raw"), sample_yaml.parent)
        svc = resolve_path(files.get("svc"), sample_yaml.parent)
        if require_inputs:
            absent = [str(path) for path in (raw, svc) if not path.is_file()]
            if absent:
                raise FileNotFoundError(f"sample objects are absent: {absent}")
    return {
        "config_path": config_path, "config": config, "parameters": parameters,
        "output_dir": output_dir, "sample_yaml": sample_yaml, "sample_id": sample_id,
        "raw": raw, "svc": svc, "geneset_path": geneset_path,
    }


def input_identities(preflight: dict) -> dict:
    files = {
        "project_config": preflight["config_path"],
        "sample_config": preflight["sample_yaml"],
        "raw": preflight["raw"], "svc": preflight["svc"],
        "geneset": preflight["geneset_path"],
    }
    return {name: file_identity(path) for name, path in files.items()}


def registered_artifacts(preflight: dict) -> dict:
    output_dir = preflight["output_dir"]
    sample_dir = output_dir / preflight["sample_id"]
    analysis_dir = sample_dir / EXPECTED_ANALYSIS
    result_path = analysis_dir / "result.json"
    result = json.loads(result_path.read_text(encoding="utf-8"))
    candidates = {result_path, sample_dir / "index.json", output_dir / "batch_result.json"}
    for relative in result.get("outputs", {}).values():
        path = (analysis_dir / relative).resolve()
        path.relative_to(analysis_dir.resolve())
        candidates.add(path)
    rows, missing = [], []
    for path in sorted(candidates):
        if path.is_file():
            identity = file_identity(path)
            identity["relative_to_output"] = str(path.relative_to(output_dir))
            rows.append(identity)
        else:
            missing.append(str(path.relative_to(output_dir)))
    return {"schema_version": 1, "generated_at": now(), "artifacts": rows,
            "missing_registered_artifacts": missing}


def stage_statuses(preflight: dict) -> dict:
    result_path = (preflight["output_dir"] / preflight["sample_id"]
                   / EXPECTED_ANALYSIS / "result.json")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    return {
        "analysis_status": result.get("status"),
        "configured_scopes": list(preflight["parameters"]["scopes"]),
        "stages": result.get("stages", []),
        "unavailable": result.get("unavailable", []),
        "stage_errors": result.get("stage_errors", []),
        "cohort_manifest": (f"{preflight['sample_id']}/{EXPECTED_ANALYSIS}/"
                            "tables/cohort_manifest.json"),
    }


def base_receipt(args, preflight: dict, repository: dict) -> dict:
    return {
        "schema_version": 1,
        "task_id": TASK_ID,
        "iteration": args.iteration,
        "producer_iteration": args.producer_iteration,
        "owner": "REVISE_Analysis_Agent",
        "upstream": "REVISE producer P2CRC_Xenium/random",
        "downstream": "independent scientific review gate",
        "status": "preparing",
        "acceptance_level": "scientific_review_pending",
        "scientific_acceptance": False,
        "started_at": now(),
        "observed_at": now(),
        "execution_location": args.execution_location,
        "qz_notebook": args.qz_notebook,
        "command": [sys.executable, *sys.argv],
        "execution_model": "single_process_no_detached_child",
        "resubmission_guard": "receipt or Impact result in this iteration blocks execution",
        "repository": repository,
        "environment": environment_identity(),
        "resource_limits": cgroup_identity(),
        "project_config": str(preflight["config_path"]),
        "output_dir": str(preflight["output_dir"]),
        "evidence_entry": str(preflight["output_dir"] / preflight["sample_id"]
                              / EXPECTED_ANALYSIS / "result.json"),
    }


def execute(args) -> dict:
    repo = Path(__file__).resolve().parents[1]
    config_path = Path(args.config) if args.config else repo / EXPECTED_CONFIG
    preflight = validate_project(config_path, repo, require_inputs=True,
                                 iteration=args.iteration,
                                 producer_iteration=args.producer_iteration)
    repository = git_identity(repo)
    require_clean_worktree(repository)
    run_dir = preflight["output_dir"] / "_run"
    receipt_path = run_dir / "run_receipt.json"
    artifacts_path = run_dir / "artifact_hashes.json"
    analysis_dir = preflight["output_dir"] / preflight["sample_id"] / EXPECTED_ANALYSIS
    if receipt_path.exists() or analysis_dir.exists():
        raise FileExistsError("iteration already contains a receipt or Impact result; use a new iteration")
    run_dir.mkdir(parents=True, exist_ok=True)
    receipt = base_receipt(args, preflight, repository)
    write_json(receipt_path, receipt)
    stop_telemetry = threading.Event()
    receipt_lock = threading.Lock()

    def persist():
        with receipt_lock:
            write_json(receipt_path, receipt)

    def telemetry():
        while not stop_telemetry.wait(30):
            with receipt_lock:
                receipt["observed_at"] = now()
                receipt["peak_rss_kib"] = int(
                    resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                write_json(receipt_path, receipt)

    telemetry_thread = threading.Thread(target=telemetry, name="run-receipt-telemetry", daemon=True)
    telemetry_thread.start()
    try:
        receipt["status"] = "hashing_inputs"
        receipt["observed_at"] = now()
        persist()
        receipt["input_identity"] = input_identities(preflight)
        receipt["status"] = "running"
        receipt["observed_at"] = now()
        persist()
        batch = run_batch(preflight["config_path"])
        manifest = registered_artifacts(preflight)
        write_json(artifacts_path, manifest)
        statuses = stage_statuses(preflight)
        receipt.update(
            status=execution_status(statuses["analysis_status"],
                                    manifest["missing_registered_artifacts"]),
            batch_status=batch.get("status"),
            result=statuses,
            artifact_manifest=str(artifacts_path),
            finished_at=now(), observed_at=now(),
            peak_rss_kib=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        )
    except BaseException as exc:
        stop_telemetry.set()
        telemetry_thread.join()
        receipt.update(status="execution_failed", finished_at=now(), observed_at=now(),
                       error={"type": type(exc).__name__, "message": str(exc),
                              "traceback": traceback.format_exc()},
                       peak_rss_kib=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        persist()
        raise
    stop_telemetry.set()
    telemetry_thread.join()
    persist()
    return receipt


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("preflight", "run"))
    parser.add_argument("--config")
    parser.add_argument("--iteration", default=DEFAULT_ITERATION)
    parser.add_argument("--producer-iteration", default=DEFAULT_ITERATION)
    parser.add_argument("--execution-location", default=socket.gethostname())
    parser.add_argument("--qz-notebook")
    args = parser.parse_args(argv)
    repo = Path(__file__).resolve().parents[1]
    config_path = Path(args.config) if args.config else repo / EXPECTED_CONFIG
    if args.command == "preflight":
        value = validate_project(config_path, repo, require_inputs=True,
                                 iteration=args.iteration,
                                 producer_iteration=args.producer_iteration)
        print(json.dumps({key: str(item) if isinstance(item, Path) else item
                          for key, item in value.items() if key != "config"},
                         ensure_ascii=False, indent=2))
        return 0
    receipt = execute(args)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0 if receipt["status"] == "execution_complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
