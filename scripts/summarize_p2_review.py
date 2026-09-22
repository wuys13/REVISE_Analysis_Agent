"""Read registered Impact evidence; never recompute or grant scientific acceptance."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite_json(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: finite_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [finite_json(item) for item in value]
    return value


def summarize(result_dir: Path) -> dict[str, Any]:
    result_dir = result_dir.resolve()
    source = result_dir / "result.json"
    result = json.loads(source.read_text())
    registered = set(result["outputs"].values())
    evidence: dict[str, Any] = {}

    def read(relative: str) -> Any:
        path = (result_dir / relative).resolve()
        if not path.is_relative_to(result_dir):
            raise ValueError(f"Evidence escapes result directory: {relative}")
        if relative not in registered:
            evidence[relative] = {"status": "not_registered"}
            return None
        if not path.is_file():
            evidence[relative] = {"status": "missing_registered_file"}
            return None
        evidence[relative] = {"status": "read", "sha256": sha256(path)}
        if path.suffix == ".json":
            return json.loads(path.read_text())
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle))

    overview = read("tables/input_overview.csv")
    partitions = read("tables/partition_summary.csv")
    label_summary = read("tables/reconstruction_label_summary.csv")
    cohorts = read("tables/cohort_manifest.json")
    scopes = {}
    for scope in result["parameters"]["scopes"]:
        slug = "".join(char if char.isalnum() else "_" for char in scope)
        support = read(f"tables/window_scale_support_{slug}.csv")
        configured = result["parameters"]["parent_window_side_microns"][scope]
        selected = None if support is None else [
            row for row in support
            if math.isclose(float(row["window_side_microns"]), configured)
        ]
        membership = read(f"tables/membership_summary_{slug}.json")
        scopes[scope] = {
            "configured_window_microns": configured,
            "support_at_configured_scale": selected,
            "support_recommendation": read(f"tables/window_scale_recommendation_{slug}.json"),
            "state_threshold": read(f"tables/state_{slug}_threshold.json"),
            "gain_threshold": read(f"tables/gain_{slug}_threshold.json"),
            "shared_gene_moran": read(f"tables/moran_shared_genes_{slug}_metadata.json"),
            "membership_summary": None if membership is None else membership.get("summary"),
        }

    anatomy = {}
    for side in ("raw", "svc"):
        points = read(f"tables/{side}_point_anatomy.csv")
        anatomy[side] = None if points is None else {
            "observed_units": len(points),
            "counts": dict(Counter(row["anatomy_region"] for row in points)),
            "interpretation": "SVC-defined context; counts use this side's observed units",
        }

    pathways = {}
    for relative in sorted(registered):
        if not relative.startswith("tables/pathway_") or not relative.endswith("_metadata.json"):
            continue
        metadata = read(relative)
        if metadata is None:
            pathways[relative] = None
            continue
        pathways[relative] = {
            "coverage": {key: item for key, item in metadata["coverage"].items()
                         if not isinstance(item, list)},
            "provider": metadata["metadata"],
            "n_native_genes": len(metadata["gene_axis"]),
            "n_eligible_units": metadata["n_eligible_units"],
            "n_scored_units": metadata["n_scored_units"],
            "expression": metadata["expression"],
        }

    return finite_json({
        "schema_version": 1,
        "kind": "read_only_impact_review_summary",
        "observed_at": datetime.now(timezone.utc).isoformat(),
        "source_result": str(source),
        "source_result_sha256": sha256(source),
        "sample_id": result["sample_id"],
        "analysis": result["analysis"],
        "result_status": result["status"],
        "scientific_acceptance": False,
        "parameters": result["parameters"],
        "stages": [{key: stage.get(key) for key in ("stage", "status")}
                   for stage in result["stages"]],
        "stage_errors": result["stage_errors"],
        "unavailable": result["unavailable"],
        "input_overview": overview,
        "partitions": partitions,
        "reconstruction_labels": label_summary,
        "cohorts": cohorts,
        "scopes": scopes,
        "anatomy": anatomy,
        "pathways": pathways,
        "evidence": evidence,
        "limitations": [
            "This summary is an evidence index, not an independent biological validation.",
            "Native AUCell scores do not establish cross-side improvement when backgrounds differ.",
            "Membership is partition correspondence, not a change of biological identity.",
            "Non-finite source JSON values are preserved as null, never zero.",
            "Unregistered or missing evidence is explicit; consult result.unavailable and stage_errors.",
        ],
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.resolve().is_relative_to(args.result_dir.resolve()):
        parser.error("Place the review summary outside the immutable result directory")
    summary = summarize(args.result_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
