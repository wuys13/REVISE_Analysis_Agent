"""Staged reconstruction-impact analysis.

Supplied SVC reconstruction labels are the scientific object of interest.
Expression-derived partitions are optional baselines and never replace those
labels or become prerequisites for native State analysis.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any
import traceback

import numpy as np
import pandas as pd

from ._shared import DEFAULT_SCOPES, coordinates, deterministic_subset, save_json, save_table, unavailable
from revise_analysis.methods.errors import PrerequisiteUnavailable

from .impact_tables import changed_unit_summary, common_valid_delta, label_composition


STAGE_ORDER = ("input", "baseline", "support", "diversity", "regions", "anatomy",
               "molecular", "membership", "integration", "figures")
SECTION_DEFINITIONS = (
    ("input", "输入与重建标签", "输入能力和已提供的重建状态标签。"),
    ("baseline", "表达 baseline", "Raw 独立表达分群和已提供的 Raw Level2 标签。"),
    ("support", "空间尺度与支持", "独立于最终结果的支持曲线和尺度推荐。"),
    ("diversity", "局部标签多样性", "原生重建标签多样性及可用的 baseline 字段。"),
    ("regions", "State 与差值候选", "连续场、阈值诊断和候选掩膜。"),
    ("anatomy", "Anatomy 背景", "完整 Raw Anatomy 窗口及按观测点组成汇总的 parent 窗口。"),
    ("molecular", "分子空间结果", "可用时保存原生 Moran 结果和固定 cohort 的 AUCell 分数。"),
    ("membership", "共同 ID 的 membership", "仅在明确共同 ID cohort 上生成的单位变化证据。"),
    ("integration", "集成空间证据", "带有明确支持数和分母的派生事实表。"),
)


def _slug(value: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in str(value)).strip("_") or "scope"


def effective_parameters(sample: Any, overrides: dict | None = None) -> dict:
    """Resolve formal sample defaults and explicit overrides once."""
    allowed = {
        "scopes", "resolution", "n_top_genes", "sample_n_units", "random_state",
        "membership_same_units", "raw_k_control", "matched_k_resolutions",
        "anatomy_window_side_microns", "parent_window_side_microns", "window_side_microns",
        "window_scale_candidates_microns", "min_window_units", "n_window_draws",
        "region_n_bootstrap", "gene_sets", "geneset_path", "gene_set_names",
        "pathway_auc_threshold", "moran_n_neighbors", "anatomy_tumor_label",
        "anatomy_normal_label", "svc_leiden_baseline",
    }
    supplied = dict(sample.config.get("analysis_parameters", {}).get("reconstruction_impact", {}) or {})
    supplied.update(dict(overrides or {}))
    if unknown := sorted(set(supplied) - allowed):
        raise ValueError(f"Unsupported analysis parameters: {', '.join(unknown)}")
    scopes_value = supplied.get("scopes", DEFAULT_SCOPES)
    if isinstance(scopes_value, (str, bytes)) or not isinstance(scopes_value, (list, tuple)):
        raise ValueError("scopes must be a sequence of names, not a scalar")
    scopes = list(scopes_value)
    if (not scopes or any(not isinstance(x, str) or not x.strip() for x in scopes)
            or len(scopes) != len(set(scopes)) or len(scopes) != len({_slug(x) for x in scopes})):
        raise ValueError("scopes collide or are not unique non-empty names after artifact slug conversion")
    old_window = supplied.get("window_side_microns")
    parent = supplied.get("parent_window_side_microns", old_window if old_window is not None else 40.0)
    for key in ("membership_same_units", "raw_k_control", "svc_leiden_baseline"):
        if key in supplied and type(supplied[key]) is not bool:
            raise ValueError(f"{key} must be a boolean")
    for key, minimum in (("n_top_genes", 1), ("random_state", 0), ("min_window_units", 1),
                         ("n_window_draws", 1), ("region_n_bootstrap", 1), ("moran_n_neighbors", 1)):
        if key in supplied and (type(supplied[key]) is not int or supplied[key] < minimum):
            raise ValueError(f"{key} must be an integer >= {minimum}")
    if ("sample_n_units" in supplied and supplied["sample_n_units"] is not None
            and (type(supplied["sample_n_units"]) is not int or supplied["sample_n_units"] < 1)):
        raise ValueError("sample_n_units must be a positive integer or null")
    for key in ("resolution", "anatomy_window_side_microns", "window_side_microns", "pathway_auc_threshold"):
        if key in supplied and (type(supplied[key]) not in (int, float) or isinstance(supplied[key], bool)
                                or not np.isfinite(float(supplied[key]))):
            raise ValueError(f"{key} must be finite numeric")
    for key in ("matched_k_resolutions", "window_scale_candidates_microns"):
        if key in supplied:
            values = supplied[key]
            if (isinstance(values, (str, bytes)) or not isinstance(values, (list, tuple))
                    or any(type(x) not in (int, float) or isinstance(x, bool) or not np.isfinite(float(x)) for x in values)):
                raise ValueError(f"{key} must be a finite numeric sequence")
    parent_values = list(parent.values()) if isinstance(parent, dict) else [parent]
    if any(type(v) not in (int, float) or isinstance(v, bool) or not np.isfinite(float(v)) for v in parent_values):
        raise ValueError("parent_window_side_microns must contain finite numeric values")
    parent = ({str(k): float(v) for k, v in parent.items()} if isinstance(parent, dict)
              else {scope: float(parent) for scope in scopes})
    if missing := [scope for scope in scopes if scope not in parent]:
        raise ValueError(f"parent_window_side_microns missing scopes: {missing}")
    result = {
        "scopes": scopes, "resolution": float(supplied.get("resolution", .5)),
        "n_top_genes": int(supplied.get("n_top_genes", 2000)),
        "sample_n_units": supplied.get("sample_n_units"),
        "random_state": int(supplied.get("random_state", 42)),
        "membership_same_units": bool(supplied.get("membership_same_units", False)),
        "raw_k_control": bool(supplied.get("raw_k_control", False)),
        "matched_k_resolutions": sorted({float(x) for x in supplied.get("matched_k_resolutions", (.2, .4, .6, .8, 1.0))}),
        "anatomy_window_side_microns": float(supplied.get("anatomy_window_side_microns", 40.0)),
        "parent_window_side_microns": parent,
        "window_scale_candidates_microns": sorted({float(x) for x in supplied.get("window_scale_candidates_microns", (16, 24, 32, 40, 56, 80))}),
        "min_window_units": int(supplied.get("min_window_units", 4)),
        "n_window_draws": int(supplied.get("n_window_draws", 200)),
        "region_n_bootstrap": int(supplied.get("region_n_bootstrap", 500)),
        "gene_sets": supplied.get("gene_sets"), "geneset_path": supplied.get("geneset_path"),
        "gene_set_names": supplied.get("gene_set_names"),
        "pathway_auc_threshold": float(supplied.get("pathway_auc_threshold", .05)),
        "moran_n_neighbors": int(supplied.get("moran_n_neighbors", 6)),
        "anatomy_tumor_label": str(supplied.get("anatomy_tumor_label", "Tumor")),
        "anatomy_normal_label": str(supplied.get("anatomy_normal_label", "Intestinal Epithelial")),
        "svc_leiden_baseline": bool(supplied.get("svc_leiden_baseline", False)),
    }
    if result["geneset_path"]:
        path = Path(result["geneset_path"])
        source = getattr(sample, "source", None)
        if not path.is_absolute() and source is not None:
            result["geneset_path"] = str((Path(source).parent / path).resolve())
    positive = ("resolution", "n_top_genes", "anatomy_window_side_microns", "min_window_units",
                "n_window_draws", "region_n_bootstrap", "moran_n_neighbors")
    if (any(not np.isfinite(float(result[k])) or result[k] <= 0 for k in positive)
            or any(not np.isfinite(v) or v <= 0 for v in parent.values())):
        raise ValueError("partition, anatomy, parent-window and sampling parameters must be positive")
    if (not result["matched_k_resolutions"] or not result["window_scale_candidates_microns"]
            or any(x <= 0 for x in (*result["matched_k_resolutions"], *result["window_scale_candidates_microns"]))):
        raise ValueError("resolution and window-scale candidate lists must be nonempty and positive")
    return result


class ImpactWorkflow:
    """Ordered workflow shared by notebook and batch execution."""

    def __init__(self, sample: Any, output_dir: Path, parameters: dict | None = None,
                 continue_on_error: bool = False):
        self.sample, self.output_dir = sample, Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parameters = effective_parameters(sample, parameters)
        self.continue_on_error = continue_on_error
        self.outputs: dict[str, str] = {}
        self.missing: list[dict] = []
        self.stage_errors: list[dict] = []
        self.stage_records: list[dict] = []
        self.stage_artifacts: dict[str, list[str]] = {name: [] for name in STAGE_ORDER}
        self.state: dict[str, Any] = {}

    def run(self) -> dict:
        for name in STAGE_ORDER:
            self.run_stage(name)
        return self.result()

    def run_stage(self, name: str) -> dict:
        """Run one named stage for notebook use and return its stage record."""
        if name not in STAGE_ORDER:
            raise ValueError(f"Unknown impact stage {name!r}")
        before = set(self.outputs)
        missing_before = len(self.missing)
        try:
            getattr(self, f"stage_{name}")()
            status = "completed"
        except Exception as exc:
            if not self.continue_on_error:
                raise
            status = "error"
            record = {"stage": name, "error_type": type(exc).__name__, "message": str(exc),
                      "traceback": traceback.format_exc()}
            self.stage_errors.append(record)
            relative = f"errors/{name}.txt"
            path = self.output_dir / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(record["traceback"], encoding="utf-8")
            self.outputs[relative] = relative
        artifacts = sorted(set(self.outputs) - before)
        if status == "completed" and len(self.missing) > missing_before:
            status = "partial" if artifacts else "unavailable"
        self.stage_artifacts[name].extend(artifacts)
        record = {"stage": name, "status": status, "artifacts": artifacts}
        self.stage_records.append(record)
        return record

    def result(self) -> dict:
        """Return the current notebook/batch result without running more stages."""
        reading_artifacts = {key: list(value) for key, value in self.stage_artifacts.items()}
        for artifact in ("tables/raw_anatomy_context.csv", "tables/raw_anatomy_windows.csv"):
            if artifact in reading_artifacts["support"]:
                reading_artifacts["support"].remove(artifact)
                if artifact not in reading_artifacts["anatomy"]:
                    reading_artifacts["anatomy"].insert(0, artifact)
        sections = [{"id": i, "title": t, "description": d,
                     "artifacts": reading_artifacts.get(i, [])}
                    for i, t, d in SECTION_DEFINITIONS]
        status = ("failed" if self.stage_errors else
                  ("partial" if self.missing else ("succeeded" if self.outputs else "skipped")))
        return {"status": status, "parameters": self.parameters, "outputs": self.outputs,
                "unavailable": self.missing, "stage_errors": self.stage_errors,
                "stages": self.stage_records, "sections": sections}

    def _save_table(self, table, relative):
        save_table(table, self.output_dir, relative, self.outputs)

    def _save_json(self, value, relative):
        save_json(value, self.output_dir, relative, self.outputs)

    def _scope(self, side: str, scope: str, *, adata=None):
        data = getattr(self.sample, side) if adata is None else adata
        if scope == "All":
            return data
        if self.sample.broad_key not in data.obs:
            raise PrerequisiteUnavailable(f"{side}: broad-label column is unavailable for {scope!r}")
        broad = self.sample.labels(side, self.sample.broad_key).reindex(data.obs_names)
        ids = broad.index[broad == scope]
        if ids.empty:
            raise PrerequisiteUnavailable(f"no {scope!r} units on {side}")
        return data[ids]

    def _expression_reason(self, side):
        method = getattr(self.sample, "expression_unavailable", None)
        return method(side) if callable(method) else None

    def _expression(self, side):
        method = getattr(self.sample, "expression", None)
        if callable(method):
            return method(side)
        configured = self.sample.config.get("expression", {})
        return configured.get(side, configured) if isinstance(configured, dict) else {"identity": "unknown", "scale": "unknown"}

    def _micron_scale(self):
        value = self.sample.config.get("spatial", {}).get("microns_per_coordinate")
        if value is None or not np.isfinite(float(value)) or float(value) <= 0:
            raise ValueError("spatial.microns_per_coordinate is required for micrometer windows")
        return float(value)

    def _parent_coords(self, side, scope):
        return coordinates(self._scope(side, scope), self.sample.spatial_key)

    def _svc_parent(self, scope):
        """Return the requested SVC parent restricted to valid main labels."""
        cohort = self._scope("svc", scope)
        labels = self.state.get("svc_labels")
        if labels is None:
            raise ValueError("SVC reconstruction labels unavailable")
        ids = cohort.obs_names.intersection(labels.index)
        if ids.empty:
            raise ValueError(f"no valid reconstructed labels in SVC scope {scope!r}")
        return cohort[ids]

    def stage_input(self):
        source = Path(self.sample.source).resolve()
        files = self.sample.config.get("files", {})
        resolved_files = {}
        for side in ("raw", "svc"):
            value = files.get(side)
            if value is not None:
                path = Path(value).expanduser()
                resolved_files[side] = str((source.parent / path).resolve() if not path.is_absolute() else path.resolve())
        self._save_json({
            "sample_id": self.sample.sample_id,
            "sample_config": str(source),
            "resolved_files": resolved_files,
            "columns": {
                "broad": self.sample.broad_key, "subtype": self.sample.subtype_key,
                "reconstruction": getattr(self.sample, "reconstruction_key", "SVC_cluster"),
            },
            "expression": {side: self._expression(side) for side in ("raw", "svc")},
            "spatial": self.sample.config.get("spatial", {}),
            "metadata": self.sample.config.get("metadata", {}),
            "provenance": self.sample.config.get("provenance", {}),
            "effective_parameters": self.parameters,
        }, "tables/input_manifest.json")
        key = getattr(self.sample, "reconstruction_key", "SVC_cluster")
        if key not in self.sample.svc.obs:
            unavailable("svc_reconstruction_labels", f"SVC reconstruction label column {key!r} is unavailable", self.missing)
            self.state["svc_labels"] = None
            return
        # Reconstruction labels are an input scientific identity.  Broad-label
        # aliases must not silently rewrite their values.
        labels = self.sample.svc.obs[key].copy()
        labels.index = self.sample.svc.obs_names
        if labels.empty or labels.notna().sum() == 0:
            unavailable("svc_reconstruction_labels", f"SVC reconstruction label column {key!r} has no valid labels", self.missing)
            self.state["svc_labels"] = None
            return
        if labels.isna().any():
            unavailable("svc_reconstruction_labels", f"{int(labels.isna().sum())} SVC units have missing reconstruction labels and are excluded", self.missing)
        labels = labels.dropna().astype(str).rename("reconstruction_label")
        self.state["svc_labels"] = labels
        rows = []
        for side in ("raw", "svc"):
            expression = self._expression(side)
            rows.append({"side": side, "n_units": int(getattr(self.sample, side).n_obs),
                         "n_genes": int(getattr(self.sample, side).n_vars),
                         "expression_identity": expression.get("identity", "unknown"),
                         "transformation_state": expression.get("scale", "unknown"),
                         "expression_available": self._expression_reason(side) is None})
        self._save_table(pd.DataFrame(rows), "tables/input_overview.csv")
        counts = labels.value_counts().rename_axis("reconstruction_label").rename("n_units").reset_index()
        counts["label_role"] = "reconstructed_state_label"
        self._save_table(counts, "tables/reconstruction_label_summary.csv")
        try:
            units = coordinates(self.sample.svc, self.sample.spatial_key).reindex(labels.index)
        except (KeyError, ValueError) as exc:
            units = pd.DataFrame(index=labels.index, columns=["x", "y"], dtype=float)
            unavailable("svc_reconstruction_coordinates", str(exc), self.missing)
        units["reconstruction_label"] = labels
        try:
            units["broad_label"] = self.sample.labels("svc", self.sample.broad_key).reindex(labels.index)
        except KeyError:
            units["broad_label"] = pd.NA
            unavailable("svc_broad_labels", f"SVC broad-label column {self.sample.broad_key!r} is unavailable", self.missing)
        units["label_role"] = "reconstructed_state_label"
        self._save_table(units, "tables/reconstruction_labels_units.csv")

    def stage_baseline(self):
        partitions, summaries, prepared_graphs, partition_cohorts = {}, [], {}, {}
        self.state["partitions"] = partitions
        self.state["prepared_graphs"], self.state["partition_cohorts"] = prepared_graphs, partition_cohorts
        for side in ("raw", "svc"):
            if side == "svc" and not self.parameters["svc_leiden_baseline"]:
                continue
            if reason := self._expression_reason(side):
                unavailable(f"partition:{side}", reason, self.missing)
                continue
            sampled = deterministic_subset(getattr(self.sample, side), self.parameters["sample_n_units"],
                                           self.parameters["random_state"] + (side == "svc"))
            for scope in self.parameters["scopes"]:
                try:
                    cohort = self._scope(side, scope, adata=sampled)
                    result = compute_partitions(
                        cohort, resolution=self.parameters["resolution"], n_top_genes=self.parameters["n_top_genes"],
                        random_state=self.parameters["random_state"],
                        transformation_state=self._expression(side).get("scale", "unknown"))
                    labels = result["labels"].astype(str).rename("partition")
                    partitions[side, scope] = labels
                    partition_cohorts[side, scope] = cohort
                    if result.get("_prepared") is not None:
                        prepared_graphs[side, scope] = result["_prepared"]
                    self._save_table(labels.to_frame().assign(side=side, scope=scope),
                                     f"tables/partitions_{side}_{_slug(scope)}.csv")
                    summaries.append({"side": side, "scope": scope, **result["summary"]})
                except (ImportError, PrerequisiteUnavailable) as exc:
                    unavailable(f"partition:{side}:{scope}", str(exc), self.missing)
        self.state["partitions"] = partitions
        self.state["prepared_graphs"], self.state["partition_cohorts"] = prepared_graphs, partition_cohorts
        if summaries:
            self._save_table(pd.DataFrame(summaries), "tables/partition_summary.csv")
        key = self.sample.subtype_key
        if key not in self.sample.raw.obs:
            unavailable("raw_level2_baseline", f"Raw column {key!r} is unavailable", self.missing)
        else:
            labels = self.sample.labels("raw", key)
            if labels.isna().any():
                unavailable("raw_level2_baseline", f"Raw column {key!r} contains missing values", self.missing)
            else:
                self.state["raw_level2"] = labels.astype(str)
                self._save_table(labels.rename("raw_level2"), "tables/raw_level2_baseline.csv")
                self._save_table(labels.value_counts().rename_axis("raw_level2").rename("n_units").reset_index(),
                                 "tables/raw_level2_baseline_summary.csv")

    def stage_support(self):
        from revise_analysis.methods.regions import assign_anatomy_candidates, assign_square_windows, select_window_scale
        try:
            scale = self._micron_scale()
            full_raw = coordinates(self.sample.raw, self.sample.spatial_key)
        except (KeyError, ValueError) as exc:
            unavailable("spatial_geometry", str(exc), self.missing)
            return
        origin = (float(full_raw.x.min()), float(full_raw.y.min()))
        anatomy_side = self.parameters["anatomy_window_side_microns"] / scale
        self.state.update(origin=origin, microns_per_coordinate=scale, anatomy_side=anatomy_side,
                          anatomy_windows=None, raw_anatomy_points=None)
        try:
            anatomy_assignments = assign_square_windows(full_raw, window_side_length=anatomy_side, origin=origin)
            raw_broad = self.sample.labels("raw", self.sample.broad_key)
            anatomy_windows = assign_anatomy_candidates(
                anatomy_assignments, raw_broad, tumor_label=self.parameters["anatomy_tumor_label"],
                normal_source_label=self.parameters["anatomy_normal_label"])
            anatomy_points = full_raw.join(anatomy_assignments[["window_id"]]).rename(columns={"window_id": "anatomy_window_id"})
            anatomy_points = anatomy_points.join(anatomy_windows.set_index("window_id")[["level1_region"]], on="anatomy_window_id")
            anatomy_points["broad_label"] = raw_broad.reindex(anatomy_points.index)
            self.state.update(anatomy_windows=anatomy_windows, raw_anatomy_points=anatomy_points)
            self._save_table(anatomy_windows, "tables/raw_anatomy_windows.csv")
            self._save_table(anatomy_points, "tables/raw_anatomy_context.csv")
        except (KeyError, ValueError) as exc:
            unavailable("raw_anatomy_context", str(exc), self.missing)
        if self.state.get("svc_labels") is None:
            unavailable("window_support", "SVC reconstruction labels unavailable", self.missing)
            return
        candidates = [x / scale for x in self.parameters["window_scale_candidates_microns"]]
        for scope in self.parameters["scopes"]:
            try:
                svc_coords = coordinates(self._svc_parent(scope), self.sample.spatial_key)
            except (KeyError, ValueError) as exc:
                unavailable(f"window_support:svc:{scope}", str(exc), self.missing)
                continue
            recommendation, svc_curve = select_window_scale(svc_coords, candidate_window_sides=candidates,
                min_parent_units=self.parameters["min_window_units"], origin=origin)
            svc_curve["window_side_microns"], svc_curve["support_kind"] = svc_curve["window_side_length"] * scale, "svc_parent"
            raw_labels = self.state.get("partitions", {}).get(("raw", scope))
            if raw_labels is not None:
                raw_coords = coordinates(self.sample.raw, self.sample.spatial_key).reindex(raw_labels.index)
            else:
                try:
                    raw_coords = self._parent_coords("raw", scope)
                except (KeyError, ValueError):
                    raw_coords = coordinates(self.sample.raw, self.sample.spatial_key).iloc[0:0]
            _, raw_curve = select_window_scale(raw_coords, candidate_window_sides=candidates,
                min_parent_units=self.parameters["min_window_units"], origin=origin)
            raw_curve["window_side_microns"] = raw_curve["window_side_length"] * scale
            raw_curve["support_kind"] = "raw_baseline" if raw_labels is not None else "raw_coordinate_potential"
            common_rows = []
            for coord_side, micron_side in zip(candidates, self.parameters["window_scale_candidates_microns"], strict=True):
                raw_counts = assign_square_windows(raw_coords, window_side_length=coord_side, origin=origin).groupby("window_id").size()
                svc_counts = assign_square_windows(svc_coords, window_side_length=coord_side, origin=origin).groupby("window_id").size()
                raw_valid = set(raw_counts.index[raw_counts >= self.parameters["min_window_units"]])
                svc_valid = set(svc_counts.index[svc_counts >= self.parameters["min_window_units"]])
                common_ids = raw_valid & svc_valid
                common_rows.append({"window_side_length": coord_side, "window_side_microns": micron_side,
                    "support_kind": "common_valid" if raw_labels is not None else "common_coordinate_potential",
                    "n_valid_windows": len(common_ids), "n_svc_valid_windows": len(svc_valid),
                    "common_fraction_of_svc_valid": len(common_ids) / len(svc_valid) if svc_valid else np.nan,
                    "raw_units_in_common_valid": int(raw_counts.reindex(list(common_ids), fill_value=0).sum()),
                    "svc_units_in_common_valid": int(svc_counts.reindex(list(common_ids), fill_value=0).sum()),
                    "raw_common_retained_fraction": float(raw_counts.reindex(list(common_ids), fill_value=0).sum() / len(raw_coords)) if len(raw_coords) else np.nan,
                    "svc_common_retained_fraction": float(svc_counts.reindex(list(common_ids), fill_value=0).sum() / len(svc_coords)) if len(svc_coords) else np.nan})
            self._save_table(pd.concat([svc_curve, raw_curve, pd.DataFrame(common_rows)], ignore_index=True, sort=False),
                             f"tables/window_scale_support_{_slug(scope)}.csv")
            recommendation = {**recommendation,
                "window_side_microns": None if recommendation.get("window_side_length") is None else recommendation["window_side_length"] * scale,
                "configured_window_side_microns": self.parameters["parent_window_side_microns"][scope],
                "configured_value_overridden": False, "basis": "svc_parent_support_only"}
            self._save_json(recommendation, f"tables/window_scale_recommendation_{_slug(scope)}.json")

    def stage_diversity(self):
        from revise_analysis.methods.regions import compute_window_diversity
        windows = {}
        if self.state.get("svc_labels") is None or "origin" not in self.state:
            unavailable("window_diversity", "SVC reconstruction labels or spatial geometry unavailable", self.missing)
            self.state["windows"] = windows
            return
        for scope in self.parameters["scopes"]:
            side_length = self.parameters["parent_window_side_microns"][scope] / self.state["microns_per_coordinate"]
            try:
                svc = self._svc_parent(scope)
            except (KeyError, ValueError) as exc:
                unavailable(f"window_diversity:svc:{scope}", str(exc), self.missing)
                continue
            svc_labels = self.state["svc_labels"].reindex(svc.obs_names)
            table = compute_window_diversity(coordinates(svc, self.sample.spatial_key), svc_labels,
                window_side_length=side_length, origin=self.state["origin"], min_units=self.parameters["min_window_units"],
                n_draws=self.parameters["n_window_draws"], random_state=self.parameters["random_state"] + 1)
            table["label_system"] = "svc_reconstruction_label"
            windows["svc", scope] = table
            self._save_table(table, f"tables/window_diversity_svc_{_slug(scope)}.csv")
            raw_partition = self.state.get("partitions", {}).get(("raw", scope))
            if raw_partition is not None:
                raw_table = compute_window_diversity(coordinates(self.sample.raw, self.sample.spatial_key).reindex(raw_partition.index), raw_partition,
                    window_side_length=side_length, origin=self.state["origin"], min_units=self.parameters["min_window_units"],
                    n_draws=self.parameters["n_window_draws"], random_state=self.parameters["random_state"])
                raw_table["label_system"] = "raw_leiden"
                windows["raw", scope] = raw_table
                self._save_table(raw_table, f"tables/window_diversity_raw_{_slug(scope)}.csv")
            svc_partition = self.state.get("partitions", {}).get(("svc", scope))
            if svc_partition is not None:
                svc_baseline = compute_window_diversity(
                    coordinates(self.sample.svc, self.sample.spatial_key).reindex(svc_partition.index), svc_partition,
                    window_side_length=side_length, origin=self.state["origin"],
                    min_units=self.parameters["min_window_units"], n_draws=self.parameters["n_window_draws"],
                    random_state=self.parameters["random_state"] + 31)
                svc_baseline["label_system"] = "svc_leiden_optional_baseline"
                windows["svc_leiden", scope] = svc_baseline
                self._save_table(svc_baseline, f"tables/window_diversity_svc_leiden_baseline_{_slug(scope)}.csv")
            if "raw_level2" in self.state:
                raw = self._scope("raw", scope)
                labels = self.state["raw_level2"].reindex(raw.obs_names)
                level2 = compute_window_diversity(coordinates(raw, self.sample.spatial_key), labels,
                    window_side_length=side_length, origin=self.state["origin"], min_units=self.parameters["min_window_units"],
                    n_draws=self.parameters["n_window_draws"], random_state=self.parameters["random_state"] + 17)
                level2["label_system"] = "raw_level2"
                windows["raw_level2", scope] = level2
                self._save_table(level2, f"tables/window_diversity_raw_level2_baseline_{_slug(scope)}.csv")
        self.state["windows"] = windows

    def _write_region(self, table, *, name, value_column, side_microns, positive_only=False):
        from revise_analysis.methods.regions import flag_region_windows, select_region_threshold
        values = pd.to_numeric(table[value_column], errors="coerce")
        info, audit = select_region_threshold(values[values > 0] if positive_only else values,
            n_bootstrap=self.parameters["region_n_bootstrap"], random_state=self.parameters["random_state"])
        self._save_json({**info, "window_side_microns": side_microns}, f"tables/{name}_threshold.json")
        if not audit.empty:
            self._save_table(audit, f"tables/{name}_threshold_bootstrap.csv")
        if info.get("threshold") is None:
            flagged = table.copy()
            flagged["in_region"], flagged["region_status"] = pd.NA, info.get("status", "no_stable_threshold")
            unavailable(name, f"stable region threshold unavailable: {info.get('status')}", self.missing)
        else:
            summary, flagged = flag_region_windows(table, threshold=float(info["threshold"]),
                value_column=value_column, window_side_length=side_microns)
            summary["window_side_microns"] = side_microns
            self._save_table(summary, f"tables/{name}_region_extent.csv")
            valid = flagged["valid_window"].astype(bool)
            flagged["in_region"] = flagged["in_region"].astype("boolean")
            flagged.loc[~valid, "in_region"] = pd.NA
            flagged["region_status"] = np.where(valid, "ok", "insufficient_support")
        self._save_table(flagged, f"tables/{name}_region_windows.csv")
        return flagged

    def stage_regions(self):
        region_tables, gains = {}, {}
        for scope in self.parameters["scopes"]:
            side_microns = self.parameters["parent_window_side_microns"][scope]
            if ("svc", scope) not in self.state["windows"]:
                unavailable(f"state:{scope}", "SVC parent local diversity unavailable", self.missing)
                continue
            state = self._write_region(self.state["windows"]["svc", scope], name=f"state_{_slug(scope)}",
                value_column="neff", side_microns=side_microns)
            region_tables["state", scope] = state
            raw = self.state["windows"].get(("raw", scope))
            if raw is not None:
                gain = common_valid_delta(raw, self.state["windows"]["svc", scope], delta_column="delta_neff_vs_raw_leiden")
                if gain.empty:
                    unavailable(f"gain:{scope}", "no common valid spatial windows", self.missing)
                else:
                    gain["interpretation"] = "delta_vs_raw_leiden_gain_candidate"
                    self._save_table(gain, f"tables/gain_{_slug(scope)}_common_valid_windows.csv")
                    gains[scope] = self._write_region(gain, name=f"gain_{_slug(scope)}",
                        value_column="delta_neff_vs_raw_leiden", side_microns=side_microns, positive_only=True)
            else:
                unavailable(f"gain:{scope}", "Raw Leiden baseline unavailable", self.missing)
            level2 = self.state["windows"].get(("raw_level2", scope))
            if level2 is not None:
                comparison = common_valid_delta(level2, self.state["windows"]["svc", scope], delta_column="descriptive_delta_neff")
                if not comparison.empty:
                    comparison = comparison.rename(columns={"raw_neff": "raw_level2_neff", "svc_neff": "svc_reconstruction_neff"})
                    self._save_table(comparison, f"tables/raw_level2_baseline_vs_svc_{_slug(scope)}_common_valid_windows.csv")
        self.state["region_tables"], self.state["gains"] = region_tables, gains
        if self.parameters["raw_k_control"]:
            self._raw_k_control()

    def _raw_k_control(self):
        from revise_analysis.methods.errors import PrerequisiteUnavailable
        from revise_analysis.methods.partition import select_matched_k_partition
        from revise_analysis.methods.regions import compute_window_diversity
        if reason := self._expression_reason("raw"):
            unavailable("raw_k_control", reason, self.missing)
            return
        controls = {}
        for scope in self.parameters["scopes"]:
            try:
                raw = self.state.get("partition_cohorts", {}).get(("raw", scope))
                prepared = self.state.get("prepared_graphs", {}).get(("raw", scope))
                if raw is None or prepared is None:
                    unavailable(f"raw_k_control:{scope}", "main Raw baseline graph/cohort unavailable", self.missing)
                    continue
                target = int(self.state["svc_labels"].reindex(self._svc_parent(scope).obs_names).nunique())
                result = select_matched_k_partition(prepared, candidate_resolutions=self.parameters["matched_k_resolutions"],
                    target_k=target, main_resolution=self.parameters["resolution"])
                self._save_table(result["candidates"], f"tables/raw_k_control_sweep_{_slug(scope)}.csv")
                self._save_json(result["selection"], f"tables/raw_k_control_selection_{_slug(scope)}.json")
                selected = result["selected"]
                labels = selected["labels"]
                side_length = self.parameters["parent_window_side_microns"][scope] / self.state["microns_per_coordinate"]
                table = compute_window_diversity(coordinates(raw, self.sample.spatial_key).reindex(labels.index), labels,
                    window_side_length=side_length, origin=self.state["origin"], min_units=self.parameters["min_window_units"],
                    n_draws=self.parameters["n_window_draws"], random_state=self.parameters["random_state"])
                self._save_table(table, f"tables/window_diversity_raw_k_control_{_slug(scope)}.csv")
                delta = common_valid_delta(table, self.state["windows"]["svc", scope], delta_column="delta_neff_vs_raw_k_control")
                if not delta.empty:
                    delta["control_status"] = result["selection"]["status"]
                    self._save_table(delta, f"tables/raw_k_control_delta_{_slug(scope)}.csv")
                controls[scope] = labels
            except (ImportError, PrerequisiteUnavailable) as exc:
                unavailable(f"raw_k_control:{scope}", str(exc), self.missing)
        self.state["raw_k_controls"] = controls

    def stage_anatomy(self):
        from revise_analysis.methods.regions import aggregate_parent_anatomy, assign_points_to_anatomy, assign_square_windows
        point_anatomy, parent_anatomy = {}, {}
        if "origin" not in self.state:
            unavailable("anatomy_mapping", "spatial geometry unavailable", self.missing)
            self.state["point_anatomy"], self.state["parent_anatomy"] = point_anatomy, parent_anatomy
            return
        for side in ("raw", "svc"):
            side_coords = coordinates(getattr(self.sample, side), self.sample.spatial_key)
            if self.state.get("anatomy_windows") is None:
                point_anatomy[side] = side_coords.copy()
                point_anatomy[side]["anatomy_window_id"] = pd.NA
                point_anatomy[side]["anatomy_covered"] = False
                point_anatomy[side]["anatomy_region"] = "Unknown"
            else:
                point_anatomy[side] = assign_points_to_anatomy(side_coords,
                    self.state["anatomy_windows"], anatomy_window_side_length=self.state["anatomy_side"], origin=self.state["origin"])
            self._save_table(point_anatomy[side], f"tables/{side}_point_anatomy.csv")
        for side in ("raw", "svc"):
            for scope in self.parameters["scopes"]:
                try:
                    coords = self._parent_coords(side, scope)
                except (KeyError, ValueError) as exc:
                    unavailable(f"parent_anatomy:{side}:{scope}", str(exc), self.missing)
                    continue
                side_length = self.parameters["parent_window_side_microns"][scope] / self.state["microns_per_coordinate"]
                assignments = assign_square_windows(coords, window_side_length=side_length, origin=self.state["origin"])
                composition = aggregate_parent_anatomy(assignments, point_anatomy[side].reindex(assignments.index))
                composition["side"], composition["scope"] = side, scope
                parent_anatomy[side, scope] = composition
                self._save_table(composition, f"tables/parent_anatomy_{side}_{_slug(scope)}.csv")
                if side == "raw":
                    baseline = self.state.get("partitions", {}).get(("raw", scope))
                    if baseline is not None:
                        baseline_assignments = assign_square_windows(
                            coordinates(self.sample.raw, self.sample.spatial_key).reindex(baseline.index),
                            window_side_length=side_length, origin=self.state["origin"])
                        baseline_context = aggregate_parent_anatomy(
                            baseline_assignments, point_anatomy["raw"].reindex(baseline.index))
                        baseline_context["side"], baseline_context["scope"] = "raw", scope
                        baseline_context["population"] = "raw_leiden_cohort"
                        parent_anatomy["raw_baseline", scope] = baseline_context
                        self._save_table(baseline_context, f"tables/parent_anatomy_raw_baseline_{_slug(scope)}.csv")
        self.state["point_anatomy"], self.state["parent_anatomy"] = point_anatomy, parent_anatomy

    def stage_molecular(self):
        from .impact_molecular import run_molecular
        run_molecular(self)

    def stage_membership(self):
        if not self.parameters["membership_same_units"]:
            return
        from revise_analysis.methods.partition import compare_membership
        from revise_analysis.methods.regions import assign_square_windows
        membership = {}
        if self.state.get("svc_labels") is None:
            unavailable("membership", "SVC reconstruction labels unavailable", self.missing)
            self.state["membership"] = membership
            return
        for scope in self.parameters["scopes"]:
            raw = self.state.get("partitions", {}).get(("raw", scope))
            if raw is None:
                unavailable(f"membership:{scope}", "Raw baseline partition unavailable", self.missing)
                continue
            # Raw defines the cohort; SVC broad labels do not silently redefine it.
            svc = self.state["svc_labels"]
            common = raw.index.intersection(svc.index)
            if common.empty:
                unavailable(f"membership:{scope}", "no common unit IDs", self.missing)
                continue
            comparison = compare_membership(raw.loc[common], svc.loc[common])
            coords = coordinates(self.sample.raw, self.sample.spatial_key).reindex(comparison.assignments.index)
            changed = coords.join(comparison.assignments)
            if "raw" in self.state.get("point_anatomy", {}):
                changed = changed.join(self.state["point_anatomy"]["raw"][["anatomy_region"]])
            state = self.state.get("region_tables", {}).get(("state", scope))
            if state is not None:
                side = self.parameters["parent_window_side_microns"][scope] / self.state["microns_per_coordinate"]
                assigned = assign_square_windows(coords, window_side_length=side, origin=self.state["origin"])
                changed["state_region"] = assigned["window_id"].map(state.set_index("window_id")["in_region"])
            membership[scope] = changed
            self._save_table(changed, f"tables/spatial_changed_map_{_slug(scope)}.csv")
            self._save_table(comparison.assignments, f"tables/membership_{_slug(scope)}.csv")
            self._save_table(comparison.contingency, f"tables/membership_contingency_{_slug(scope)}.csv")
            self._save_json({"summary": comparison.summary, "mapping": comparison.mapping}, f"tables/membership_summary_{_slug(scope)}.json")
        self.state["membership"] = membership

    def stage_integration(self):
        from revise_analysis.methods.regions import assign_square_windows
        if "origin" not in self.state or not self.state.get("point_anatomy"):
            unavailable("integrated_spatial_evidence", "spatial geometry or anatomy mapping unavailable", self.missing)
            return
        for scope in self.parameters["scopes"]:
            if ("state", scope) not in self.state.get("region_tables", {}):
                continue
            slug, svc = _slug(scope), self._svc_parent(scope)
            coords = coordinates(svc, self.sample.spatial_key)
            side = self.parameters["parent_window_side_microns"][scope] / self.state["microns_per_coordinate"]
            assigned = assign_square_windows(coords, window_side_length=side, origin=self.state["origin"])
            units = assigned[["window_id"]].join(self.state["point_anatomy"]["svc"][["anatomy_region"]])
            state = self.state["region_tables"]["state", scope].set_index("window_id")["in_region"]
            units["state_region"] = units["window_id"].map(state)
            gain = self.state.get("gains", {}).get(scope)
            units["gain_region"] = (units["window_id"].map(gain.set_index("window_id")["in_region"])
                                    if gain is not None else pd.NA)
            units["reconstruction_label"] = self.state["svc_labels"].reindex(units.index)
            units["scope"] = scope
            self._save_table(units, f"tables/integrated_label_units_{slug}.csv")
            composition = label_composition(units, self.state["svc_labels"], group_columns=["state_region", "anatomy_region"])
            composition["scope"] = scope
            composition["evidence_table"] = f"tables/integrated_label_units_{slug}.csv"
            self._save_table(composition, f"tables/integrated_label_composition_{slug}.csv")
            summary_rows = []
            for kind in ("state", "gain"):
                column = f"{kind}_region"
                for (region, anatomy_region), group in units.groupby([column, "anatomy_region"], dropna=False):
                    summary_rows.append({"scope": scope, "region_kind": kind, "in_region": region,
                                         "anatomy_region": anatomy_region, "n_units": int(group.shape[0]),
                                         "denominator_units": int(units.shape[0]),
                                         "fraction_of_parent_units": float(group.shape[0] / units.shape[0]),
                                         "region_definition_status": "available" if units[column].notna().any() else "unavailable",
                                         "evidence_table": f"tables/integrated_label_units_{slug}.csv"})
            self._save_table(pd.DataFrame(summary_rows), f"tables/integrated_region_anatomy_summary_{slug}.csv")
            anatomy = self.state["parent_anatomy"]["svc", scope]
            facts = self.state["region_tables"]["state", scope].merge(anatomy, on="window_id", how="left", suffixes=("", "_anatomy"))
            facts["evidence_type"], facts["scope"] = "state_x_anatomy", scope
            facts["region_definition_status"] = facts["region_status"]
            facts["evidence_table"] = f"tables/state_{slug}_region_windows.csv"
            if gain is not None:
                gain_facts = gain.merge(anatomy, on="window_id", how="left", suffixes=("", "_anatomy"))
                gain_facts["evidence_type"], gain_facts["scope"] = "delta_vs_raw_leiden_x_anatomy", scope
                gain_facts["region_definition_status"] = gain_facts["region_status"]
                gain_facts["evidence_table"] = f"tables/gain_{slug}_region_windows.csv"
                facts = pd.concat([facts, gain_facts], ignore_index=True, sort=False)
            self._save_table(facts, f"tables/integrated_spatial_evidence_{slug}.csv")
            changed = self.state.get("membership", {}).get(scope)
            if changed is not None:
                groups = [x for x in ("state_region", "anatomy_region") if x in changed]
                if groups:
                    summary = changed_unit_summary(changed, group_columns=groups)
                    summary["scope"] = scope
                    self._save_table(summary, f"tables/integrated_changed_units_{slug}.csv")
        from .impact_molecular import aggregate_programs
        aggregate_programs(self)

    def stage_figures(self):
        """Render declared current-run tables and attach explicit reading sections."""
        from revise_analysis.plotting import (
            plot_anatomy_context, plot_gain, plot_membership_change, plot_moran_distribution,
            plot_partition_sizes, plot_pathway_scores, plot_window_field,
            plot_support_curves, plot_threshold_bootstrap,
        )
        from revise_analysis.plotting.impact import plot_reconstruction_labels, plot_diversity_distribution
        tables = {Path(key).name: key for key in list(self.outputs)
                  if key.startswith("tables/") and key.endswith(".csv")}
        def read(name):
            return pd.read_csv(self.output_dir / tables[name])
        def draw(name, section, callback):
            relative = f"figures/{name}.png"
            callback(self.output_dir / relative)
            self.outputs[relative] = relative
            if relative not in self.stage_artifacts[section]:
                self.stage_artifacts[section].append(relative)
        if "partition_summary.csv" in tables:
            draw("partition_sizes", "baseline", lambda p: plot_partition_sizes(read("partition_summary.csv"), p))
        if "raw_anatomy_context.csv" in tables:
            draw("raw_anatomy_context", "anatomy", lambda p: plot_anatomy_context(read("raw_anatomy_context.csv"), p))
        if "reconstruction_labels_units.csv" in tables:
            labels = read("reconstruction_labels_units.csv")
            if {"x", "y"}.issubset(labels):
                draw("reconstruction_labels", "input", lambda p: plot_reconstruction_labels(labels, p))
        for name in tables:
            if not (name.startswith(("window_scale_support_", "gain_", "spatial_changed_map_", "moran_", "pathway_"))
                    or name.endswith(("_threshold_bootstrap.csv", "_region_windows.csv"))):
                continue
            frame = read(name)
            stem = Path(name).stem
            if name.startswith("window_scale_support_"):
                draw(stem, "support", lambda p: plot_support_curves(frame, p, title="候选尺度的有效窗口支持"))
                retention = frame.copy()
                retention["retained_unit_fraction"] = retention["retained_parent_unit_fraction"].fillna(retention.get("svc_common_retained_fraction"))
                draw(stem + "_retention", "support", lambda p: plot_support_curves(retention, p, y_column="retained_unit_fraction", title="单位保留比例（共同曲线使用 SVC 分母）"))
                coverage = frame.loc[frame.support_kind.str.startswith("common_")]
                draw(stem + "_common_coverage", "support", lambda p: plot_support_curves(coverage, p, y_column="common_fraction_of_svc_valid", title="共同支持 / SVC 有效窗口"))
            elif name.endswith("_threshold_bootstrap.csv"):
                draw(stem, "regions", lambda p: plot_threshold_bootstrap(frame, p, title=_figure_title(stem)))
            elif name.startswith("gain_") and name.endswith("_common_valid_windows.csv"):
                if "delta_neff_vs_raw_leiden" in frame:
                    frame["gain_neff"] = frame["delta_neff_vs_raw_leiden"]
                elif "gain_neff" not in frame:
                    continue
                draw(stem, "regions", lambda p: plot_gain(frame, p))
            elif name.endswith("_region_windows.csv"):
                value = "delta_neff_vs_raw_leiden" if "delta_neff_vs_raw_leiden" in frame else "neff"
                title = _figure_title(stem.removesuffix("_region_windows"))
                draw(stem, "regions", lambda p: plot_window_field(frame, p, value_column=value, title=title))
                draw(stem + "_distribution", "regions", lambda p: plot_diversity_distribution(frame, p, value_column=value, title=title + " 分布"))
                if "n_units" in frame:
                    support = frame.drop(columns=["in_region"], errors="ignore")
                    draw(stem + "_support", "support", lambda p: plot_window_field(support, p, value_column="n_units", title=title + " 单位支持"))
            elif name.startswith("spatial_changed_map_"):
                draw(stem, "membership", lambda p: plot_membership_change(frame, p, title=_figure_title(stem)))
            elif name.startswith("moran_") and "shared_genes" not in name:
                draw(stem, "molecular", lambda p: plot_moran_distribution(frame, p, title=_figure_title(stem)))
            elif name.startswith("pathway_") and "score" in frame:
                draw(stem, "molecular", lambda p: plot_pathway_scores(frame["score"], p, title=_figure_title(stem)))

        # Relationship views reuse the saved summaries and fixed unit scores.
        # They add presentation artifacts, never another scientific stage.
        from revise_analysis.plotting.relationships import render_relationship_figures
        for relative in render_relationship_figures(self.output_dir, self.outputs):
            self.outputs[relative] = relative
            section = "molecular" if Path(relative).stem.startswith("relationships_moran_shared_genes_") else "integration"
            if relative not in self.stage_artifacts[section]:
                self.stage_artifacts[section].append(relative)


def _figure_title(stem: str) -> str:
    """Give saved Impact figures a Chinese reader-facing title."""
    value = stem.replace("_", " ")
    replacements = (
        ("state ", "State "),
        ("gain ", "Gain "),
        ("moran ", "Moran "),
        ("pathway ", "通路 "),
        ("spatial changed map ", "空间 membership 变化图 "),
        ("threshold bootstrap", "阈值 bootstrap"),
        ("region windows", "区域窗口"),
        ("raw", "Raw"),
        ("svc", "SVC"),
    )
    for source, target in replacements:
        value = value.replace(source, target)
    return value.strip()


def raw_anatomy_context(sample: Any, *, tumor_label="Tumor", normal_label="Intestinal Epithelial"):
    """Raw source labels; authoritative Interface is assigned at window level."""
    frame = coordinates(sample.raw, sample.spatial_key)
    frame["broad_label"] = sample.labels("raw", sample.broad_key).reindex(frame.index)
    frame["anatomy_source_class"] = np.select(
        [frame["broad_label"] == tumor_label, frame["broad_label"] == normal_label],
        ["Tumor", "Normal"], default="Other")
    return frame


def native_window_diversity(sample, side, labels, *, adata=None, window_side_length, origin,
                            min_units, n_draws, random_state):
    from revise_analysis.methods.regions import compute_window_diversity
    adata = getattr(sample, side) if adata is None else adata
    return compute_window_diversity(coordinates(adata, sample.spatial_key).reindex(labels.index), labels,
        window_side_length=window_side_length, origin=origin, min_units=min_units,
        n_draws=n_draws, random_state=random_state)


def compute_partitions(adata, *, resolution, n_top_genes, random_state,
                       transformation_state="untransformed_nonnegative"):
    from revise_analysis.methods.partition import prepare_partition_graph, partition_prepared_graph
    prepared = prepare_partition_graph(adata, n_top_genes=n_top_genes, random_state=random_state,
        min_genes=0, min_cells=0, transformation_state=transformation_state)
    result = partition_prepared_graph(prepared, resolution=resolution)
    # Private workflow handle: every K-control candidate uses this exact graph.
    result["_prepared"] = prepared
    return result


def partition_scope(adata, labels, scope, *, resolution, n_top_genes, random_state,
                    transformation_state="untransformed_nonnegative"):
    if scope != "All":
        ids = labels.index[labels == scope]
        if ids.empty:
            raise ValueError(f"no units labelled {scope!r}")
        adata = adata[ids]
    return compute_partitions(adata, resolution=resolution, n_top_genes=n_top_genes,
        random_state=random_state, transformation_state=transformation_state)


def prepare_side(sample, side, *, sample_n_units, random_state):
    adata = deterministic_subset(getattr(sample, side), sample_n_units, random_state)
    return adata, sample.labels(side).reindex(adata.obs_names)


def compute_impact(sample, output_dir, parameters=None):
    return run(sample, output_dir, parameters)


def run(sample, output_dir, parameters=None, *, continue_on_error=False):
    return ImpactWorkflow(sample, output_dir, parameters, continue_on_error=continue_on_error).run()
