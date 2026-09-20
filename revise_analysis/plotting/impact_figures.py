"""Section-scoped reconstruction-impact figures from saved artifacts only."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .impact import (
    plot_anatomy_context,
    plot_anatomy_neff_state,
    plot_baseline_window_metrics,
    plot_diversity_distribution,
    plot_gain,
    plot_membership_change,
    plot_partition_sizes,
    plot_reconstruction_labels,
    plot_support_curves,
    plot_window_support_grid,
    plot_threshold_bootstrap,
    plot_window_field,
)
from .programs import plot_pathway_scores
from .relationships import _registered_files, render_relationship_figures
from .spatial import plot_moran_distribution


SECTIONS = {
    "input", "baseline", "support", "diversity", "regions", "anatomy",
    "molecular", "membership", "integration",
}


def render_impact_figures(
    output_dir: Path,
    outputs: dict,
    parameters: dict,
    section: str | None = None,
) -> dict[str, str]:
    """Render registered saved tables for one reading section or the full batch.

    The returned mapping is ``relative figure path -> owning section``.  This
    function never scans the result directory, opens H5AD, or changes a saved
    scientific table.
    """
    if section is not None and section not in SECTIONS:
        raise ValueError(f"Unknown impact figure section {section!r}")
    root = Path(output_dir).resolve()
    registered = _registered_files(root, outputs)
    paths = {relative.as_posix(): path for relative, path in registered}
    tables = {
        relative.name: (relative.as_posix(), path)
        for relative, path in registered
        if relative.suffix.lower() == ".csv" and relative.parent.as_posix() == "tables"
    }
    rendered: dict[str, str] = {}

    def wanted(owner: str) -> bool:
        return section is None or section == owner

    def read(name: str) -> pd.DataFrame | None:
        record = tables.get(name)
        if record is None:
            return None
        return pd.read_csv(record[1])

    def draw(name: str, owner: str, callback: Callable[[Path], object]) -> None:
        if not wanted(owner):
            return
        relative = f"figures/{name}.png"
        callback(root / relative)
        rendered[relative] = owner

    if wanted("input"):
        frame = read("reconstruction_labels_units.csv")
        if frame is not None and {"x", "y"}.issubset(frame):
            draw("reconstruction_labels", "input", lambda path: plot_reconstruction_labels(frame, path))

    if wanted("baseline"):
        frame = read("partition_summary.csv")
        if frame is not None:
            draw("partition_sizes", "baseline", lambda path: plot_partition_sizes(frame, path))

    if wanted("support"):
        for name in sorted(tables):
            if not name.startswith("window_support_grid_"):
                continue
            frame = read(name)
            if frame is None:
                continue
            stem = Path(name).stem
            draw(stem, "support", lambda path, frame=frame, stem=stem: plot_window_support_grid(
                frame, path, title=_figure_title(stem),
            ))
        for name in sorted(tables):
            if not name.startswith("window_scale_support_"):
                continue
            frame = read(name)
            if frame is None:
                continue
            stem = Path(name).stem
            draw(stem, "support", lambda path, frame=frame: plot_support_curves(
                frame, path, title="候选尺度的有效窗口支持",
            ))
            retention = frame.copy()
            if "retained_parent_unit_fraction" in retention:
                retention["retained_unit_fraction"] = retention["retained_parent_unit_fraction"]
                if "svc_common_retained_fraction" in retention:
                    retention["retained_unit_fraction"] = retention["retained_unit_fraction"].fillna(
                        retention["svc_common_retained_fraction"]
                    )
                draw(stem + "_retention", "support", lambda path, frame=retention: plot_support_curves(
                    frame, path, y_column="retained_unit_fraction",
                    title="单位保留比例（共同曲线使用 SVC 分母）",
                ))
            if "support_kind" in frame and "common_fraction_of_svc_valid" in frame:
                coverage = frame.loc[frame.support_kind.astype(str).str.startswith("common_")].copy()
                if not coverage.empty:
                    draw(stem + "_common_coverage", "support", lambda path, frame=coverage: plot_support_curves(
                        frame, path, y_column="common_fraction_of_svc_valid",
                        title="共同支持 / SVC 有效窗口",
                    ))
        for name in sorted(tables):
            if not name.endswith("_region_windows.csv"):
                continue
            frame = read(name)
            if frame is None or "n_units" not in frame:
                continue
            support = frame.drop(columns=["in_region"], errors="ignore")
            stem = Path(name).stem
            draw(stem + "_support", "support", lambda path, frame=support, stem=stem: plot_window_field(
                frame, path, value_column="n_units", title=_figure_title(stem.removesuffix("_region_windows")) + " 单位支持",
            ))

    if wanted("diversity"):
        diversity: dict[str, dict[str, pd.DataFrame]] = {}
        for name in sorted(tables):
            parsed = _diversity_table_role(name)
            if parsed is None:
                continue
            role, scope = parsed
            frame = read(name)
            if frame is None:
                continue
            diversity.setdefault(scope, {})[role] = frame
            if {"window_x", "window_y", "neff"}.issubset(frame):
                stem = Path(name).stem
                draw(stem, "diversity", lambda path, frame=frame, stem=stem: plot_window_field(
                    frame, path, value_column="neff", title=_figure_title(stem),
                ))
        for scope, frames in diversity.items():
            if "svc_reconstruction" not in frames:
                continue
            draw(f"baseline_window_metrics_{scope}", "diversity", lambda path, frames=frames, scope=scope:
                 plot_baseline_window_metrics(frames, path, title=f"{scope}：四个局部多样性指标与多 baseline"))

    if wanted("regions"):
        for name in sorted(tables):
            frame = read(name)
            if frame is None:
                continue
            stem = Path(name).stem
            if name.endswith("_threshold_bootstrap.csv"):
                draw(stem, "regions", lambda path, frame=frame, stem=stem: plot_threshold_bootstrap(
                    frame, path, title=_figure_title(stem),
                ))
            elif name.startswith("gain_") and name.endswith("_common_valid_windows.csv"):
                gain = frame.copy()
                if "delta_neff_vs_raw_leiden" in gain:
                    gain["gain_neff"] = gain["delta_neff_vs_raw_leiden"]
                if "gain_neff" in gain:
                    draw(stem, "regions", lambda path, frame=gain: plot_gain(frame, path))
            elif name.startswith(("state_", "gain_")) and name.endswith("_region_windows.csv"):
                value = "delta_neff_vs_raw_leiden" if "delta_neff_vs_raw_leiden" in frame else "neff"
                if value not in frame:
                    continue
                title = _figure_title(stem.removesuffix("_region_windows"))
                draw(stem, "regions", lambda path, frame=frame, value=value, title=title: plot_window_field(
                    frame, path, value_column=value, title=title,
                ))
                draw(stem + "_distribution", "regions", lambda path, frame=frame, value=value, title=title:
                     plot_diversity_distribution(frame, path, value_column=value, title=title + " 分布"))

    if wanted("anatomy"):
        anatomy = read("svc_anatomy_context.csv")
        if anatomy is not None:
            draw("svc_anatomy_context", "anatomy", lambda path: plot_anatomy_context(anatomy, path))
        windows = read("svc_anatomy_windows.csv")
        if windows is not None:
            micron_scale = _microns_per_coordinate(paths.get("tables/input_manifest.json"))
            anatomy_microns = _finite_positive(parameters.get("anatomy_window_side_microns"))
            anatomy_side = anatomy_microns / micron_scale if anatomy_microns is not None and micron_scale is not None else None
            configured = parameters.get("parent_window_side_microns", {})
            for name in sorted(tables):
                if not name.startswith("state_") or not name.endswith("_region_windows.csv"):
                    continue
                scope = name[len("state_"):-len("_region_windows.csv")]
                state_frame = read(name)
                if state_frame is None:
                    continue
                parent_microns = _scope_parameter(configured, scope)
                parent_microns = _finite_positive(parent_microns)
                parent_side = parent_microns / micron_scale if parent_microns is not None and micron_scale is not None else None
                draw(f"anatomy_neff_state_{scope}", "anatomy", lambda path, state_frame=state_frame,
                     scope=scope, anatomy_side=anatomy_side, parent_side=parent_side: plot_anatomy_neff_state(
                         windows, state_frame, path, anatomy_window_side=anatomy_side,
                         parent_window_side=parent_side, title=f"{scope}：SVC Anatomy、SVC Neff 与 State",
                     ))

    if wanted("molecular"):
        for name in sorted(tables):
            frame = read(name)
            if frame is None:
                continue
            stem = Path(name).stem
            if name.startswith("moran_") and "shared_genes" not in name:
                draw(stem, "molecular", lambda path, frame=frame, stem=stem: plot_moran_distribution(
                    frame, path, title=_figure_title(stem),
                ))
            elif name.startswith("pathway_") and "score" in frame:
                draw(stem, "molecular", lambda path, frame=frame, stem=stem: plot_pathway_scores(
                    frame["score"], path, title=_figure_title(stem),
                ))

    if wanted("membership"):
        for name in sorted(tables):
            if not name.startswith("spatial_changed_map_"):
                continue
            frame = read(name)
            if frame is None:
                continue
            stem = Path(name).stem
            draw(stem, "membership", lambda path, frame=frame, stem=stem: plot_membership_change(
                frame, path, title=_figure_title(stem),
            ))

    relationship_sections = ("molecular", "integration") if section is None else (section,)
    for owner in relationship_sections:
        if owner not in {"molecular", "integration"}:
            continue
        for relative in render_relationship_figures(root, outputs, section=owner):
            rendered[relative] = owner
    return rendered


def _diversity_table_role(name: str) -> tuple[str, str] | None:
    patterns = (
        ("window_diversity_svc_leiden_baseline_", "svc_leiden"),
        ("window_diversity_raw_level2_baseline_", "raw_level2"),
        ("window_diversity_raw_k_control_", "raw_k_control"),
        ("window_diversity_svc_", "svc_reconstruction"),
        ("window_diversity_raw_", "raw_leiden"),
    )
    for prefix, role in patterns:
        if name.startswith(prefix) and name.endswith(".csv"):
            return role, name[len(prefix):-4]
    return None


def _microns_per_coordinate(manifest_path: Path | None) -> float | None:
    if manifest_path is None:
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    candidates = [
        manifest.get("spatial", {}).get("microns_per_coordinate") if isinstance(manifest.get("spatial"), dict) else None,
        manifest.get("config", {}).get("spatial", {}).get("microns_per_coordinate")
        if isinstance(manifest.get("config"), dict) and isinstance(manifest.get("config", {}).get("spatial"), dict) else None,
    ]
    for value in candidates:
        parsed = _finite_positive(value)
        if parsed is not None:
            return parsed
    return None


def _finite_positive(value) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if np.isfinite(parsed) and parsed > 0 else None


def _scope_parameter(configured, artifact_scope: str):
    if not isinstance(configured, dict):
        return configured
    if artifact_scope in configured:
        return configured[artifact_scope]
    return next(
        (value for scope, value in configured.items() if _slug(scope) == artifact_scope),
        None,
    )


def _slug(value: str) -> str:
    return "".join(character if character.isalnum() else "_" for character in str(value)).strip("_") or "scope"


def _figure_title(stem: str) -> str:
    value = stem.replace("_", " ")
    for source, target in (
        ("state ", "State "), ("gain ", "Gain "), ("moran ", "Moran "),
        ("pathway ", "通路 "), ("spatial changed map ", "空间 membership 变化图 "),
        ("threshold bootstrap", "阈值 bootstrap"), ("region windows", "区域窗口"),
        ("window diversity", "窗口多样性"), ("raw", "Raw"), ("svc", "SVC"),
    ):
        value = value.replace(source, target)
    return value.strip()
