"""Relationship figures rendered only from explicitly registered saved tables.

The renderers in this module are deliberately descriptive.  They do not join
back to H5AD inputs, derive regions, rescore programs, or combine Raw and SVC
into a shared scientific object.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .impact import _pyplot, _save


Renderer = Callable[[pd.DataFrame, Path, dict[str, Path]], bool]


def render_relationship_figures(output_dir: Path, outputs: dict) -> dict[str, str]:
    """Render relationship figures from registered, contained CSV artifacts.

    ``outputs`` is the workflow artifact registry.  Only its relative-path
    values are considered; the output directory is never scanned.  Missing,
    escaped, unsupported, or scientifically empty tables are ignored.  The
    returned mapping follows the workflow registry convention.
    """
    root = Path(output_dir).resolve()
    registered = _registered_files(root, outputs)
    metadata = {
        relative.name: path
        for relative, path in registered
        if relative.suffix.lower() == ".json"
    }
    rendered: dict[str, str] = {}
    for relative, source in registered:
        if relative.suffix.lower() != ".csv":
            continue
        renderer = _renderer_for(relative.name)
        if renderer is None:
            continue
        try:
            table = pd.read_csv(source)
        except (OSError, pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeError):
            continue
        destination_relative = Path("figures") / f"relationships_{relative.stem}.png"
        destination = root / destination_relative
        if renderer(table, destination, metadata):
            value = destination_relative.as_posix()
            rendered[value] = value
    return rendered


def _registered_files(root: Path, outputs: dict) -> list[tuple[Path, Path]]:
    files: list[tuple[Path, Path]] = []
    seen: set[Path] = set()
    for value in outputs.values():
        if not isinstance(value, (str, Path)):
            continue
        relative = Path(value)
        if relative.is_absolute() or relative in seen:
            continue
        candidate = (root / relative).resolve()
        try:
            candidate.relative_to(root)
        except ValueError:
            continue
        if not candidate.is_file():
            continue
        seen.add(relative)
        files.append((relative, candidate))
    return sorted(files, key=lambda item: item[0].as_posix())


def _renderer_for(name: str) -> Renderer | None:
    if name.startswith("integrated_region_anatomy_summary_"):
        return _render_region_anatomy
    if name.startswith("integrated_label_composition_"):
        return _render_label_composition
    if name.startswith("integrated_program_") and name.endswith("_units.csv"):
        return _render_program_units
    if name.startswith("integrated_program_") and name.endswith(("_anatomy.csv", "_region.csv")):
        return _render_program_summary
    if name.startswith("integrated_changed_units_"):
        return _render_changed_units
    if name.startswith("moran_shared_genes_"):
        return _render_moran_shared
    return None


def _render_region_anatomy(table: pd.DataFrame, destination: Path, metadata: dict[str, Path]) -> bool:
    required = {
        "region_kind", "in_region", "anatomy_region", "n_units",
        "denominator_units", "fraction_of_parent_units",
    }
    if not required.issubset(table):
        return False
    frame = table.copy()
    frame["_fraction"] = pd.to_numeric(frame["fraction_of_parent_units"], errors="coerce")
    frame["_n_units"] = pd.to_numeric(frame["n_units"], errors="coerce")
    frame["_denominator"] = pd.to_numeric(frame["denominator_units"], errors="coerce")
    frame["_region"] = frame["in_region"].map(_region_bucket)
    frame = frame.loc[
        frame["_fraction"].notna() & frame["_n_units"].notna()
        & frame["_denominator"].gt(0)
    ].copy()
    panels: list[tuple[str, pd.DataFrame]] = []
    for kind, subset in frame.groupby("region_kind", sort=True, dropna=False):
        if subset.empty or _all_region_definitions_unavailable(subset):
            continue
        panels.append((str(kind), subset))
    if not panels:
        return False

    plt = _pyplot()
    figure, axes = plt.subplots(1, len(panels), figsize=(6.2 * len(panels), 4.8), squeeze=False)
    colors = {"inside": "#C44E52", "outside": "#4C78A8", "unknown": "#B8BEC7"}
    region_labels = {
        "inside": "区域内（有效窗口）",
        "outside": "区域外（有效窗口）",
        "unknown": "区域状态未知",
    }
    for axis, (kind, subset) in zip(axes[0], panels, strict=True):
        anatomy = _ordered_strings(subset["anatomy_region"])
        bottoms = np.zeros(len(anatomy), dtype=float)
        for region in ("inside", "outside", "unknown"):
            values = _saved_group_values(subset, "anatomy_region", "_region", region, "_fraction", anatomy)
            axis.bar(anatomy, values, bottom=bottoms, color=colors[region], label=region_labels[region])
            bottoms += values
        denominator = _denominator_label(subset["_denominator"])
        axis.set(
            title=f"{kind} × Anatomy（分母 {denominator}）",
            xlabel="Anatomy", ylabel="占 parent units 的保存比例",
        )
        axis.tick_params(axis="x", rotation=35)
        axis.set_ylim(0, max(1.0, float(bottoms.max(initial=0))))
    axes[0, 0].legend(title="区域状态")
    figure.suptitle("区域与 Anatomy 组成（各 region kind 独立分母）")
    _save(figure, destination)
    return True


def _render_label_composition(table: pd.DataFrame, destination: Path, metadata: dict[str, Path]) -> bool:
    required = {
        "state_region", "anatomy_region", "reconstruction_label",
        "n_units", "fraction", "denominator_units",
    }
    if not required.issubset(table):
        return False
    frame = table.copy()
    frame["_fraction"] = pd.to_numeric(frame["fraction"], errors="coerce")
    frame["_n_units"] = pd.to_numeric(frame["n_units"], errors="coerce")
    frame["_denominator"] = pd.to_numeric(frame["denominator_units"], errors="coerce")
    frame["_state"] = frame["state_region"].map(_region_bucket)
    frame["_anatomy"] = frame["anatomy_region"].astype("string").fillna("Unknown").astype(str)
    frame["_label"] = frame["reconstruction_label"].astype("string").fillna("Unknown").astype(str)
    frame = frame.loc[
        frame["_fraction"].ge(0) & frame["_n_units"].ge(0) & frame["_denominator"].gt(0)
    ].copy()
    if frame.empty:
        return False

    states = [state for state in ("inside", "outside", "unknown") if (frame["_state"] == state).any()]
    labels = _ordered_strings(frame["_label"])
    groups = [
        (state, anatomy)
        for state in states
        for anatomy in _ordered_strings(frame.loc[frame["_state"] == state, "_anatomy"])
    ]
    plt = _pyplot()
    figure, axis = plt.subplots(figsize=(8.8, max(4.8, .48 * len(groups) + 1.8)))
    palette = plt.get_cmap("tab20")
    y = np.arange(len(groups))
    left = np.zeros(len(groups), dtype=float)
    for number, label in enumerate(labels):
        fractions, counts = [], []
        for state, anatomy in groups:
            subset = frame.loc[
                (frame["_state"] == state) & (frame["_anatomy"] == anatomy) & (frame["_label"] == label)
            ]
            fractions.append(float(subset["_fraction"].sum()) if not subset.empty else 0.0)
            counts.append(float(subset["_n_units"].sum()) if not subset.empty else 0.0)
        values = np.asarray(fractions)
        axis.barh(y, values, left=left, color=palette(number % 20), label=label)
        for row, (start, width, count) in enumerate(zip(left, values, counts, strict=True)):
            if width >= .08:
                axis.text(start + width / 2, row, f"{int(count)}", ha="center", va="center", fontsize=7)
        left += values
    state_labels = {"inside": "区域内（有效窗口）", "outside": "区域外（有效窗口）", "unknown": "区域状态未知"}
    row_labels = []
    for state, anatomy in groups:
        denominator = _denominator_label(frame.loc[
            (frame["_state"] == state) & (frame["_anatomy"] == anatomy), "_denominator"
        ])
        row_labels.append(f"State {state_labels[state]} · {anatomy} (n={denominator})")
    axis.set_yticks(y, row_labels)
    axis.invert_yaxis()
    axis.set(
        title="State × Anatomy 的重建标签组成（柱内为单位数）",
        xlabel="保存的 label fraction", ylabel="State 状态 · Anatomy（分母）",
    )
    axis.set_xlim(0, max(1.0, float(left.max(initial=0))))
    axis.legend(title="重建标签", bbox_to_anchor=(1.02, 1), loc="upper left")
    _save(figure, destination)
    return True


def _render_program_units(table: pd.DataFrame, destination: Path, metadata: dict[str, Path]) -> bool:
    if not {"x", "y", "score"}.issubset(table):
        return False
    frame = table.copy()
    frame["_x"] = pd.to_numeric(frame["x"], errors="coerce")
    frame["_y"] = pd.to_numeric(frame["y"], errors="coerce")
    frame["_score"] = pd.to_numeric(frame["score"], errors="coerce")
    coordinates = frame["_x"].notna() & frame["_y"].notna()
    scored = coordinates & frame["_score"].notna()
    if not scored.any():
        return False
    sides = frame["side"].dropna().astype(str).unique() if "side" in frame else np.array([])
    if len(sides) > 1:
        return False

    figure, axis = _pyplot().subplots(figsize=(6.4, 5.2))
    missing = coordinates & frame["_score"].isna()
    if missing.any():
        axis.scatter(
            frame.loc[missing, "_x"], frame.loc[missing, "_y"],
            s=8, color="#CDD2D9", marker="x", alpha=.65,
            label="未评分（保持缺失）", rasterized=True,
        )
    points = axis.scatter(
        frame.loc[scored, "_x"], frame.loc[scored, "_y"],
        c=frame.loc[scored, "_score"], s=9, cmap="viridis", alpha=.85, rasterized=True,
    )
    figure.colorbar(points, ax=axis, label="保存的 program score")
    side = sides[0] if len(sides) == 1 else ""
    program = _single_value(frame, "program")
    scope = _single_value(frame, "scope")
    details = [
        f"n_scored={int(frame['_score'].notna().sum())}/{len(frame)}",
        f"plot_count={int(scored.sum())}",
    ]
    coverage = _program_coverage(metadata, side, program)
    if coverage is not None:
        details.append(f"gene coverage={coverage}")
    axis.set(
        title=" | ".join(value for value in (side, scope, program) if value) + "\n" + ", ".join(details),
        xlabel="x", ylabel="y", aspect="equal",
    )
    if missing.any():
        axis.legend(loc="best")
    _save(figure, destination)
    return True


def _render_program_summary(table: pd.DataFrame, destination: Path, metadata: dict[str, Path]) -> bool:
    required = {"n_units", "n_scored", "score_mean", "score_median"}
    if not required.issubset(table):
        return False
    frame = table.copy()
    for source, target in (
        ("n_units", "_n_units"), ("n_scored", "_n_scored"),
        ("score_mean", "_mean"), ("score_median", "_median"),
    ):
        frame[target] = pd.to_numeric(frame[source], errors="coerce")
    frame = frame.loc[frame["_n_units"].gt(0) & frame["_n_scored"].ge(0)].copy()
    valid_summary = frame["_n_scored"].gt(0) & (frame["_mean"].notna() | frame["_median"].notna())
    if frame.empty or not valid_summary.any():
        return False
    if "side" in frame and frame["side"].dropna().astype(str).nunique() > 1:
        return False

    frame["_group"] = frame.apply(_program_group_label, axis=1)
    figure, axis = _pyplot().subplots(figsize=(7.2, max(3.8, .52 * len(frame) + 1.5)))
    y = np.arange(len(frame))
    for row, (_, values) in enumerate(frame.iterrows()):
        if pd.notna(values["_mean"]) and pd.notna(values["_median"]):
            axis.plot([values["_mean"], values["_median"]], [row, row], color="#AAB2BC", linewidth=1.2)
        if values["_n_scored"] == 0:
            axis.text(
                .01, row, f"未评分 0/{int(values['_n_units'])}",
                transform=axis.get_yaxis_transform(), va="center", color="#6B7280", fontsize=8,
            )
        elif pd.isna(values["_mean"]) and pd.isna(values["_median"]):
            axis.text(
                .01, row, "评分汇总缺失", transform=axis.get_yaxis_transform(),
                va="center", color="#6B7280", fontsize=8,
            )
    axis.scatter(
        frame.loc[valid_summary, "_mean"], y[valid_summary.to_numpy()],
        marker="o", color="#4C78A8", label="mean",
    )
    axis.scatter(
        frame.loc[valid_summary, "_median"], y[valid_summary.to_numpy()],
        marker="s", facecolors="white", edgecolors="#C44E52", label="median",
    )
    labels = [
        f"{group} (n={int(scored)}/{int(total)})"
        for group, scored, total in zip(frame["_group"], frame["_n_scored"], frame["_n_units"], strict=True)
    ]
    axis.set_yticks(y, labels)
    axis.set(
        title=" | ".join(value for value in (
            _single_value(frame, "side"), _single_value(frame, "scope"), _single_value(frame, "program"),
        ) if value),
        xlabel="保存的 program score 汇总", ylabel="保存的分组（评分数/单位数）",
    )
    axis.legend()
    _save(figure, destination)
    return True


def _render_changed_units(table: pd.DataFrame, destination: Path, metadata: dict[str, Path]) -> bool:
    required = {"n_compared", "n_changed", "changed_fraction"}
    if not required.issubset(table):
        return False
    frame = table.copy()
    frame["_compared"] = pd.to_numeric(frame["n_compared"], errors="coerce")
    frame["_changed"] = pd.to_numeric(frame["n_changed"], errors="coerce")
    frame["_fraction"] = pd.to_numeric(frame["changed_fraction"], errors="coerce")
    frame = frame.loc[
        frame["_compared"].gt(0) & frame["_changed"].ge(0) & frame["_fraction"].notna()
    ].copy()
    if frame.empty:
        return False
    frame["_group"] = frame.apply(_changed_group_label, axis=1)
    figure, axis = _pyplot().subplots(figsize=(7, max(3.6, .5 * len(frame) + 1.4)))
    y = np.arange(len(frame))
    bars = axis.barh(y, frame["_fraction"], color="#B55A69")
    axis.set_yticks(y, frame["_group"])
    for bar, changed, compared in zip(bars, frame["_changed"], frame["_compared"], strict=True):
        axis.text(
            bar.get_width(), bar.get_y() + bar.get_height() / 2,
            f" {int(changed)}/{int(compared)}", va="center", fontsize=8,
        )
    axis.set(
        title=f"Membership 变化比例 | {_single_value(frame, 'scope')}",
        xlabel="changed fraction（标注为 changed/compared）", ylabel="保存的分组",
        xlim=(0, max(1.0, float(frame["_fraction"].max()))),
    )
    _save(figure, destination)
    return True


def _render_moran_shared(table: pd.DataFrame, destination: Path, metadata: dict[str, Path]) -> bool:
    required = {"raw_moran", "svc_moran", "raw_n_units", "svc_n_units", "n_common_genes"}
    if not required.issubset(table):
        return False
    frame = table.copy()
    frame["_raw"] = pd.to_numeric(frame["raw_moran"], errors="coerce")
    frame["_svc"] = pd.to_numeric(frame["svc_moran"], errors="coerce")
    available = frame["_raw"].notna() & frame["_svc"].notna()
    if "comparison_available" in frame:
        available &= frame["comparison_available"].map(_saved_true)
    frame = frame.loc[available].copy()
    if frame.empty:
        return False

    figure, axis = _pyplot().subplots(figsize=(6, 5.4))
    axis.scatter(frame["_raw"], frame["_svc"], s=18, alpha=.65, color="#4C78A8", rasterized=True)
    low = float(min(frame["_raw"].min(), frame["_svc"].min()))
    high = float(max(frame["_raw"].max(), frame["_svc"].max()))
    if low == high:
        low, high = low - .5, high + .5
    axis.plot([low, high], [low, high], color="#999", linewidth=1, linestyle="--")
    common = f"common genes={_integer_values(table['n_common_genes'])}; plotted genes={len(frame)}"
    raw = f"Raw units={_integer_values(table['raw_n_units'])}"
    svc = f"SVC units={_integer_values(table['svc_n_units'])}"
    if "raw_n_native_genes" in table:
        raw += f"; native genes={_integer_values(table['raw_n_native_genes'])}"
    if "svc_n_native_genes" in table:
        svc += f"; native genes={_integer_values(table['svc_n_native_genes'])}"
    axis.set(
        title=(f"Raw vs SVC 共享基因 Moran I | {_single_value(table, 'scope')}\n"
               f"{raw}\n{svc}\n{common}"),
        xlabel="Raw Moran I", ylabel="SVC Moran I", xlim=(low, high), ylim=(low, high), aspect="equal",
    )
    _save(figure, destination)
    return True


def _all_region_definitions_unavailable(frame: pd.DataFrame) -> bool:
    if "region_definition_status" not in frame:
        return (frame["_region"] == "unknown").all()
    statuses = frame["region_definition_status"].astype("string").str.lower().dropna()
    return not statuses.empty and statuses.eq("unavailable").all()


def _saved_group_values(
    frame: pd.DataFrame,
    row_column: str,
    series_column: str,
    series_value: str,
    value_column: str,
    row_order: list[str],
) -> np.ndarray:
    selected = frame.loc[
        frame[series_column].astype(str) == str(series_value), [row_column, value_column]
    ].copy()
    selected[row_column] = selected[row_column].astype(str)
    saved = selected.groupby(row_column, dropna=False)[value_column].sum(min_count=1)
    return saved.reindex(row_order, fill_value=0).to_numpy(dtype=float)


def _ordered_strings(values: pd.Series) -> list[str]:
    return sorted(values.astype("string").fillna("Unknown").astype(str).unique().tolist())


def _region_bucket(value) -> str:
    if pd.isna(value):
        return "unknown"
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "inside"}:
        return "inside"
    if normalized in {"false", "0", "outside"}:
        return "outside"
    return "unknown"


def _saved_true(value) -> bool:
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"true", "1"}


def _denominator_label(values: pd.Series) -> str:
    return _integer_values(values)


def _integer_values(values: pd.Series) -> str:
    numeric = pd.to_numeric(values, errors="coerce").dropna().astype(int).unique().tolist()
    return "/".join(str(value) for value in sorted(numeric)) or "unknown"


def _single_value(frame: pd.DataFrame, column: str) -> str:
    if column not in frame:
        return ""
    values = frame[column].dropna().astype(str).unique()
    return values[0] if len(values) == 1 else ""


def _program_group_label(row: pd.Series) -> str:
    parts: list[str] = []
    if "state_region" in row.index:
        parts.append(f"state={_region_bucket(row['state_region'])}")
    if "gain_region" in row.index:
        parts.append(f"gain={_region_bucket(row['gain_region'])}")
    if "anatomy_region" in row.index:
        parts.append(f"anatomy={row['anatomy_region'] if pd.notna(row['anatomy_region']) else 'Unknown'}")
    if "window_id" in row.index:
        parts.append(f"window={row['window_id']}")
    return ", ".join(parts) or "saved group"


def _changed_group_label(row: pd.Series) -> str:
    parts: list[str] = []
    if "state_region" in row.index:
        parts.append(f"state={_region_bucket(row['state_region'])}")
    if "anatomy_region" in row.index:
        parts.append(f"anatomy={row['anatomy_region'] if pd.notna(row['anatomy_region']) else 'Unknown'}")
    return ", ".join(parts) or "all compared units"


def _program_coverage(metadata: dict[str, Path], side: str, program: str) -> str | None:
    if not side or not program:
        return None
    name = f"pathway_{side}_{_slug(program)}_metadata.json"
    path = metadata.get(name)
    if path is None:
        return None
    try:
        saved = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    coverage = saved.get("coverage", {})
    if not isinstance(coverage, dict):
        return None
    value = coverage.get("fraction_present")
    if not isinstance(value, (int, float)) or not np.isfinite(value):
        return None
    return str(value)


def _slug(value: str) -> str:
    return "".join(character if character.isalnum() else "_" for character in str(value)).strip("_") or "scope"
