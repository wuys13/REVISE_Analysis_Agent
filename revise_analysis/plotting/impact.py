"""Reusable impact figures from tabular workflow outputs."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd


def plot_partition_sizes(summary: pd.DataFrame, destination: Path) -> Path:
    """Plot saved partition summary counts by side and scope."""
    plt = _pyplot()
    if not {"side", "scope"}.issubset(summary):
        raise ValueError("partition summary requires side and scope columns")
    count_column = next((column for column in ("n_clusters", "n_units") if column in summary), None)
    if count_column is None:
        raise ValueError("partition summary needs n_units or n_clusters")
    pivot = summary.pivot(index="scope", columns="side", values=count_column)
    axis = pivot.plot.bar(rot=0)
    axis.set_ylabel({"n_clusters": "分群数", "n_units": "单位数"}.get(count_column, count_column.replace("_", " ")))
    axis.set_title("两侧独立分群")
    return _save(axis.figure, destination)


def plot_gain(gain: pd.DataFrame, destination: Path) -> Path:
    """Plot only valid windows already selected by the workflow."""
    plt = _pyplot()
    if "gain_neff" not in gain:
        raise ValueError("gain table requires gain_neff")
    values = pd.to_numeric(gain["gain_neff"], errors="coerce").dropna()
    figure, axis = _pyplot().subplots()
    if values.empty:
        _no_values(axis, "没有可计算的共同有效窗口 Gain")
    else:
        axis.hist(values, bins=min(30, max(5, len(values))))
    axis.set_xlabel("SVC − Raw 有效多样性")
    axis.set_title("共同有效空间窗口中的 ΔNeff")
    return _save(axis.figure, destination)


def plot_anatomy_context(context: pd.DataFrame, destination: Path) -> Path:
    """Draw the saved SVC-defined Anatomy context in native coordinates."""
    category = "level1_region" if "level1_region" in context else "anatomy_class"
    required = {"x", "y", category}
    if missing := required - set(context.columns):
        raise ValueError(f"anatomy context missing columns: {sorted(missing)}")
    plt = _pyplot()
    colors = {"Tumor": "#c43c39", "Normal": "#39835b", "Interface": "#e89b20", "Other": "#87929e"}
    figure, axis = plt.subplots(figsize=(6, 5))
    for label, frame in context.groupby(category, dropna=False, sort=True):
        axis.scatter(frame["x"], frame["y"], s=7, alpha=.75,
                     color=colors.get(str(label), "#87929e"), label=str(label), rasterized=True)
    axis.set(title="SVC-defined Anatomy 窗口背景", xlabel="x", ylabel="y", aspect="equal")
    axis.legend(title="窗口背景", markerscale=2)
    return _save(figure, destination)


def plot_membership_change(assignments: pd.DataFrame, destination: Path, *, title: str) -> Path:
    """Map the explicit Raw-defined cohort membership comparison only."""
    required = {"x", "y", "unit_changed"}
    if missing := required - set(assignments.columns):
        raise ValueError(f"membership map missing columns: {sorted(missing)}")
    plt = _pyplot()
    figure, axis = plt.subplots(figsize=(6, 5))
    changed = assignments["unit_changed"].astype(bool)
    axis.scatter(assignments.loc[~changed, "x"], assignments.loc[~changed, "y"], s=7,
                 color="#aeb7c2", alpha=.55, label="membership 一致", rasterized=True)
    axis.scatter(assignments.loc[changed, "x"], assignments.loc[changed, "y"], s=9,
                 color="#b43c55", alpha=.9, label="membership 变化", rasterized=True)
    axis.set(title=title, xlabel="Raw x", ylabel="Raw y", aspect="equal")
    axis.legend(markerscale=2)
    return _save(figure, destination)


def plot_window_field(windows: pd.DataFrame, destination: Path, *, value_column: str, title: str) -> Path:
    """Map a continuous field beside its accepted region/support mask.

    A threshold failure is represented by missing in_region values. Those
    values stay visible in the continuous panel and are never coerced to true.
    The second panel shows the accepted mask and marks valid-but-unclassified
    windows in grey.
    """
    required = {"window_x", "window_y", value_column}
    if missing := required - set(windows.columns):
        raise ValueError(f"window field missing columns: {sorted(missing)}")
    plt = _pyplot()
    values = pd.to_numeric(windows[value_column], errors="coerce")
    display_title = title.replace("_", " ").title() if "_" in title else title
    if "in_region" not in windows:
        figure, axis = plt.subplots(figsize=(6, 5))
        _plot_continuous_field(axis, figure, windows, values, value_column, display_title)
        return _save(figure, destination)

    figure, axes = plt.subplots(1, 2, figsize=(11, 5), squeeze=False)
    figure.suptitle(display_title)
    continuous_axis, mask_axis = axes[0]
    _plot_continuous_field(
        continuous_axis, figure, windows, values, value_column,
        "连续场",
    )

    selected = _true_mask(windows["in_region"])
    unsupported = windows["in_region"].isna()
    if "valid_window" in windows:
        supported = _true_mask(windows["valid_window"])
    else:
        supported = values.notna()
    if supported.any():
        mask_axis.scatter(
            windows.loc[supported, "window_x"], windows.loc[supported, "window_y"],
            s=38, color="#b7bec8", alpha=.8, label="支持充分的窗口", rasterized=True,
        )
    if selected.any():
        mask_axis.scatter(
            windows.loc[selected, "window_x"], windows.loc[selected, "window_y"],
            s=60, color="#c43c39", alpha=.9, label="已接受区域", rasterized=True,
        )
    unknown = unsupported & supported
    if unknown.any():
        mask_axis.scatter(
            windows.loc[unknown, "window_x"], windows.loc[unknown, "window_y"],
            s=78, color="#6b7280", marker="x", linewidths=1.2,
            label="阈值不可用", rasterized=True,
        )
    if not selected.any():
        message = "没有已接受区域"
        if unknown.any():
            message += "\n阈值不可用"
        _no_values(mask_axis, message)
    mask_axis.set(
        title="已接受区域掩膜",
        xlabel="窗口 x", ylabel="窗口 y", aspect="equal",
    )
    if supported.any() or selected.any() or unknown.any():
        mask_axis.legend(markerscale=1)
    return _save(figure, destination)


def plot_support_curves(
    support: pd.DataFrame,
    destination: Path,
    *,
    x_column: str | None = None,
    y_column: str | None = None,
    group_column: str | None = None,
    title: str = "候选尺度的窗口支持",
) -> Path:
    """Plot saved support curves without selecting a scale.

    The workflow may call the x/y columns window_side_microns and
    retained_fraction or use equivalent names. Raw, SVC, and common-valid
    curves are kept as separate series; the SVC series is visually primary
    when it is present. A flat or unavailable curve is still drawn when its
    values are saved and never receives an inferred knee here.
    """
    x_column = x_column or _first_column(
        support, ("window_side_microns", "window_side_length", "side_length", "scale", "candidate")
    )
    y_column = y_column or _first_column(
        support, (
            "n_valid_windows",
            "retained_fraction", "support_fraction", "fraction_retained",
            "valid_fraction", "n_valid_fraction", "retained_parent_unit_fraction",
            "valid_window_fraction", "common_fraction_of_svc_valid",
        )
    )
    if x_column is None or y_column is None:
        raise ValueError("support table requires candidate-scale and retained-fraction columns")
    group_column = group_column or _first_column(
        support, ("curve", "basis", "support_kind", "side", "source", "comparison", "label")
    )
    frame = support.copy()
    frame["_support_x"] = pd.to_numeric(frame[x_column], errors="coerce")
    frame["_support_y"] = pd.to_numeric(frame[y_column], errors="coerce")
    frame = frame.loc[frame["_support_x"].notna() & frame["_support_y"].notna()].copy()
    plt = _pyplot()
    figure, axis = plt.subplots(figsize=(7, 4.5))
    if frame.empty:
        _no_values(axis, "没有可计算的支持曲线值")
    else:
        if group_column is None:
            groups = [("", frame)]
        else:
            groups = list(frame.groupby(group_column, dropna=False, sort=True))
        for group, values in groups:
            raw_label = str(group) if group != "" else "support"
            normalized = raw_label.lower().replace("_", "-").replace(" ", "-")
            label = _support_label(raw_label)
            primary = normalized in {"svc", "svc-parent", "svc-parent-support", "common-svc"}
            values = values.sort_values("_support_x")
            axis.plot(
                values["_support_x"], values["_support_y"],
                marker="s" if "common" in normalized else "o",
                markersize=9 if "common" in normalized else 6,
                markerfacecolor="none" if "common" in normalized else None,
                linewidth=2.5 if primary else 1.6,
                alpha=1.0 if primary else .82,
                label=label,
                linestyle="-" if primary else "--",
            )
    axis.set(title=title, xlabel=_axis_label(x_column), ylabel=_axis_label(y_column))
    if not frame.empty and group_column is not None:
        axis.legend(title="支持来源")
    return _save(figure, destination)


def plot_support_curve(support: pd.DataFrame, destination: Path, **kwargs) -> Path:
    """Singular compatibility alias for plot_support_curves."""
    return plot_support_curves(support, destination, **kwargs)


def plot_window_support_grid(
    windows: pd.DataFrame,
    destination: Path,
    *,
    title: str = "当前 SVC parent 网格支持",
) -> Path:
    """Plot the saved configured grid with its actual coordinate side length."""
    required = {"window_x", "window_y", "n_units", "valid_window", "window_side_length"}
    if missing := required - set(windows):
        raise ValueError(f"window support grid missing columns: {sorted(missing)}")
    frame = windows.copy()
    for column in ("window_x", "window_y", "n_units", "window_side_length"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.loc[
        frame[["window_x", "window_y", "n_units", "window_side_length"]].notna().all(axis=1)
        & frame["window_side_length"].gt(0)
    ].copy()
    if frame.empty:
        raise ValueError("window support grid has no finite positive saved windows")
    sides = frame["window_side_length"].unique()
    if len(sides) != 1:
        raise ValueError("window support grid must use one configured coordinate side length")
    side = float(sides[0])
    figure, axis = _pyplot().subplots(figsize=(6.5, 5.5))
    _square_values(axis, figure, frame, frame["n_units"], side, "每窗单位数")
    valid = _true_mask(frame["valid_window"])
    invalid = ~valid
    if invalid.any():
        axis.scatter(
            frame.loc[invalid, "window_x"], frame.loc[invalid, "window_y"], marker="x",
            color="#4B5563", linewidths=.8, s=18, label="支持不足",
        )
        axis.legend()
    x = frame["window_x"]
    y = frame["window_y"]
    axis.set(
        title=title, xlabel="窗口 x", ylabel="窗口 y", aspect="equal",
        xlim=(float(x.min() - side / 2), float(x.max() + side / 2)),
        ylim=(float(y.min() - side / 2), float(y.max() + side / 2)),
    )
    micron_values = pd.to_numeric(frame.get("window_side_microns", pd.Series(dtype=float)), errors="coerce").dropna().unique()
    if len(micron_values) == 1:
        axis.text(.01, .01, f"方窗边长 {micron_values[0]:g} μm", transform=axis.transAxes, fontsize=8, color="#4B5563")
    return _save(figure, destination)


def plot_threshold_bootstrap(
    thresholds: pd.DataFrame,
    destination: Path,
    *,
    value_column: str | None = None,
    title: str = "阈值 bootstrap 分布",
) -> Path:
    """Plot a saved threshold-bootstrap table, including unstable attempts."""
    value_column = value_column or _first_column(
        thresholds, ("bootstrap_threshold", "threshold", "breakpoint", "value")
    )
    if value_column is None:
        raise ValueError("threshold table requires a bootstrap value column")
    values = pd.to_numeric(thresholds[value_column], errors="coerce").dropna()
    figure, axis = _pyplot().subplots()
    if values.empty:
        _no_values(axis, "没有可计算的 bootstrap 阈值")
    else:
        axis.hist(values, bins=min(30, max(5, len(values))), color="#4C78A8", alpha=.85)
    axis.set(title=title, xlabel=_axis_label(value_column), ylabel="bootstrap 抽样次数")
    return _save(figure, destination)


def plot_anatomy_composition(
    context: pd.DataFrame,
    destination: Path,
    *,
    category_column: str | None = None,
    title: str = "Anatomy 窗口组成",
) -> Path:
    """Plot counts of saved anatomy context categories.

    This is a descriptive count of the saved point/window assignments. It does
    not infer labels or turn a missing category into zero.
    """
    category_column = category_column or _first_column(
        context, ("level1_region", "anatomy_class", "broad_label")
    )
    if category_column is None:
        raise ValueError("anatomy context requires a category column")
    counts = context[category_column].astype("string").fillna("<missing>").value_counts(sort=True)
    figure, axis = _pyplot().subplots(figsize=(6, 4))
    if counts.empty:
        _no_values(axis, "没有保存 Anatomy 组成")
    else:
        counts.sort_values().plot.barh(ax=axis, color="#ECA82C")
    axis.set(title=title, xlabel="已保存点/窗口数", ylabel=_axis_label(category_column))
    return _save(figure, destination)


def plot_baseline_window_metrics(
    tables: dict[str, pd.DataFrame],
    destination: Path,
    *,
    title: str,
) -> Path:
    """Compare four saved diversity metrics without forcing one shared cohort.

    The left column uses each label system's own valid native windows.  The
    right column computes a separate common-valid window set for every
    baseline against the SVC reconstruction table.  It never requires all
    baselines to share one intersection.
    """
    if "svc_reconstruction" not in tables:
        raise ValueError("baseline metric view requires an SVC reconstruction table")
    metrics = ("k_obs", "entropy", "neff", "evenness")
    if not all(metric in tables["svc_reconstruction"] for metric in metrics):
        raise ValueError("SVC reconstruction table lacks the four diversity metrics")
    colors = {
        "svc_reconstruction": "#C44E52", "raw_leiden": "#4C78A8",
        "raw_level2": "#59A14F", "svc_leiden": "#F28E2B", "raw_k_control": "#B279A2",
    }
    labels = {
        "svc_reconstruction": "SVC reconstruction", "raw_leiden": "Raw Leiden",
        "raw_level2": "Raw Level2", "svc_leiden": "SVC Leiden", "raw_k_control": "Raw K-control",
    }
    svc = tables["svc_reconstruction"]
    comparison_available = any(
        key != "svc_reconstruction"
        and any(metric in frame and not _common_valid_metric(frame, svc, metric).empty for metric in metrics)
        for key, frame in tables.items()
    )
    columns = 2 if comparison_available else 1
    figure, axes = _pyplot().subplots(
        len(metrics), columns, figsize=(12 if columns == 2 else 6.5, 3.2 * len(metrics)), squeeze=False,
    )
    for row, metric in enumerate(metrics):
        native_values, native_labels, native_colors = [], [], []
        for key, frame in tables.items():
            if metric not in frame:
                continue
            mask = _valid_metric_mask(frame, metric)
            values = pd.to_numeric(frame.loc[mask, metric], errors="coerce").dropna().to_numpy()
            if values.size:
                native_values.append(values)
                native_labels.append(labels.get(key, key))
                native_colors.append(colors.get(key, "#7F8791"))
        if native_values:
            boxes = axes[row, 0].boxplot(native_values, tick_labels=native_labels, patch_artist=True, showfliers=False)
            for box, color in zip(boxes["boxes"], native_colors, strict=True):
                box.set_facecolor(color)
                box.set_alpha(.65)
        else:
            _no_values(axes[row, 0], f"没有可计算的 {_axis_label(metric)}")
        axes[row, 0].set(title=f"{_axis_label(metric)}：各自 native 有效窗口", ylabel=_axis_label(metric))
        axes[row, 0].tick_params(axis="x", rotation=18)

        if not comparison_available:
            continue
        deltas, delta_labels, delta_colors = [], [], []
        for key, frame in tables.items():
            if key == "svc_reconstruction" or metric not in frame:
                continue
            paired = _common_valid_metric(frame, svc, metric)
            if paired.empty:
                continue
            deltas.append((paired["svc"] - paired["baseline"]).to_numpy())
            delta_labels.append(labels.get(key, key))
            delta_colors.append(colors.get(key, "#7F8791"))
        if deltas:
            boxes = axes[row, 1].boxplot(deltas, tick_labels=delta_labels, patch_artist=True, showfliers=False)
            for box, color in zip(boxes["boxes"], delta_colors, strict=True):
                box.set_facecolor(color)
                box.set_alpha(.65)
            axes[row, 1].axhline(0, color="#777", linewidth=1, linestyle="--")
        else:
            _no_values(axes[row, 1], f"没有与 SVC 共同有效的 {_axis_label(metric)} 窗口")
        axes[row, 1].set(
            title=f"{_axis_label(metric)}：各 baseline 分别与 SVC 取共同有效窗口",
            ylabel=f"SVC − baseline {_axis_label(metric)}",
        )
        axes[row, 1].tick_params(axis="x", rotation=18)
    figure.suptitle(title)
    figure.tight_layout()
    return _save(figure, destination)


def plot_anatomy_neff_state(
    anatomy: pd.DataFrame,
    state: pd.DataFrame,
    destination: Path,
    *,
    anatomy_window_side: float | None,
    parent_window_side: float | None,
    title: str,
) -> Path:
    """Draw SVC-defined Anatomy, SVC Neff, and State on one shared coordinate extent."""
    anatomy_required = {"window_x", "window_y", "level1_region"}
    state_required = {"window_x", "window_y", "neff"}
    if missing := anatomy_required - set(anatomy):
        raise ValueError(f"anatomy windows missing columns: {sorted(missing)}")
    if missing := state_required - set(state):
        raise ValueError(f"state windows missing columns: {sorted(missing)}")
    figure, axes = _pyplot().subplots(1, 3, figsize=(15, 5), squeeze=False)
    axes = axes[0]
    colors = {"Tumor": "#C44E52", "Normal": "#59A14F", "Interface": "#F2A541", "Other": "#87929E"}
    _square_categories(axes[0], anatomy, "level1_region", anatomy_window_side, colors)
    axes[0].set_title("SVC Anatomy 方窗")

    neff = pd.to_numeric(state["neff"], errors="coerce")
    _square_values(axes[1], figure, state, neff, parent_window_side, "Neff")
    missing_neff = neff.isna()
    if missing_neff.any():
        _draw_squares_or_centers(
            axes[1], state.loc[missing_neff], parent_window_side,
            color="#D1D5DB", label="Neff 缺失（支持不足）", values=None,
        )
        from matplotlib.lines import Line2D
        axes[1].legend(handles=[Line2D(
            [0], [0], marker="s", color="none", markerfacecolor="#D1D5DB",
            markeredgecolor="none", label="Neff 缺失（支持不足）", markersize=9,
        )])
    axes[1].set_title("SVC Neff 连续场")

    if "in_region" in state:
        known = state["in_region"].notna()
        selected = _true_mask(state["in_region"])
        categories = pd.Series("支持不足", index=state.index, dtype="string")
        categories.loc[known & ~selected] = "区域外"
        categories.loc[selected] = "区域内"
        categories.loc[~known & neff.notna()] = "区域状态未知"
        _square_status(axes[2], state, categories, parent_window_side)
    else:
        _no_values(axes[2], "未保存 State 区域状态")
    axes[2].set_title("State 区域掩膜")

    all_x = pd.concat([pd.to_numeric(anatomy["window_x"], errors="coerce"), pd.to_numeric(state["window_x"], errors="coerce")]).dropna()
    all_y = pd.concat([pd.to_numeric(anatomy["window_y"], errors="coerce"), pd.to_numeric(state["window_y"], errors="coerce")]).dropna()
    half = max([value / 2 for value in (anatomy_window_side, parent_window_side) if value is not None] or [0.0])
    if not all_x.empty and not all_y.empty:
        extent = (float(all_x.min() - half), float(all_x.max() + half), float(all_y.min() - half), float(all_y.max() + half))
        for axis in axes:
            axis.set(xlim=extent[:2], ylim=extent[2:], xlabel="x", ylabel="y", aspect="equal")
    unknown = []
    if anatomy_window_side is None:
        unknown.append("Anatomy 方窗边长")
    if parent_window_side is None:
        unknown.append("parent 方窗边长")
    note = f"；尺度未知：{'、'.join(unknown)}（仅显示保存中心）" if unknown else ""
    figure.suptitle(title + note)
    figure.tight_layout()
    return _save(figure, destination)


def _valid_metric_mask(frame: pd.DataFrame, metric: str) -> pd.Series:
    mask = pd.to_numeric(frame[metric], errors="coerce").notna()
    if "valid_window" in frame:
        mask &= _true_mask(frame["valid_window"])
    return mask


def _common_valid_metric(baseline: pd.DataFrame, svc: pd.DataFrame, metric: str) -> pd.DataFrame:
    if "window_id" not in baseline or "window_id" not in svc:
        return pd.DataFrame(columns=["baseline", "svc"])
    left = baseline.loc[_valid_metric_mask(baseline, metric), ["window_id", metric]].copy()
    right = svc.loc[_valid_metric_mask(svc, metric), ["window_id", metric]].copy()
    left["window_id"], right["window_id"] = left["window_id"].astype(str), right["window_id"].astype(str)
    return left.merge(right, on="window_id", suffixes=("_baseline", "_svc")).rename(
        columns={f"{metric}_baseline": "baseline", f"{metric}_svc": "svc"}
    )[["baseline", "svc"]]


def _square_categories(axis, frame, column, side, colors) -> None:
    from matplotlib.lines import Line2D
    handles = []
    for label, subset in frame.groupby(column, dropna=False, sort=True):
        color = colors.get(str(label), "#87929E")
        _draw_squares_or_centers(axis, subset, side, color=color, label=str(label), values=None)
        handles.append(Line2D([0], [0], marker="s", color="none", markerfacecolor=color,
                              markeredgecolor="none", label=str(label), markersize=9))
    if handles:
        axis.legend(handles=handles, title="Anatomy")


def _square_values(axis, figure, frame, values, side, label, *, cmap="viridis") -> None:
    valid = values.notna()
    if not valid.any():
        _no_values(axis, f"没有可计算的 {label}")
        return
    if side is None:
        points = axis.scatter(frame.loc[valid, "window_x"], frame.loc[valid, "window_y"], c=values.loc[valid], cmap=cmap, s=42)
        _inset_colorbar(axis, figure, points, label)
        return
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    patches = [Rectangle((float(x) - side / 2, float(y) - side / 2), side, side)
               for x, y in zip(frame.loc[valid, "window_x"], frame.loc[valid, "window_y"], strict=True)]
    collection = PatchCollection(patches, cmap=cmap, edgecolor="none")
    collection.set_array(values.loc[valid].to_numpy(dtype=float))
    axis.add_collection(collection)
    _inset_colorbar(axis, figure, collection, label)


def _draw_squares_or_centers(axis, frame, side, *, color, label, values) -> None:
    if side is None:
        axis.scatter(frame["window_x"], frame["window_y"], s=42, color=color, label=label)
        return
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    patches = [Rectangle((float(x) - side / 2, float(y) - side / 2), side, side)
               for x, y in zip(frame["window_x"], frame["window_y"], strict=True)]
    collection = PatchCollection(patches, facecolor=color, edgecolor="white", linewidth=.15, label=label)
    axis.add_collection(collection)


def _square_status(axis, frame, categories, side) -> None:
    from matplotlib.lines import Line2D
    colors = {"区域内": "#C44E52", "区域外": "#4C78A8", "区域状态未知": "#7C8796", "支持不足": "#D1D5DB"}
    handles = []
    for category in ("区域内", "区域外", "区域状态未知", "支持不足"):
        selected = categories.eq(category)
        if not selected.any():
            continue
        color = colors[category]
        _draw_squares_or_centers(axis, frame.loc[selected], side, color=color, label=category, values=None)
        handles.append(Line2D([0], [0], marker="s", color="none", markerfacecolor=color,
                              markeredgecolor="none", label=category, markersize=9))
    if handles:
        axis.legend(handles=handles, title="State 状态")


def _inset_colorbar(axis, figure, mappable, label: str) -> None:
    """Keep map axes identical by placing a compact colorbar inside the axes."""
    color_axis = axis.inset_axes([.12, .035, .76, .035])
    colorbar = figure.colorbar(mappable, cax=color_axis, orientation="horizontal")
    colorbar.set_label(label, fontsize=8)
    colorbar.ax.tick_params(labelsize=7, pad=1)


def _first_column(frame: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    return next((column for column in candidates if column in frame.columns), None)


def _true_mask(values: pd.Series) -> pd.Series:
    """Parse bool columns written by pandas, including nullable/CSV forms."""
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    return (
        values.astype("string").str.strip().str.lower().eq("true")
        .fillna(False).astype(bool)
    )


def _plot_continuous_field(axis, figure, windows, values, value_column, title) -> None:
    if values.notna().any():
        scatter = axis.scatter(
            windows["window_x"], windows["window_y"], c=values, s=35,
            cmap="viridis", alpha=.9, rasterized=True,
        )
        figure.colorbar(scatter, ax=axis, label=_axis_label(value_column))
    else:
        _no_values(axis, f"没有可计算的 {_axis_label(value_column)}")
    axis.set(title=title, xlabel="窗口 x", ylabel="窗口 y", aspect="equal")


def _no_values(axis, message: str) -> None:
    axis.text(.5, .5, message, transform=axis.transAxes, ha="center", va="center", color="#667")


def _pyplot():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("matplotlib is required to create figures") from exc
    plt.rcParams["font.sans-serif"] = ["PingFang SC", "Heiti SC", "Arial Unicode MS", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    return plt


def _save(figure, destination: Path) -> Path:
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(destination, dpi=160, bbox_inches="tight")
    _pyplot().close(figure)
    return destination


def plot_reconstruction_labels(units: pd.DataFrame, destination: Path) -> Path:
    """Show supplied labels without assigning biological names."""
    figure, axis = _pyplot().subplots(figsize=(8, 6))
    palette = _pyplot().get_cmap('tab20')
    for number, (label, frame) in enumerate(units.groupby('reconstruction_label', sort=True)):
        axis.scatter(frame.x, frame.y, s=4, alpha=.7, color=palette(number % 20),
                     label=str(label), rasterized=True)
    axis.set(title='输入的重建状态标签', xlabel='x', ylabel='y', aspect='equal')
    axis.legend(title='重建标签', bbox_to_anchor=(1.02, 1), loc='upper left', markerscale=3)
    return _save(figure, destination)


def plot_diversity_distribution(windows: pd.DataFrame, destination: Path, *, value_column: str, title: str) -> Path:
    values = pd.to_numeric(windows[value_column], errors='coerce').dropna()
    figure, axis = _pyplot().subplots(figsize=(6, 4))
    if values.empty:
        _no_values(axis, '没有支持充分的数值')
    else:
        axis.hist(values, bins=min(40, max(5, len(values))), color='#4c78a8')
    axis.set(title=title, xlabel=_axis_label(value_column), ylabel='支持充分的窗口数')
    return _save(figure, destination)


def _axis_label(value: str) -> str:
    labels = {
        "window_side_microns": "窗口边长（μm）",
        "window_side_length": "窗口边长（坐标单位）",
        "retained_fraction": "保留比例",
        "retained_parent_unit_fraction": "保留 parent 单位比例",
        "valid_window_fraction": "有效窗口比例",
        "common_fraction_of_svc_valid": "占 SVC 有效窗口比例",
        "bootstrap_threshold": "bootstrap 阈值",
        "n_valid_windows": "有效窗口数",
        "neff": "Neff",
        "delta_neff_vs_raw_leiden": "相对 Raw Leiden 的 ΔNeff",
    }
    return labels.get(value, value.replace("_", " "))


def _support_label(value: str) -> str:
    normalized = value.lower().replace("_", "-").replace(" ", "-")
    labels = {
        "svc": "SVC parent",
        "svc-parent": "SVC parent",
        "svc-parent-support": "SVC parent",
        "common-svc": "共同 SVC 支持",
        "raw-baseline": "Raw baseline",
        "raw-coordinate-potential": "Raw 坐标潜在支持",
        "common-valid": "共同有效窗口",
        "common-coordinate-potential": "共同坐标潜在支持",
        "support": "支持",
    }
    return labels.get(normalized, value)
