"""Reusable impact figures from tabular workflow outputs."""
from __future__ import annotations

from pathlib import Path
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
    """Draw the saved full-Raw Level1 anatomy context in native coordinates."""
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
    axis.set(title="Raw Level1 Anatomy 窗口背景", xlabel="x", ylabel="y", aspect="equal")
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
