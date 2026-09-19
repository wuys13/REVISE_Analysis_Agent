"""Deterministic reader-facing facts for saved reconstruction-impact results.

This module only reads ``result.json`` and output files declared by it.  The
facts are deliberately modest: they expose counts, denominators, and signed
comparisons already present in saved artifacts without opening AnnData or
recomputing an analysis.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def artifact_anchor(relative: str) -> str:
    """Return the stable HTML anchor used for one declared artifact."""
    value = "".join(char if char.isalnum() else "-" for char in str(relative)).strip("-")
    return f"artifact-{value or 'saved-output'}"


def scope_slug(scope: str) -> str:
    """Mirror the analysis artifact slug contract without importing it."""
    return "".join(char if char.isalnum() else "_" for char in str(scope)).strip("_") or "scope"


def configured_scopes(result: dict) -> list[str]:
    """Return scopes in configuration order, with ``All`` left independent."""
    parameters = result.get("parameters") if isinstance(result.get("parameters"), dict) else {}
    values = parameters.get("scopes", [])
    if isinstance(values, (str, bytes)) or not isinstance(values, (list, tuple)):
        return []
    return [str(value) for value in values]


def declared_output_paths(result: dict) -> set[str]:
    outputs = result.get("outputs") if isinstance(result.get("outputs"), dict) else {}
    return {Path(value).as_posix() for value in outputs.values() if isinstance(value, str)}


def reading_summary(result_dir: Path, result: dict | None = None) -> list[dict]:
    """Return at most five deterministic, evidence-linked saved facts.

    Every item has exactly ``text`` and ``anchor``.  A fact is omitted when its
    required file or columns are unavailable; absence is never described as no
    change.
    """
    result_dir = Path(result_dir)
    if result is None:
        result = json.loads((result_dir / "result.json").read_text(encoding="utf-8"))
    declared = declared_output_paths(result)
    facts: list[dict[str, str]] = []

    overview = _read_csv(result_dir, declared, "tables/input_overview.csv")
    if overview is not None and not overview.empty and {"side", "n_units", "n_genes"}.issubset(overview):
        parts = []
        for row in overview.itertuples(index=False):
            if _finite(row.n_units) and _finite(row.n_genes):
                parts.append(f"{str(row.side).upper()} {int(row.n_units)} 个单位、{int(row.n_genes)} 个原生基因")
        if parts:
            _append(facts, "；".join(parts) + "。对象范围独立，总数差不表示单位丢失。", "tables/input_overview.csv")

    labels = _read_csv(result_dir, declared, "tables/reconstruction_label_summary.csv")
    if labels is not None and not labels.empty and "n_units" in labels:
        counts = pd.to_numeric(labels["n_units"], errors="coerce")
        coverage = (f"覆盖 {int(counts.sum())} 个单位，" if counts.notna().all() else "单位总数未知，")
        _append(
            facts,
            f"SVC 重建状态标签{coverage}分为 {len(labels)} 类。",
            "tables/reconstruction_label_summary.csv",
        )

    scopes = configured_scopes(result)
    for scope in scopes:
        relative = f"tables/state_{scope_slug(scope)}_region_extent.csv"
        extent = _read_csv(result_dir, declared, relative)
        if extent is None or extent.empty or not {"n_valid_windows", "n_region_windows"}.issubset(extent):
            continue
        row = extent.iloc[0]
        if _finite(row.n_valid_windows) and _finite(row.n_region_windows):
            _append(facts, f"{scope} 的 {int(row.n_valid_windows)} 个有效窗口中，{int(row.n_region_windows)} 个被定义为 State；支持不足单位单独保留为未知。", relative)
            break
    for scope in scopes:
        relative = f"tables/membership_summary_{scope_slug(scope)}.json"
        payload = _read_json(result_dir, declared, relative)
        summary = payload.get("summary") if isinstance(payload, dict) else None
        row = summary[0] if isinstance(summary, list) and summary and isinstance(summary[0], dict) else None
        if row and _finite(row.get("n_units")) and _finite(row.get("st_unit_change_fraction")):
            _append(
                facts,
                f"{scope} 的共同 ID cohort 比较了 {int(row['n_units'])} 个单位；匹配后的 membership 变化比例为 {float(row['st_unit_change_fraction']):.1%}。",
                relative,
            )
            break

    for scope in scopes:
        relative = f"tables/moran_shared_genes_{scope_slug(scope)}.csv"
        table = _read_csv(result_dir, declared, relative)
        if table is None or table.empty or "comparison_available" not in table or "delta_moran" not in table:
            continue
        available = table.loc[_as_bool(table["comparison_available"])]
        values = pd.to_numeric(available["delta_moran"], errors="coerce").dropna()
        if not values.empty:
            _append(
                facts,
                f"{scope} 有 {len(values)} 个共同可计算基因；保存的 SVC−Raw Moran 差值中位数为 {values.median():.3g}。",
                relative,
            )
            break

    for scope in scopes:
        relative = f"tables/integrated_region_anatomy_summary_{scope_slug(scope)}.csv"
        table = _read_csv(result_dir, declared, relative)
        if table is None or table.empty or not {"region_kind", "in_region", "n_units", "denominator_units"}.issubset(table):
            continue
        state = table.loc[(table["region_kind"].astype(str) == "state") & _as_bool(table["in_region"])]
        if state.empty:
            continue
        units = pd.to_numeric(state["n_units"], errors="coerce")
        denominators = pd.to_numeric(state["denominator_units"], errors="coerce").dropna()
        if not units.notna().all() or denominators.empty:
            continue
        n_units = int(units.sum())
        denominator = int(denominators.max())
        _append(
            facts,
            f"{scope} 的 State 区域包含 {n_units}/{denominator} 个 scope 内 SVC 单位；其 Anatomy 组成可在关联章节展开。",
            relative,
        )
        break

    return facts[:5]


def _append(facts: list[dict[str, str]], text: str, relative: str) -> None:
    facts.append({"text": text, "anchor": artifact_anchor(relative)})


def _inside(root: Path, relative: str) -> Path | None:
    candidate = (root / relative).resolve()
    resolved = root.resolve()
    if resolved not in candidate.parents or not candidate.is_file():
        return None
    return candidate


def _read_csv(root: Path, declared: set[str], relative: str) -> pd.DataFrame | None:
    if relative not in declared or (path := _inside(root, relative)) is None:
        return None
    try:
        return pd.read_csv(path)
    except (OSError, pd.errors.ParserError, UnicodeDecodeError):
        return None


def _read_json(root: Path, declared: set[str], relative: str) -> object | None:
    if relative not in declared or (path := _inside(root, relative)) is None:
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return None


def _as_bool(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values.dtype):
        return values.fillna(False)
    return values.astype(str).str.lower().isin({"true", "1", "yes"})


def _finite(value: object) -> bool:
    try:
        return bool(pd.notna(value)) and float(value) not in (float("inf"), float("-inf"))
    except (TypeError, ValueError):
        return False


__all__ = [
    "artifact_anchor",
    "configured_scopes",
    "declared_output_paths",
    "reading_summary",
    "scope_slug",
]
