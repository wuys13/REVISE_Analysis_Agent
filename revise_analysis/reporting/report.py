"""A deliberately modest static report; it never loads an AnnData object."""
from __future__ import annotations

from html import escape
from pathlib import Path
import json

import pandas as pd


def render_report(result_dir: Path) -> Path:
    """Render ``report.html`` from an already-written ``result.json`` and artifacts."""
    result_dir = Path(result_dir)
    result_path = result_dir / "result.json"
    result = json.loads(result_path.read_text(encoding="utf-8"))
    outputs = result.get("outputs", {})
    if not isinstance(outputs, dict):
        raise ValueError("result.json outputs must be a mapping")
    artifacts = _artifacts(result_dir, outputs)
    report_path = result_dir / "report.html"
    report_path.write_text(_document(result, artifacts), encoding="utf-8")
    return report_path


def _artifacts(result_dir: Path, outputs: dict) -> list[dict]:
    result = []
    for label, relative in sorted(outputs.items()):
        if not isinstance(relative, str):
            continue
        if Path(relative).as_posix() == "report.html":
            continue
        path = (result_dir / relative).resolve()
        if result_dir.resolve() not in path.parents or not path.exists():
            continue
        entry = {"label": label, "relative": relative, "suffix": path.suffix.lower(), "path": path}
        if path.suffix.lower() == ".csv":
            entry["preview"] = _csv_preview(path)
        result.append(entry)
    return result


def _csv_preview(path: Path) -> str:
    try:
        table = pd.read_csv(path, nrows=8)
    except (OSError, pd.errors.ParserError, UnicodeDecodeError):
        return ""
    if table.empty:
        return "<p class=\"muted\">Empty table.</p>"
    headers = "".join(f"<th>{escape(str(value))}</th>" for value in table.columns)
    rows = "".join(
        "<tr>" + "".join(f"<td>{escape(str(value))}</td>" for value in row) + "</tr>"
        for row in table.itertuples(index=False, name=None)
    )
    return f"<table><thead><tr>{headers}</tr></thead><tbody>{rows}</tbody></table>"


def _document(result: dict, artifacts: list[dict]) -> str:
    status = escape(str(result.get("status", "unknown")))
    sample_id = escape(str(result.get("sample_id", "unspecified sample")))
    analysis = escape(str(result.get("analysis", "saved analysis")))
    parameters = escape(json.dumps(result.get("parameters", {}), indent=2, sort_keys=True))
    unavailable = result.get("unavailable", [])
    limitations = "".join(
        f"<li><code>{escape(str(item.get('component', 'component')))}</code>: {escape(str(item.get('reason', 'unavailable')))}</li>"
        for item in unavailable if isinstance(item, dict)
    ) or "<li>None recorded.</li>"
    groups = [
        ("overall", "What representation does each side contain?"),
        ("membership", "Which Raw-defined cohorts changed membership?"),
        ("spatial", "Where are anatomy, local State, and positive Gain observed?"),
        ("genes", "What native spatial autocorrelation and pathway activity were saved?"),
    ]
    sections, navigation = [], []
    for key, heading in groups:
        selected = [artifact for artifact in artifacts if _group(artifact["relative"]) == key]
        if not selected:
            continue
        anchor = f"section-{key}"
        navigation.append(f'<li><a href="#{anchor}">{escape(heading)}</a></li>')
        figures = [artifact for artifact in selected if artifact["suffix"] in {".png", ".jpg", ".jpeg", ".svg"}]
        evidence = "".join(_evidence(artifact) for artifact in selected if artifact not in figures)
        figure_html = "".join(f'<figure><img src="{escape(item["relative"])}" alt="{escape(item["label"])}"><figcaption><a href="{escape(item["relative"])}">{escape(item["label"])}</a></figcaption></figure>' for item in figures)
        sections.append(f'<section id="{anchor}"><h2>{escape(heading)}</h2>{figure_html}<details><summary>Saved evidence</summary>{evidence}</details></section>')
    observations = "".join(f"<li>{escape(item)}</li>" for item in _observations(artifacts))
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{sample_id} — {analysis} report</title>
<style>body{{font:16px/1.5 system-ui,sans-serif;max-width:1100px;margin:auto;padding:2rem;color:#20242a}}nav{{position:sticky;top:0;background:#fff;padding:.5rem 0;border-bottom:1px solid #ddd}}section{{margin:2.5rem 0}}table{{border-collapse:collapse;max-width:100%;display:block;overflow:auto}}th,td{{border:1px solid #ddd;padding:.35rem .55rem;text-align:left}}th{{background:#f6f7f9}}pre{{background:#f6f7f9;padding:1rem;overflow:auto}}img{{max-width:100%;height:auto}}figure{{margin:1rem 0}}.muted{{color:#667}}</style>
</head><body><main><h1>{sample_id}: {analysis}</h1><p>Status: <strong>{status}</strong>. This page reports saved computations and does not establish biological interpretation or causal conclusions.</p>
<section><h2>Saved factual observations</h2><ul>{observations or '<li>No readable summary table was saved.</li>'}</ul></section>
<nav><strong>Questions</strong><ol>{''.join(navigation) or '<li>No readable artifacts.</li>'}</ol></nav>
<details><summary>Effective parameters</summary><pre>{parameters}</pre></details>
<details><summary>Unavailable components and limits</summary><ul>{limitations}</ul></details>
{''.join(sections)}</main></body></html>"""


def _group(relative: str) -> str:
    value = relative.lower()
    if "membership" in value or "changed_map" in value:
        return "membership"
    if any(term in value for term in ("moran", "pathway", "program")):
        return "genes"
    if any(term in value for term in ("anatomy", "window", "state", "gain")):
        return "spatial"
    return "overall"


def _evidence(artifact: dict) -> str:
    label, relative = escape(artifact["label"]), escape(artifact["relative"])
    body = artifact.get("preview", "") or f'<p><a href="{relative}">Open saved artifact</a></p>'
    return f'<article><h3><a href="{relative}">{label}</a></h3>{body}</article>'


def _observations(artifacts: list[dict]) -> list[str]:
    """State only values directly present in compact saved tables."""
    messages = []
    for artifact in artifacts:
        relative, path = artifact["relative"], artifact["path"]
        if relative.endswith("overview.csv"):
            table = _read_table(path)
            if table is not None and {"side", "n_units", "n_genes"}.issubset(table):
                detail = ", ".join(
                    f"{row.side}: {int(row.n_units)} units and {int(row.n_genes)} genes"
                    for row in table.itertuples(index=False)
                )
                messages.append(f"The independent overview records {detail}.")
        elif relative.endswith("partition_summary.csv"):
            table = _read_table(path)
            if table is not None and {"side", "scope", "n_clusters"}.issubset(table):
                messages.append("The partition summary records native Leiden cluster counts for each available side and scope.")
        elif "region_extent" in relative:
            table = _read_table(path)
            if table is not None and not table.empty and {"n_region_windows", "n_valid_windows"}.issubset(table):
                row = table.iloc[0]
                messages.append(
                    f"{artifact['label']} records {int(row['n_region_windows'])} selected windows among "
                    f"{int(row['n_valid_windows'])} valid windows."
                )
    return list(dict.fromkeys(messages))


def _read_table(path: Path) -> pd.DataFrame | None:
    try:
        return pd.read_csv(path)
    except (OSError, pd.errors.ParserError, UnicodeDecodeError):
        return None
