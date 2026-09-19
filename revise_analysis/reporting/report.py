"""Render a static report from saved workflow results only.

The report is deliberately a reader for result.json and its declared
artifacts. It never opens an input H5AD or reruns a scientific computation.
New reconstruction-impact results can declare their own ordered sections; the
filename-based grouping remains as a compatibility path for older records.
"""
from __future__ import annotations

from html import escape
from pathlib import Path
import json

import pandas as pd

from .impact_reading import artifact_anchor, configured_scopes, reading_summary, scope_slug


_DEFAULT_GROUPS = (
    ("overall", "两侧分别保存了什么表示？"),
    ("membership", "哪些由 Raw 定义的 cohort 发生了 membership 变化？"),
    ("spatial", "Anatomy、局部 State 和正向 Gain 出现在哪里？"),
    ("genes", "保存了哪些原生空间自相关和通路活性结果？"),
)


# Older result manifests keep English section metadata. Translate those known
# workflow sections at render time so re-rendering is enough to update the
# reader-facing page without changing scientific tables or rerunning analysis.
_SECTION_LOCALIZATION = {
    "input": {
        "Input and reconstruction labels": ("输入与重建标签", "输入能力和已提供的重建状态标签。"),
        "Input capabilities and supplied reconstructed-state labels.": "输入能力和已提供的重建状态标签。",
    },
    "baseline": {
        "Expression baselines": ("表达 baseline", "Raw 独立表达分群和已提供的 Raw Level2 标签。"),
        "Independent Raw expression partitions and supplied Raw Level2 labels.": "Raw 独立表达分群和已提供的 Raw Level2 标签。",
    },
    "support": {
        "Spatial scales and support": ("空间尺度与支持", "与结果无关的支持曲线和尺度推荐。"),
        "Outcome-independent support curves and scale recommendations.": "独立于最终结果的支持曲线和尺度推荐。",
    },
    "diversity": {
        "Local label diversity": ("局部标签多样性", "原生重建标签多样性及可用的 baseline 字段。"),
        "Native reconstructed-label diversity and available baseline fields.": "原生重建标签多样性及可用的 baseline 字段。",
    },
    "regions": {
        "State and delta candidates": ("State 与差值候选", "连续场、阈值诊断和候选掩膜。"),
        "Continuous fields, threshold diagnostics, and candidate masks.": "连续场、阈值诊断和候选掩膜。",
    },
    "anatomy": {
        "Anatomy context": ("Anatomy 背景", "完整 Raw Anatomy 窗口及按观测点组成汇总的 parent 窗口。"),
        "Full-Raw anatomy windows and point-assigned parent-window composition.": "完整 Raw Anatomy 窗口及按观测点组成汇总的 parent 窗口。",
    },
    "molecular": {
        "Molecular spatial results": ("分子空间结果", "可用时保存原生 Moran 结果和固定 cohort 的 AUCell 分数。"),
        "Native Moran results and fixed-cohort AUCell scores when available.": "可用时保存原生 Moran 结果和固定 cohort 的 AUCell 分数。",
    },
    "membership": {
        "Shared-ID membership": ("共同 ID 的 membership", "仅在明确共同 ID cohort 上生成的单位变化证据。"),
        "Optional changed-unit evidence on the explicit common-ID cohort.": "仅在明确共同 ID cohort 上生成的单位变化证据。",
    },
    "integration": {
        "Integrated spatial evidence": ("集成空间证据", "带有明确支持数和分母的派生事实表。"),
        "Derived fact tables with explicit support and denominators.": "带有明确支持数和分母的派生事实表。",
    },
}

_STATUS_LABELS = {
    "succeeded": "已完成",
    "completed": "已完成",
    "partial": "部分完成",
    "skipped": "已跳过",
    "unavailable": "不可用",
    "failed": "失败",
    "error": "执行错误",
    "unknown": "未知",
}

_STAGE_LABELS = {
    "input": "输入",
    "baseline": "baseline",
    "support": "支持诊断",
    "diversity": "局部多样性",
    "regions": "区域",
    "anatomy": "Anatomy",
    "molecular": "分子结果",
    "membership": "membership",
    "integration": "集成证据",
    "figures": "结果图",
}

_ANALYSIS_LABELS = {
    "reconstruction_impact": "重建影响分析",
    "spatial_autocorrelation": "空间自相关分析",
    "pathway_activity": "通路活性分析",
}


def render_report(result_dir: Path) -> Path:
    """Render report.html from an already-written result.json.

    Only files named by result.outputs and present inside result_dir are
    linked. This makes re-rendering safe when a result directory also contains
    notebooks, source inputs, or user notes.
    """
    result_dir = Path(result_dir)
    result_path = result_dir / "result.json"
    result = json.loads(result_path.read_text(encoding="utf-8"))
    outputs = result.get("outputs", {})
    if not isinstance(outputs, dict):
        raise ValueError("result.json outputs must be a mapping")
    artifacts = _artifacts(result_dir, outputs)
    report_path = result_dir / "report.html"
    report_path.write_text(_document(result, artifacts, result_dir), encoding="utf-8")
    return report_path


def _artifacts(result_dir: Path, outputs: dict) -> list[dict]:
    """Build an index of existing, in-result output files.

    outputs is the provenance boundary for the page. A missing declared file
    is retained in the section's missing-reference list, but is never guessed
    from a directory scan.
    """
    root = result_dir.resolve()
    entries = []
    for label, relative in sorted(outputs.items(), key=lambda item: str(item[0])):
        if not isinstance(relative, str):
            continue
        relative_path = Path(relative)
        if relative_path.as_posix() == "report.html":
            continue
        path = (result_dir / relative_path).resolve()
        if root not in path.parents or not path.is_file():
            continue
        entry = {
            "label": str(label),
            "relative": relative_path.as_posix(),
            "suffix": path.suffix.lower(),
            "path": path,
            "size": path.stat().st_size,
        }
        if path.suffix.lower() == ".csv":
            entry["preview"] = _csv_preview(path)
        entries.append(entry)
    return entries


def _csv_preview(path: Path) -> str:
    try:
        table = pd.read_csv(path, nrows=8)
    except (OSError, pd.errors.ParserError, UnicodeDecodeError):
        return ""
    if table.empty:
        return '<p class="muted">表格为空。</p>'
    headers = "".join(f"<th>{escape(str(value))}</th>" for value in table.columns)
    rows = "".join(
        "<tr>" + "".join(f"<td>{escape(str(value))}</td>" for value in row) + "</tr>"
        for row in table.itertuples(index=False, name=None)
    )
    return f"<table><thead><tr>{headers}</tr></thead><tbody>{rows}</tbody></table>"


def _document(result: dict, artifacts: list[dict], result_dir: Path | None = None) -> str:
    if result.get("analysis") == "reconstruction_impact" and result_dir is not None:
        return _impact_document(result_dir, result, artifacts)
    raw_status = str(result.get("status", "unknown"))
    status = escape(raw_status)
    status_label = escape(_STATUS_LABELS.get(raw_status, raw_status))
    sample_id = escape(str(result.get("sample_id", "未指定样本")))
    raw_analysis = str(result.get("analysis", "saved_analysis"))
    analysis = escape(_ANALYSIS_LABELS.get(raw_analysis, raw_analysis))
    parameters = escape(json.dumps(result.get("parameters", {}), indent=2, sort_keys=True, ensure_ascii=False))
    sections = _declared_sections(result, artifacts)
    navigation = "".join(
        f'<li><a href="#section-{escape(section["id"])}">{escape(section["title"])}</a></li>'
        for section in sections
    )

    unavailable = _unavailable_items(result)
    limitations = "".join(
        f"<li><code>{escape(item['component'])}</code>: {escape(item['reason'])}</li>"
        for item in unavailable
    ) or "<li>没有记录。</li>"
    stage_status = _stage_rows(result)
    stage_html = _status_table(stage_status, "没有保存阶段状态记录。")
    file_index = _file_index(artifacts)
    rendered_sections = "".join(_render_section(section) for section in sections)
    observations = "".join(f"<li>{escape(item)}</li>" for item in _observations(artifacts))
    vocabulary = _vocabulary(result)
    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{sample_id} — {analysis} 报告</title>
<style>
body{{font:16px/1.5 system-ui,sans-serif;max-width:1180px;margin:auto;padding:2rem;color:#20242a}}
nav{{background:#fff;padding:.6rem 0;border-bottom:1px solid #ddd}}
nav ol{{display:flex;flex-wrap:wrap;gap:.35rem 1.2rem;padding-left:1.3rem;margin:.5rem 0}}
section{{margin:2.5rem 0}} table{{border-collapse:collapse;max-width:100%;display:block;overflow:auto}}
th,td{{border:1px solid #ddd;padding:.35rem .55rem;text-align:left;vertical-align:top}}
th{{background:#f6f7f9}} pre{{background:#f6f7f9;padding:1rem;overflow:auto}}
img{{max-width:100%;height:auto}} figure{{margin:1rem 0}}.muted{{color:#667}}
.status{{font-weight:650}}.status-partial{{color:#9a6700}}.status-failed{{color:#b42318}}
.file-index td:nth-child(2){{font-family:ui-monospace,monospace;font-size:.9em}}
.boundary{{background:#f6f7f9;border-left:4px solid #6b7280;padding:.8rem 1rem}}
.section-description{{color:#4b5563}} details article{{margin:1rem 0}}
</style>
</head><body><main><h1>{sample_id}：{analysis} 报告</h1>
<p>状态：<strong class="status status-{status}">{status_label}</strong>（{status}）。本页只展示已保存的计算结果，不单独建立生物学解释或因果结论。</p>
{vocabulary}
<section><h2>已保存的事实观察</h2><ul>{observations or '<li>没有保存可阅读的摘要表。</li>'}</ul></section>
<nav><strong>页面分节</strong><ol>{navigation or '<li>没有声明分节。</li>'}</ol></nav>
<section><h2>阶段状态</h2>{stage_html}</section>
<details><summary>最终生效参数</summary><pre>{parameters}</pre></details>
<details><summary>不可用组件与限制</summary><ul>{limitations}</ul></details>
{rendered_sections}
<details><summary>已保存文件索引（Saved file index）</summary>{file_index}</details>
</main></body></html>"""


_IMPACT_QUESTIONS = (
    ("features", "重建后特征与两侧差异", (
        ("population", "分析群体与重建标签"),
        ("rawbaseline", "Raw baseline"),
        ("membership", "共同 ID 的 membership 变化"),
        ("genes", "原生 Moran 与共同基因比较"),
        ("program-overall", "分子 program 总体结果"),
    )),
    ("spatial-state", "状态与差异的空间位置", (
        ("anatomy", "Anatomy 背景"),
        ("scale", "尺度与窗口支持"),
        ("diversity", "局部标签多样性"),
        ("state", "State 连续场与可靠性"),
        ("gain", "Gain 与 Raw K-control"),
    )),
    ("associations", "空间位置上的分子与成员关联", (
        ("region-anatomy", "Region × Anatomy"),
        ("label-composition", "重建标签组成"),
        ("program-location", "Program 位置与汇总"),
        ("changed-location", "变化单位的位置与汇总"),
    )),
)

_TOPIC_CONTEXT = {
    "population": "Raw 与 SVC 的对象范围可不同；这里描述各自表示，不把单位数或标签数的差异解释为改善。重建标签尚不等于生物学亚型。",
    "rawbaseline": "Raw Leiden、已有 Raw Level2 和可选 SVC Leiden 分别保留来源；它们不覆盖输入的重建主标签。",
    "membership": "仅比较明确共同 ID cohort 的分区对应差异，不表示生物学身份转变，也不是其他分析的前置条件。",
    "genes": "Moran 是每侧原生空间图上的并行证据，不依赖 State 区域。完整基因与共享基因结果分别阅读，差值不是配对重建误差。",
    "program-overall": "每侧在固定 cohort 与基因轴上评分一次。查看覆盖率和实际 rank cutoff 后再判断比较条件；两侧分数不自动作差。",
    "anatomy": "完整 Raw 的组织背景独立于 State 网格；Interface 来自窗口内 Tumor 与指定 Normal 来源共存。无覆盖与 Other 分开。",
    "scale": "推荐只依据 SVC 标签支持，不依据 State/Gain 大小；Raw/common 曲线若为坐标潜在支持，不代表已有表达 baseline。显式尺度不会被推荐覆盖。",
    "diversity": "Kobs、entropy、Neff 与 evenness 描述不同方面；主 State 使用完整输入 SVC 标签，表达 baseline 保留自己的实际 cohort。",
    "state": "State 表示 SVC 自身局部标签多样性。先读连续场，再读阈值和区域；阈值不可用时保留连续结果，未知不计入区域外。",
    "gain": "ΔNeff vs Raw Leiden 保留共同有效窗口中的正负差值。正值只作 Gain candidate；Raw K-control 只检查标签数量粒度，不证明生物学改善。",
    "region-anatomy": "按 SVC 观测点落入 Anatomy 的组成阅读，Raw 背景另保留；State、差值候选各有独立分母，不把同源汇总当独立验证。",
    "label-composition": "标签组成帮助解释局部多样性，但 State 本身由这些标签定义，因此区域内外组成不是独立验证。",
    "program-location": "这里将固定单位评分与已有 window、Anatomy、State/Gain 关联，不重新评分。未评分保留缺失；共定位不等于独立验证或因果关系。",
    "changed-location": "分母只包括实际完成共同 ID 比较的单位；区域和 Anatomy 用于定位这些分区对应差异，不反向定义 State。",
}


def _impact_document(result_dir: Path, result: dict, artifacts: list[dict]) -> str:
    """Render the fixed scientific-question tree for impact results."""
    raw_status = str(result.get("status", "unknown"))
    sample_id = escape(str(result.get("sample_id", "未指定样本")))
    by_relative = {artifact["relative"]: artifact for artifact in artifacts}
    scopes = configured_scopes(result) or _scopes_from_exact_outputs(by_relative)
    capability = ""
    overview_artifact = by_relative.get("tables/input_overview.csv")
    overview = _read_table(overview_artifact["path"]) if overview_artifact else None
    if overview is not None and {"side", "expression_available"}.issubset(overview):
        unavailable_sides = overview.loc[~_bool_series(overview.expression_available), "side"].astype(str).str.upper().tolist()
        if unavailable_sides:
            capability = "；" + "、".join(unavailable_sides) + " 的表达分析前提尚未满足，具体限制在相应章节说明"
    facts = reading_summary(result_dir, result)
    fact_html = "".join(
        f'<li><a href="#{escape(item["anchor"])}">{escape(item["text"])}</a></li>'
        for item in facts
    ) or '<li class="muted">没有足够的已保存数值形成顶部事实；请从下方可用性和证据表开始阅读。</li>'
    used_images: set[str] = set()
    body_parts, nav_parts = [], []
    for question_id, question_title, topics in _IMPACT_QUESTIONS:
        topic_parts, topic_nav = [], []
        for topic_id, topic_title in topics:
            blocks = _impact_scope_blocks(topic_id, scopes, result, by_relative)
            scope_html = "".join(
                _render_impact_scope(topic_id, topic_title, block, result, by_relative, used_images,
                                     single_scope=len(scopes) == 1)
                for block in blocks
            )
            topic_anchor = f"topic-{topic_id}"
            topic_parts.append(
                f'<section class="topic" id="{topic_anchor}">{_legacy_aliases(topic_id, result)}'
                f'<h3>{escape(topic_title)}</h3>{scope_html}</section>'
            )
            scope_links = "".join(
                f'<li><a href="#{escape(block["anchor"])}">{escape(block["label"])}</a></li>' for block in blocks
            )
            nested = f"<ol>{scope_links}</ol>" if len(scopes) > 1 else ""
            topic_nav.append(f'<li><a href="#{topic_anchor}">{escape(topic_title)}</a>{nested}</li>')
        body_parts.append(
            f'<section class="question" id="question-{question_id}"><h2>{escape(question_title)}</h2>{"".join(topic_parts)}</section>'
        )
        nav_parts.append(
            f'<li><a href="#question-{question_id}">{escape(question_title)}</a><ol>{"".join(topic_nav)}</ol></li>'
        )

    status_label = escape(_STATUS_LABELS.get(raw_status, raw_status))
    parameters = escape(json.dumps(result.get("parameters", {}), indent=2, sort_keys=True, ensure_ascii=False))
    unavailable = _unavailable_items(result)
    limitations = "".join(
        f'<li><code>{escape(item["component"])}</code>：{escape(item["reason"])}</li>' for item in unavailable
    ) or '<li>没有记录不可用组件。</li>'
    stage_html = _status_table(_stage_rows(result), "没有保存阶段状态记录。")
    legacy = _legacy_metadata(result)
    file_index = _file_index(artifacts)
    manifest = by_relative.get("tables/input_manifest.json")
    manifest_html = (f'<p><a href="{escape(manifest["relative"])}">输入来源、表达声明与临时载体信息</a></p>'
                     if manifest else '<p>本次未保存输入来源清单。</p>')
    nav_parts.insert(0, '<li><a href="#sample-findings">本样本主要发现</a></li>')
    nav_parts.append('<li><a href="#comparison-basis">比较基础与完整证据</a><ol>'
                     '<li><a href="#basis-input">输入来源与参数</a></li>'
                     '<li><a href="#basis-status">执行状态与限制</a></li>'
                     '<li><a href="#basis-files">完整文件索引</a></li></ol></li>')
    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{sample_id} — 重建影响报告</title>
<style>
:root{{--ink:#20242a;--muted:#626b77;--line:#d9dde3;--soft:#f5f7f9;--accent:#22577a}}
*{{box-sizing:border-box}}html{{scroll-behavior:smooth}}body{{margin:0;color:var(--ink);font:16px/1.55 system-ui,sans-serif}}
.layout{{display:grid;grid-template-columns:minmax(230px,290px) minmax(0,980px);gap:2rem;max-width:1320px;margin:auto;padding:1.5rem}}
.toc{{position:sticky;top:1rem;align-self:start;max-height:calc(100vh - 2rem);overflow:auto;border:1px solid var(--line);border-radius:10px;background:#fff}}
.toc summary{{cursor:pointer;padding:.8rem 1rem;font-weight:700}}.toc ol{{padding-left:1.4rem;margin:.35rem 0 .75rem}}.toc ol ol{{font-size:.92em}}
.toc a{{color:#334155;text-decoration:none}}.toc a.active{{color:var(--accent);font-weight:700}}
main{{min-width:0}}a,code,h1,.scope-summary,.scope-peek,.availability{{overflow-wrap:anywhere}}h1{{line-height:1.2}}h2{{margin-top:3.5rem;border-bottom:2px solid var(--line);padding-bottom:.35rem}}
h3{{margin-top:2.5rem}}.lead,.scope-summary{{font-size:1.03rem}}.facts{{background:#fafbfc;border:1px solid var(--line);border-radius:8px;padding:.8rem 1.2rem}}
.scope{{scroll-margin-top:1rem;border:1px solid var(--line);border-radius:10px;padding:1rem;margin:1rem 0;background:#fff}}
.scope>summary{{cursor:pointer;font-size:1.08rem;font-weight:700}}.scope[open]>summary{{margin-bottom:.8rem}}
.scope-peek{{display:block;margin:.25rem 1.25rem 0 0;color:var(--muted);font-size:.88rem;font-weight:400}}.scope[open] .scope-peek{{display:none}}
.availability{{color:#7a4b00;background:#fff8e6;padding:.6rem .8rem;border-radius:6px}}.muted{{color:var(--muted)}}
.main-evidence{{margin:1rem 0}}figure{{margin:1rem 0}}img{{display:block;max-width:100%;height:auto;border:1px solid var(--line)}}
table{{border-collapse:collapse;display:block;max-width:100%;overflow:auto}}th,td{{border:1px solid var(--line);padding:.35rem .55rem;text-align:left;vertical-align:top}}th{{background:var(--soft)}}
.limits{{background:var(--soft);padding:.7rem .9rem;border-radius:6px}}.evidence-list article{{margin:1rem 0}}pre{{background:var(--soft);padding:1rem;overflow:auto}}
.basis{{margin-top:4rem}}.file-index td:nth-child(2){{font-family:ui-monospace,monospace;font-size:.9em}}
@media(max-width:760px){{.layout{{display:block;padding:1rem}}.toc{{position:relative;max-height:none;margin-bottom:1.5rem}}.toc details:not([open]) ol{{display:none}}}}
</style></head><body><div class="layout">
<nav class="toc" aria-label="三层报告目录"><details open><summary>科学问题目录</summary><ol>{"".join(nav_parts)}</ol></details></nav>
<main><h1>{sample_id}：重建影响报告</h1>
<p class="lead">状态：<strong>{status_label}</strong>（{escape(raw_status)}）{escape(capability)}。{("All 是独立 scope，不作为其他 scope 的父级或汇总替代。" if "All" in scopes else "")}</p>
<section class="facts" id="sample-findings"><h2 style="margin-top:.2rem">本样本主要发现</h2><ol>{fact_html}</ol></section>
<details><summary>State 与对照的阅读口径</summary>{_vocabulary(result)}</details>
{"".join(body_parts)}
<section class="basis" id="comparison-basis"><h2>比较依据与完整记录</h2>
<details id="basis-input"><summary>输入来源、表达声明与生效参数</summary>{manifest_html}<pre>{parameters}</pre></details>
<details id="basis-status"><summary>阶段状态与执行错误</summary>{stage_html}</details>
<details><summary>不可用组件、分母限制与原因</summary><ul>{limitations}</ul></details>
{legacy}
<details id="basis-files"><summary>完整已保存文件索引（Saved file index）</summary>{file_index}</details></section>
</main></div>
<script>
const nav=document.querySelector('.toc');
function revealHash(){{const id=decodeURIComponent(location.hash.slice(1));if(!id)return;const target=document.getElementById(id);if(!target)return;let node=target;while(node){{if(node.tagName==='DETAILS')node.open=true;node=node.parentElement;}}}}
function markActive(){{const id=location.hash.slice(1);document.querySelectorAll('.toc a').forEach(a=>a.classList.toggle('active',a.getAttribute('href')==='#'+id));}}
addEventListener('hashchange',()=>{{revealHash();markActive();}});addEventListener('popstate',()=>{{revealHash();markActive();}});revealHash();markActive();
const mobile=matchMedia('(max-width:760px)');const tocDetails=document.querySelector('.toc>details');function setTocForViewport(e){{if(tocDetails)tocDetails.open=!e.matches;}}setTocForViewport(mobile);mobile.addEventListener('change',setTocForViewport);
if('IntersectionObserver'in window){{const links=[...document.querySelectorAll('.toc a[href^="#"]')];const obs=new IntersectionObserver(es=>{{const e=es.filter(x=>x.isIntersecting).sort((a,b)=>a.boundingClientRect.top-b.boundingClientRect.top)[0];if(!e)return;links.forEach(a=>a.classList.toggle('active',a.hash==='#'+e.target.id));}},{{rootMargin:'-10% 0px -75% 0px'}});links.forEach(a=>{{const id=decodeURIComponent(a.hash.slice(1));const t=id?document.getElementById(id):null;if(t)obs.observe(t)}});}}
</script></body></html>"""


def _impact_scope_blocks(topic: str, scopes: list[str], result: dict, by_relative: dict[str, dict]) -> list[dict]:
    blocks: list[dict] = []
    global_paths = _impact_paths(topic, None, result, by_relative)
    if global_paths:
        blocks.append({"scope": None, "label": "对象概览", "anchor": f"scope-{topic}-global", "paths": global_paths})
    if topic in {"population", "program-overall"}:
        if not blocks:
            blocks.append({"scope": None, "label": "可用性", "anchor": f"scope-{topic}-availability", "paths": []})
        return blocks
    for scope in scopes:
        blocks.append({
            "scope": scope,
            "label": scope,
            "anchor": f"scope-{topic}-{scope_slug(scope)}",
            "paths": _impact_paths(topic, scope, result, by_relative),
        })
    if not blocks:
        blocks.append({"scope": None, "label": "可用性", "anchor": f"scope-{topic}-availability", "paths": []})
    return blocks


def _impact_paths(topic: str, scope: str | None, result: dict, by_relative: dict[str, dict]) -> list[str]:
    """Map a topic to exact workflow artifact templates, never loose keywords."""
    slug = scope_slug(scope) if scope is not None else None
    candidates: list[str] = []
    if topic == "population" and scope is None:
        candidates = ["tables/input_overview.csv", "tables/reconstruction_label_summary.csv", "tables/reconstruction_labels_units.csv",
                      "figures/reconstruction_labels.png"]
    elif topic == "rawbaseline" and scope is None:
        candidates = ["tables/raw_level2_baseline_summary.csv", "tables/raw_level2_baseline.csv", "tables/partition_summary.csv", "figures/partition_sizes.png"]
    elif topic == "rawbaseline" and slug:
        candidates = [f"tables/partitions_raw_{slug}.csv", f"tables/partitions_svc_{slug}.csv"]
    elif topic == "membership" and slug:
        candidates = [f"tables/membership_summary_{slug}.json", f"tables/membership_{slug}.csv",
                      f"tables/membership_contingency_{slug}.csv"]
    elif topic == "genes" and slug:
        candidates = [f"tables/moran_shared_genes_{slug}.csv", f"tables/moran_shared_genes_{slug}_metadata.json",
                      f"tables/moran_raw_{slug}.csv", f"tables/moran_svc_{slug}.csv",
                      f"figures/moran_raw_{slug}.png", f"figures/moran_svc_{slug}.png"]
    elif topic == "program-overall" and scope is None:
        for program in _program_slugs(result, by_relative):
            for side in ("raw", "svc"):
                candidates += [f"tables/pathway_{side}_{program}.csv", f"tables/pathway_{side}_{program}_metadata.json",
                               f"figures/pathway_{side}_{program}.png"]
    elif topic == "anatomy" and scope is None:
        candidates = ["tables/raw_anatomy_context.csv", "tables/raw_anatomy_windows.csv",
                      "tables/raw_point_anatomy.csv", "tables/svc_point_anatomy.csv", "figures/raw_anatomy_context.png"]
    elif topic == "anatomy" and slug:
        candidates = [f"tables/parent_anatomy_raw_{slug}.csv", f"tables/parent_anatomy_raw_baseline_{slug}.csv",
                      f"tables/parent_anatomy_svc_{slug}.csv"]
    elif topic == "scale" and slug:
        candidates = [f"tables/window_scale_recommendation_{slug}.json", f"tables/window_scale_support_{slug}.csv",
                      f"figures/window_scale_support_{slug}.png", f"figures/window_scale_support_{slug}_retention.png",
                      f"figures/window_scale_support_{slug}_common_coverage.png",
                      f"figures/state_{slug}_region_windows_support.png"]
    elif topic == "diversity" and slug:
        candidates = [f"tables/window_diversity_svc_{slug}.csv", f"tables/window_diversity_raw_{slug}.csv",
                      f"tables/window_diversity_svc_leiden_baseline_{slug}.csv", f"tables/window_diversity_raw_level2_baseline_{slug}.csv"]
    elif topic == "state" and slug:
        candidates = [f"tables/state_{slug}_region_extent.csv", f"tables/state_{slug}_threshold.json",
                      f"tables/state_{slug}_threshold_bootstrap.csv", f"tables/state_{slug}_region_windows.csv",
                      f"figures/state_{slug}_region_windows.png", f"figures/state_{slug}_region_windows_distribution.png",
                      f"figures/state_{slug}_threshold_bootstrap.png"]
    elif topic == "gain" and slug:
        candidates = [f"tables/gain_{slug}_region_extent.csv", f"tables/gain_{slug}_threshold.json",
                      f"tables/gain_{slug}_threshold_bootstrap.csv", f"tables/gain_{slug}_common_valid_windows.csv",
                      f"tables/gain_{slug}_region_windows.csv", f"tables/raw_k_control_selection_{slug}.json",
                      f"tables/raw_k_control_sweep_{slug}.csv", f"tables/raw_k_control_delta_{slug}.csv",
                      f"tables/window_diversity_raw_k_control_{slug}.csv",
                      f"tables/raw_level2_baseline_vs_svc_{slug}_common_valid_windows.csv",
                      f"figures/gain_{slug}_region_windows.png", f"figures/gain_{slug}_common_valid_windows.png",
                      f"figures/gain_{slug}_threshold_bootstrap.png", f"figures/gain_{slug}_region_windows_distribution.png",
                      f"figures/gain_{slug}_region_windows_support.png"]
    elif topic == "region-anatomy" and slug:
        candidates = [f"tables/integrated_region_anatomy_summary_{slug}.csv", f"tables/integrated_spatial_evidence_{slug}.csv"]
    elif topic == "label-composition" and slug:
        candidates = [f"tables/integrated_label_composition_{slug}.csv", f"tables/integrated_label_units_{slug}.csv"]
    elif topic == "program-location" and slug:
        for program in _program_slugs(result, by_relative):
            for side in ("raw", "svc"):
                prefix = f"tables/integrated_program_{side}_{slug}_{program}"
                candidates += [f"{prefix}_units.csv", f"{prefix}_anatomy.csv", f"{prefix}_region.csv", f"{prefix}_window.csv"]
    elif topic == "changed-location" and slug:
        candidates = [f"figures/spatial_changed_map_{slug}.png", f"tables/integrated_changed_units_{slug}.csv",
                      f"tables/spatial_changed_map_{slug}.csv"]
    present: list[str] = []
    for relative in candidates:
        if relative in by_relative and relative not in present:
            present.append(relative)
        if relative.endswith(".csv"):
            relationship = f"figures/relationships_{Path(relative).stem}.png"
            if relationship in by_relative and relationship not in present:
                present.append(relationship)
    return present


def _program_slugs(result: dict, by_relative: dict[str, dict]) -> list[str]:
    """Resolve program slugs from declared pathway outputs and configured names."""
    configured: list[str] = []
    params = result.get("parameters") if isinstance(result.get("parameters"), dict) else {}
    gene_sets = params.get("gene_sets")
    if isinstance(gene_sets, dict):
        configured.extend(scope_slug(name) for name in gene_sets)
    names = params.get("gene_set_names")
    if isinstance(names, (list, tuple)):
        configured.extend(scope_slug(name) for name in names)
    discovered = []
    for side in ("raw", "svc"):
        prefix, suffix = f"tables/pathway_{side}_", ".csv"
        for relative in by_relative:
            if relative.startswith(prefix) and relative.endswith(suffix):
                discovered.append(relative[len(prefix):-len(suffix)])
    return list(dict.fromkeys([*configured, *sorted(discovered)]))


def _scopes_from_exact_outputs(by_relative: dict[str, dict]) -> list[str]:
    prefix, suffix = "tables/moran_shared_genes_", ".csv"
    return [relative[len(prefix):-len(suffix)] for relative in by_relative
            if relative.startswith(prefix) and relative.endswith(suffix)]


def _render_impact_scope(topic: str, topic_title: str, block: dict, result: dict,
                         by_relative: dict[str, dict], used_images: set[str], *, single_scope: bool) -> str:
    paths = block["paths"]
    artifacts = [by_relative[path] for path in paths]
    table = next((item for item in artifacts if item["suffix"] == ".csv"), None)
    image = next((item for item in artifacts if item["suffix"] in {".png", ".jpg", ".jpeg", ".svg"}
                  and item["relative"] not in used_images), None)
    if image is not None:
        used_images.add(image["relative"])
    summary = _scope_fact(topic, block.get("scope"), block["label"], artifacts, result)
    availability = _topic_unavailable(topic, block.get("scope"), result)
    if not artifacts:
        reason = "；".join(availability) if availability else "本次结果没有保存该主题产物"
        content = f'<p class="availability">不可用：{escape(reason)}。此状态不表示没有变化。</p>'
    else:
        main = ""
        if image is not None:
            main = (f'<figure class="main-evidence" id="{artifact_anchor(image["relative"])}">'
                    f'<img src="{escape(image["relative"])}" alt="{escape(image["label"])}">'
                    f'<figcaption>{escape(topic_title)} · {escape(block["label"])} · <a href="{escape(image["relative"])}">查看原图</a></figcaption></figure>')
        elif table is not None:
            main = (f'<div class="main-evidence" id="{artifact_anchor(table["relative"])}">'
                    f'{table.get("preview", "")}<p><a href="{escape(table["relative"])}">打开完整表</a></p></div>')
        limits = _TOPIC_CONTEXT[topic]
        denominator = _denominator_note(table)
        if denominator:
            limits += " " + denominator
        if availability:
            limits += "；" + "；".join(availability)
        primary = image if image is not None else table
        evidence = "".join(
            _impact_evidence(item, primary=item is primary, used_images=used_images) for item in artifacts
        )
        content = (f'<p class="scope-summary">{escape(summary)}</p>{main}'
                   f'<p class="limits"><strong>分母与限制：</strong>{escape(limits or "以各证据表保存的 cohort 与有效记录为准。")}</p>'
                   f'<details><summary>全部证据与来源</summary><div class="evidence-list">{evidence}</div></details>')
    open_attr = " open" if single_scope or block.get("scope") is None else ""
    return (f'<details class="scope" id="{escape(block["anchor"])}"{open_attr}>'
            f'<summary>{escape(block["label"])}<span class="scope-peek">{escape(summary)}</span></summary>{content}</details>')


def _scope_fact(topic: str, scope: str | None, scope_label: str, artifacts: list[dict], result: dict) -> str:
    frames = {item["relative"]: _read_table(item["path"]) for item in artifacts if item["suffix"] == ".csv"}
    jsons = {item["relative"]: _read_json_artifact(item["path"]) for item in artifacts if item["suffix"] == ".json"}
    available = {path: frame for path, frame in frames.items() if frame is not None}
    slug = scope_slug(scope) if scope is not None else None
    if topic == "population":
        overview = available.get("tables/input_overview.csv")
        labels = available.get("tables/reconstruction_label_summary.csv")
        parts = []
        if overview is not None and {"side", "n_units"}.issubset(overview):
            parts.append("、".join(f"{str(row.side).upper()} {int(row.n_units)} 单位" for row in overview.itertuples()))
        if labels is not None and "n_units" in labels:
            total = _complete_sum(labels, "n_units")
            parts.append(f"SVC 标签{f'覆盖 {int(total)} 单位' if total is not None else '单位总数未知'}/{len(labels)} 类")
        if parts:
            return "；".join(parts) + "。两侧对象范围独立，数量差不表示单位丢失。"
    if topic == "rawbaseline":
        if slug:
            partition = available.get(f"tables/partitions_raw_{slug}.csv")
            if partition is not None and "partition" in partition:
                return f"{scope_label} 的 Raw baseline 使用 {len(partition)} 个实际入组单位，得到 {partition.partition.nunique(dropna=True)} 个 cluster。"
        partitions = available.get("tables/partition_summary.csv")
        if partitions is not None and {"side", "scope", "n_clusters"}.issubset(partitions):
            return "；".join(f"{str(row.side).upper()} / {row.scope}：{row.n_clusters} 类" for row in partitions.itertuples()) + "。各 scope 单独计数。"
        level2 = available.get("tables/raw_level2_baseline_summary.csv")
        if level2 is not None and "n_units" in level2:
            total = _complete_sum(level2, "n_units")
            return (f"Raw Level2 baseline 保存 {len(level2)} 类，共 {int(total)} 个单位。" if total is not None
                    else f"Raw Level2 baseline 保存 {len(level2)} 类，单位总数未知。")
    if topic == "membership" and slug:
        payload = jsons.get(f"tables/membership_summary_{slug}.json")
        summary = payload.get("summary") if isinstance(payload, dict) else None
        row = summary[0] if isinstance(summary, list) and summary else None
        if isinstance(row, dict) and pd.notna(row.get("n_units")):
            fraction = row.get("st_unit_change_fraction")
            tail = f"，匹配后变化 {float(fraction):.1%}" if pd.notna(fraction) else ""
            return f"{scope_label} 在共同 ID cohort 比较 {int(row['n_units'])} 个单位{tail}。"
    if topic == "genes" and slug:
        shared = available.get(f"tables/moran_shared_genes_{slug}.csv")
        if shared is not None and "comparison_available" in shared:
            usable = shared.loc[_bool_series(shared.comparison_available)]
            deltas = pd.to_numeric(usable.get("delta_moran", pd.Series(dtype=float)), errors="coerce").dropna()
            raw_n = _first_number(shared, "raw_n_units")
            svc_n = _first_number(shared, "svc_n_units")
            delta_text = f"，SVC−Raw 中位数 {deltas.median():.3g}" if not deltas.empty else ""
            return f"{scope_label} 共同可计算 Moran 基因 {len(usable)} 个（Raw {raw_n or 'NA'}、SVC {svc_n or 'NA'} 单位）{delta_text}。"
        native = []
        for side in ("raw", "svc"):
            frame = available.get(f"tables/moran_{side}_{slug}.csv")
            if frame is not None and "moran" in frame:
                native.append(f"{side.upper()}：{pd.to_numeric(frame.moran, errors='coerce').notna().sum()}/{len(frame)} 个基因可计算，{_first_number(frame, 'n_units')} 个单位")
        if native:
            return "；".join(native) + "。原生结果可单侧阅读，尚不足以完成两侧对照。"
    if topic == "program-overall":
        score_frames = [frame for path, frame in available.items() if path.startswith("tables/pathway_")]
        if score_frames:
            pieces = [f"{str(frame.side.iloc[0]).upper()} / {frame.program.iloc[0] if 'program' in frame else 'program'}：{frame.score.notna().sum()}/{len(frame)} 个评分 cohort 单位有分数"
                      for frame in score_frames if not frame.empty and {"side", "score"}.issubset(frame)]
            for payload in jsons.values():
                if isinstance(payload, dict):
                    coverage = payload.get("coverage", {}).get("fraction_present")
                    if pd.notna(coverage):
                        pieces.append(f"基因集覆盖率 {float(coverage):.1%}")
            if pieces:
                return "Program 总体评分：" + "；".join(pieces) + "。"
    if topic == "anatomy":
        windows = available.get("tables/raw_anatomy_windows.csv")
        if windows is not None and "level1_region" in windows:
            counts = windows.level1_region.value_counts()
            return f"完整 Raw Anatomy 共 {len(windows)} 个窗口：" + "、".join(f"{key} {value}" for key, value in counts.items()) + "。Interface 按窗口共存定义。"
        candidates = [frame for frame in available.values() if frame is not None and "n_units" in frame]
        if candidates:
            pieces = []
            for path, frame in available.items():
                if "n_units" not in frame:
                    continue
                total, known = _complete_sum(frame, "n_units"), _complete_sum(frame, "n_anatomy_covered")
                if total is not None:
                    side = str(frame.side.iloc[0]).upper() if "side" in frame and not frame.empty else scope_label
                    cohort = " baseline cohort" if "population" in frame else ""
                    pieces.append(f"{side}{cohort} {int(total)} 单位，Anatomy 覆盖 {int(known) if known is not None else '未知'}")
            if pieces:
                return "；".join(pieces) + "。各侧按实际观测点汇总。"
    if topic == "scale" and slug:
        support = available.get(f"tables/window_scale_support_{slug}.csv")
        configured = result.get("parameters", {}).get("parent_window_side_microns", {})
        configured = configured.get(scope) if isinstance(configured, dict) else configured
        if support is not None:
            recommendation = jsons.get(f"tables/window_scale_recommendation_{slug}.json") or {}
            recommended = recommendation.get("window_side_microns") if recommendation.get("status") == "ok" else None
            chosen = f"显式使用 {configured:g} μm" if configured is not None else "显式尺度未记录"
            suggestion = f"支持曲线推荐 {recommended:g} μm" if recommended is not None else "支持曲线无明确推荐"
            return f"{scope_label} {chosen}；{suggestion}。推荐不自动覆盖配置。"
    if topic == "diversity" and slug:
        frame = available.get(f"tables/window_diversity_svc_{slug}.csv")
        if frame is not None:
            valid = frame.loc[_bool_series(frame.valid_window)] if "valid_window" in frame else frame
            neff = pd.to_numeric(valid.get("neff", pd.Series(dtype=float)), errors="coerce").dropna()
            return f"{scope_label} 的 SVC 局部多样性有 {len(valid)} 个有效窗口；Neff 中位数为 {neff.median():.3g}。" if not neff.empty else f"{scope_label} 没有可读取的有效 Neff 窗口。"
    if topic in {"state", "gain"} and slug:
        extent = available.get(f"tables/{topic}_{slug}_region_extent.csv")
        if extent is not None and not extent.empty:
            row = extent.iloc[0]
            n_valid, n_region = row.get("n_valid_windows"), row.get("n_region_windows")
            threshold = jsons.get(f"tables/{topic}_{slug}_threshold.json") or {}
            status = threshold.get("status", "unknown")
            if pd.notna(n_valid) and pd.notna(n_region):
                return f"{scope_label} 的 {topic.title()} 有 {int(n_valid)} 个有效窗口、{int(n_region)} 个区域窗口；阈值状态 {status}。"
        windows = available.get(f"tables/{topic}_{slug}_region_windows.csv")
        if windows is not None:
            valid = int(_bool_series(windows.valid_window).sum()) if "valid_window" in windows else len(windows)
            known = int(windows.in_region.notna().sum()) if "in_region" in windows else 0
            return f"{scope_label} 的 {topic.title()} 连续场保存 {valid} 个有效窗口，其中 {known} 个有可判定的区域状态。"
    if topic == "region-anatomy" and slug:
        frame = available.get(f"tables/integrated_region_anatomy_summary_{slug}.csv")
        if frame is not None and {"region_kind", "in_region", "n_units", "denominator_units"}.issubset(frame):
            all_state = frame.loc[frame.region_kind.astype(str) == "state"]
            if all_state.empty or not all_state.in_region.notna().any():
                return f"{scope_label} 的 State 区域状态尚不可判定；Anatomy 记录仍保留，不能把未知写成零个区域单位。"
            state = all_state.loc[_bool_series(all_state.in_region)]
            n_value = _complete_sum(state, "n_units")
            n_value = 0 if state.empty else n_value
            denominators = pd.to_numeric(all_state.denominator_units, errors="coerce").dropna()
            if n_value is None or denominators.empty:
                return f"{scope_label} 的 Region × Anatomy 结果已保存，但单位分母未知。"
            n, denominator = int(n_value), int(denominators.max())
            unknown_rows = all_state.loc[all_state.in_region.isna()]
            unknown = 0 if unknown_rows.empty else _complete_sum(unknown_rows, "n_units")
            detail = "、".join(f"{row.anatomy_region} {int(row.n_units)}" for row in state.itertuples() if pd.notna(row.n_units))
            return f"{scope_label} 的 State 区域含 {n}/{denominator} 个 scope 内 SVC 单位（{detail or '无区域内单位'}）；区域状态未知 {int(unknown) if unknown is not None else '未记录'} 个。"
    if topic == "label-composition" and slug:
        frame = available.get(f"tables/integrated_label_composition_{slug}.csv")
        if frame is not None and not frame.empty and {"reconstruction_label", "n_units", "denominator_units"}.issubset(frame):
            counts = pd.to_numeric(frame.n_units, errors="coerce")
            if counts.notna().any():
                top = frame.loc[counts.idxmax()]
                state_value = top.get("state_region")
                region = "未知区域状态" if pd.isna(state_value) else ("State 内" if str(state_value).lower() in {"true", "1"} else "有效 State 外")
                denominator = int(top.denominator_units) if pd.notna(top.denominator_units) else "未知"
                return f"{scope_label} 共 {frame.reconstruction_label.nunique()} 类重建标签；单位数最多的组成项位于{region} / {top.get('anatomy_region', 'Unknown')}，标签 {top.reconstruction_label} 占该组 {int(top.n_units)}/{denominator}。"
    if topic == "program-location" and slug:
        summaries = [frame for path, frame in available.items() if path.endswith("_anatomy.csv")]
        if summaries:
            pieces = []
            for frame in summaries:
                scored, total = _complete_sum(frame, "n_scored"), _complete_sum(frame, "n_units")
                if scored is not None and total is not None:
                    side = str(frame.side.iloc[0]).upper() if "side" in frame else ""
                    program = str(frame.program.iloc[0]) if "program" in frame else "program"
                    pieces.append(f"{side} / {program} {int(scored)}/{int(total)} 个单位有分数")
            if pieces:
                return f"{scope_label} 空间汇总：" + "；".join(pieces) + "。各侧独立，未评分保留缺失。"
    if topic == "changed-location" and slug:
        frame = available.get(f"tables/integrated_changed_units_{slug}.csv")
        if frame is not None and {"n_compared", "n_changed"}.issubset(frame):
            compared_value, changed_value = _complete_sum(frame, "n_compared"), _complete_sum(frame, "n_changed")
            if compared_value is None or changed_value is None:
                return f"{scope_label} 的变化单位位置汇总已保存，但比较分母未知。"
            compared, changed = int(compared_value), int(changed_value)
            return f"{scope_label} 的位置汇总比较 {compared} 个共同 ID 单位，其中 {changed} 个 membership 发生变化。"
    if artifacts:
        return f"{scope_label} 保存了该主题证据，但缺少生成稳定摘要所需的字段；请展开核对表与限制。"
    return f"{scope_label} 未保存该主题产物；不可据此判断没有变化。"


def _denominator_note(table: dict | None) -> str:
    if table is None:
        return ""
    frame = _read_table(table["path"])
    if frame is None or frame.empty:
        return "主表为空或不可读取，不生成图形或无变化结论"
    denominator_columns = [column for column in ("denominator_units", "n_units", "n_compared", "n_scored",
                                                   "n_valid_windows", "n_common_genes", "n_computable_shared_genes")
                           if column in frame]
    parts = []
    for column in denominator_columns[:3]:
        values = pd.to_numeric(frame[column], errors="coerce").dropna()
        if not values.empty:
            unique = sorted(set(values.tolist()))
            shown = ", ".join(f"{value:g}" for value in unique[:5])
            parts.append(f"{column}={shown}" + ("…" if len(unique) > 5 else ""))
    return "；".join(parts)


def _impact_evidence(artifact: dict, *, primary: bool, used_images: set[str]) -> str:
    relative = escape(artifact["relative"])
    label = escape(artifact["label"])
    anchor = artifact_anchor(artifact["relative"])
    if primary:
        return f'<article><h4><a href="#{anchor}">{label}</a></h4><p class="muted">已在本 scope 的主证据位置展示。</p></article>'
    if artifact["suffix"] in {".png", ".jpg", ".jpeg", ".svg"}:
        if artifact["relative"] in used_images:
            return f'<article><h4><a href="#{anchor}">{label}</a></h4><p class="muted">图已在前文嵌入。</p></article>'
        used_images.add(artifact["relative"])
        return (f'<article id="{anchor}"><h4>{label}</h4><figure><img src="{relative}" alt="{label}">'
                f'<figcaption><a href="{relative}">打开原图</a></figcaption></figure></article>')
    return f'<article id="{anchor}"><h4>{label}</h4><p><a href="{relative}">打开完整保存产物</a></p></article>'


def _topic_unavailable(topic: str, scope: str | None, result: dict) -> list[str]:
    prefixes = {
        "population": ("svc_reconstruction", "input"), "rawbaseline": ("partition", "raw_level2"),
        "membership": ("membership",), "genes": ("moran",), "program-overall": ("pathway",),
        "anatomy": ("anatomy", "spatial_geometry", "parent_anatomy"), "scale": ("window_support",),
        "diversity": ("window_diversity",), "state": ("state",), "gain": ("gain", "raw_k_control"),
        "region-anatomy": ("integrated_spatial_evidence", "anatomy"),
        "label-composition": ("integrated_spatial_evidence",), "program-location": ("program_aggregation", "pathway"),
        "changed-location": ("membership",),
    }.get(topic, ())
    messages = []
    for item in _unavailable_items(result):
        component = item["component"]
        matched_prefix = next((prefix for prefix in prefixes
                               if component == prefix or component.startswith(prefix + ":")
                               or component.startswith(prefix + "_")), None)
        if matched_prefix is None:
            continue
        tokens = component.split(":")[1:]
        component_scopes = [candidate for candidate in configured_scopes(result)
                            if candidate in tokens or component == f"{matched_prefix}_{scope_slug(candidate)}"]
        if scope is not None and component_scopes and scope not in component_scopes:
            continue
        messages.append(f"{component}: {item['reason']}")
    return messages


def _legacy_metadata(result: dict) -> str:
    declared = result.get("sections")
    if not isinstance(declared, list):
        return ""
    rows = []
    for section in declared:
        if not isinstance(section, dict):
            continue
        title = str(section.get("title") or section.get("heading") or section.get("id") or "section")
        description = str(section.get("description") or "")
        refs = section.get("artifacts", [])
        if isinstance(refs, str):
            refs = [refs]
        refs_text = ", ".join(str(ref) for ref in refs) if isinstance(refs, list) else ""
        rows.append(
            f'<li><strong>{escape(title)}</strong>{": " + escape(description) if description else ""}'
            f'{"（" + escape(refs_text) + "）" if refs_text else ""}</li>'
        )
    return f'<details><summary>旧版显式 sections 元数据</summary><ul>{"".join(rows)}</ul></details>' if rows else ""


def _legacy_aliases(topic: str, result: dict) -> str:
    """Keep old explicit section hashes valid at their fixed-tree destination."""
    destinations = {
        "population": ("input",), "rawbaseline": ("baseline",), "scale": ("support",),
        "diversity": ("diversity",), "state": ("regions",), "anatomy": ("anatomy",),
        "genes": ("molecular",), "membership": ("membership",), "region-anatomy": ("integration",),
    }
    declared = result.get("sections")
    if not isinstance(declared, list):
        return ""
    ids = {str(section.get("id")) for section in declared if isinstance(section, dict) and section.get("id")}
    return "".join(
        f'<span id="section-{escape(section_id)}" class="legacy-anchor" aria-hidden="true"></span>'
        for section_id in destinations.get(topic, ()) if section_id in ids
    )


def _read_json_artifact(path: Path) -> object | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return None


def _bool_series(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values.dtype):
        return values.fillna(False)
    return values.astype(str).str.lower().isin({"true", "1", "yes"})


def _first_number(frame: pd.DataFrame, column: str) -> int | float | None:
    if column not in frame:
        return None
    values = pd.to_numeric(frame[column], errors="coerce").dropna()
    return values.iloc[0] if not values.empty else None


def _max_number(frame: pd.DataFrame, column: str) -> int | float | None:
    if column not in frame:
        return None
    values = pd.to_numeric(frame[column], errors="coerce").dropna()
    return values.max() if not values.empty else None


def _complete_sum(frame: pd.DataFrame, column: str) -> int | float | None:
    if column not in frame or frame.empty:
        return None
    values = pd.to_numeric(frame[column], errors="coerce")
    return values.sum() if values.notna().all() else None


def _declared_sections(result: dict, artifacts: list[dict]) -> list[dict]:
    """Resolve explicit section artifact references while preserving order."""
    declared = result.get("sections")
    if isinstance(declared, dict):
        declared = [
            dict(value, id=key) if isinstance(value, dict) else {"id": key}
            for key, value in declared.items()
        ]
    elif not isinstance(declared, list):
        declared = None
    by_ref = {}
    for artifact in artifacts:
        by_ref.setdefault(artifact["label"], artifact)
        by_ref.setdefault(artifact["relative"], artifact)
    if declared is None:
        return _legacy_sections(artifacts)

    sections = []
    for index, raw in enumerate(declared):
        if not isinstance(raw, dict):
            continue
        section_id = str(raw.get("id") or raw.get("key") or f"section-{index + 1}")
        title = str(raw.get("title") or raw.get("heading") or section_id.replace("_", " ").title())
        description = str(raw.get("description", ""))
        title, description = _localize_section(section_id, title, description)
        refs = raw.get("artifacts", [])
        if isinstance(refs, str):
            refs = [refs]
        selected, missing = [], []
        for ref in refs if isinstance(refs, list) else []:
            if isinstance(ref, dict):
                ref = ref.get("path") or ref.get("relative") or ref.get("artifact") or ref.get("key")
            artifact = by_ref.get(str(ref)) if ref is not None else None
            if artifact is None:
                if ref is not None:
                    missing.append(str(ref))
            elif artifact not in selected:
                selected.append(artifact)
        sections.append({
            "id": section_id,
            "title": title,
            "description": description,
            "artifacts": selected,
            "missing": missing,
        })
    return sections


def _legacy_sections(artifacts: list[dict]) -> list[dict]:
    sections = []
    for key, heading in _DEFAULT_GROUPS:
        selected = [artifact for artifact in artifacts if _group(artifact["relative"]) == key]
        if selected:
            sections.append({"id": key, "title": heading, "description": "", "artifacts": selected, "missing": []})
    return sections


def _localize_section(section_id: str, title: str, description: str) -> tuple[str, str]:
    mapping = _SECTION_LOCALIZATION.get(section_id, {})
    localized_title = mapping.get(title)
    localized_description = mapping.get(description)
    if isinstance(localized_title, tuple):
        title = localized_title[0]
        if not description or description in mapping:
            description = localized_title[1]
    elif isinstance(localized_title, str):
        title = localized_title
    if isinstance(localized_description, str):
        description = localized_description
    return title, description


def _render_section(section: dict) -> str:
    artifacts = section.get("artifacts", [])
    figures = [artifact for artifact in artifacts if artifact["suffix"] in {".png", ".jpg", ".jpeg", ".svg"}]
    evidence = "".join(_evidence(artifact) for artifact in artifacts if artifact not in figures)
    missing = section.get("missing", [])
    missing_html = "".join(f"<li>{escape(item)}</li>" for item in missing)
    if missing_html:
        evidence += f'<article><h3>声明的产物不可用</h3><ul>{missing_html}</ul></article>'
    if not figures and not evidence:
        evidence = '<p class="muted">本节没有可用的已保存产物。</p>'
    figure_html = "".join(
        f'<figure><img src="{escape(item["relative"])}" alt="{escape(item["label"])}">'
        f'<figcaption><a href="{escape(item["relative"])}">{escape(item["label"])}</a></figcaption></figure>'
        for item in figures
    )
    description = (
        f'<p class="section-description">{escape(section["description"])}</p>'
        if section.get("description") else ""
    )
    return (
        f'<section id="section-{escape(section["id"])}"><h2>{escape(section["title"])}</h2>'
        f'{description}{figure_html}<details open><summary>已保存证据（Saved evidence）</summary>{evidence}</details></section>'
    )


def _evidence(artifact: dict) -> str:
    label, relative = escape(artifact["label"]), escape(artifact["relative"])
    body = artifact.get("preview", "") or f'<p><a href="{relative}">打开已保存产物</a></p>'
    return f'<article><h3><a href="{relative}">{label}</a></h3>{body}</article>'


def _file_index(artifacts: list[dict]) -> str:
    if not artifacts:
        return '<p class="muted">没有索引到已保存输出文件。</p>'
    rows = "".join(
        f'<tr><td>{escape(item["label"])}</td><td><a href="{escape(item["relative"])}">{escape(item["relative"])}</a></td>'
        f'<td>{escape(item["suffix"] or "file")}</td><td>{item["size"]}</td></tr>'
        for item in artifacts
    )
    return '<table class="file-index"><thead><tr><th>输出键</th><th>保存路径</th><th>类型</th><th>字节数</th></tr></thead><tbody>' + rows + '</tbody></table>'


def _status_table(rows: list[dict], empty: str) -> str:
    if not rows:
        return f'<p class="muted">{escape(empty)}</p>'
    body = "".join(
        f'<tr><td>{escape(row["stage"])}</td><td>{escape(row["status"])}</td>'
        f'<td>{escape(row["reason"])}</td><td>{escape(row["detail"])}</td></tr>'
        for row in rows
    )
    return '<table><thead><tr><th>阶段</th><th>状态</th><th>原因</th><th>详情</th></tr></thead><tbody>' + body + '</tbody></table>'


def _stage_rows(result: dict) -> list[dict]:
    rows = []
    stages = result.get("stages", [])
    if isinstance(stages, dict):
        stages = [
            dict(value, stage=key) if isinstance(value, dict) else {"stage": key, "status": value}
            for key, value in stages.items()
        ]
    if isinstance(stages, list):
        for item in stages:
            if not isinstance(item, dict):
                continue
            rows.append({
                "stage": _localize_stage(str(item.get("stage") or item.get("id") or item.get("name") or "stage")),
                "status": _STATUS_LABELS.get(str(item.get("status", "unknown")), str(item.get("status", "unknown"))),
                "reason": _localize_reason(str(item.get("reason") or item.get("error") or "")),
                "detail": _localize_reason(str(item.get("detail") or item.get("message") or "")),
            })
    errors = result.get("stage_errors", [])
    if isinstance(errors, dict):
        errors = [
            dict(value, stage=key) if isinstance(value, dict) else {"stage": key, "reason": value}
            for key, value in errors.items()
        ]
    if isinstance(errors, list):
        for item in errors:
            if not isinstance(item, dict):
                continue
            stage = _localize_stage(str(item.get("stage") or item.get("component") or "stage"))
            if not any(row["stage"] == stage and row["status"] == _STATUS_LABELS["error"] for row in rows):
                rows.append({
                    "stage": stage,
                    "status": _STATUS_LABELS["error"],
                    "reason": _localize_reason(str(item.get("reason") or item.get("message") or item.get("error") or "")),
                    "detail": _localize_reason(str(item.get("detail") or "")),
                })
    return rows


def _localize_stage(value: str) -> str:
    """Translate workflow stage names while retaining scope/program names."""
    parts = value.split(":", 1)
    head = parts[0]
    labels = {
        **_STAGE_LABELS,
        "partition": "分群",
        "moran": "Moran",
        "pathway": "AUCell",
        "raw_k_control": "Raw K-control",
        "state": "State",
        "gain": "Gain",
        "program_aggregation": "program 聚合",
        "anatomy_mapping": "Anatomy 映射",
        "window_support": "窗口支持",
        "window_diversity": "窗口多样性",
    }
    translated = labels.get(head, head)
    return translated if len(parts) == 1 else f"{translated}：{parts[1]}"


def _localize_reason(value: str) -> str:
    """Translate common scientific availability messages at the UI boundary."""
    if not value:
        return ""
    exact = {
        "OmicVerse unavailable": "AUCell provider 不可用（OmicVerse unavailable）",
        "provider unavailable": "provider 不可用（provider unavailable）",
        "Raw Leiden baseline unavailable": "Raw Leiden baseline 不可用",
        "both native Moran results required": "需要两侧原生 Moran 结果",
        "no common valid spatial windows": "没有共同有效空间窗口",
        "no common unit IDs": "没有共同单位 ID",
        "SVC reconstruction labels unavailable": "SVC 重建标签不可用",
        "Raw baseline partition unavailable": "Raw baseline 分群不可用",
        "nonempty expression cohort and gene axis required": "需要非空表达 cohort 和基因轴",
        "expression must be finite and nonnegative": "表达值必须有限且非负",
        "coordinates must be finite": "坐标必须为有限值",
        "no jointly computable shared genes": "没有可共同计算的共享基因",
        "stable region threshold unavailable: no_stable_threshold": "稳定区域阈值不可用：no_stable_threshold",
    }
    if value in exact:
        return exact[value]
    result = value
    result = result.replace("expression matrix identity is unknown", "表达矩阵身份未知")
    result = result.replace("expression transformation state is unknown", "表达变换状态未知")
    result = result.replace("Raw column", "Raw 列")
    result = result.replace(" is unavailable", " 不可用")
    result = result.replace(" contains missing values", " 含有缺失值")
    result = result.replace("SVC parent local diversity unavailable", "SVC parent 局部多样性不可用")
    result = result.replace("SVC parent local diversity", "SVC parent 局部多样性")
    result = result.replace("spatial geometry unavailable", "空间几何信息不可用")
    result = result.replace("broad-label column is unavailable", "broad 标签列不可用")
    result = result.replace("no valid reconstructed labels in SVC scope", "SVC scope 中没有有效重建标签")
    result = result.replace("stable region threshold unavailable", "稳定区域阈值不可用")
    result = result.replace("not computable", "不可计算")
    result = result.replace("provider unavailable", "provider 不可用")
    result = result.replace("skipped:", "已跳过：")
    # Keep traceback and unusual provider messages intact; they are evidence,
    # while the surrounding page still tells the reader what kind of record it is.
    return result


def _unavailable_items(result: dict) -> list[dict[str, str]]:
    items = []
    unavailable = result.get("unavailable", [])
    if isinstance(unavailable, dict):
        unavailable = [
            dict(value, component=key) if isinstance(value, dict) else {"component": key, "reason": value}
            for key, value in unavailable.items()
        ]
    if isinstance(unavailable, list):
        for item in unavailable:
            if isinstance(item, dict):
                items.append({
                    "component": str(item.get("component") or item.get("stage") or "component"),
                    "reason": _localize_reason(str(item.get("reason") or item.get("message") or "unavailable")),
                })
    return items


def _vocabulary(result: dict) -> str:
    """Keep reconstruction-impact reading terms visible at report entry."""
    if result.get("analysis") not in (None, "reconstruction_impact") and not result.get("sections"):
        return ""
    return (
        '<p class="boundary"><strong>阅读口径：</strong>State 使用声明的 SVC parent 分群中的 '
        '<em>重建状态标签（reconstructed-state labels）</em>。'
        '<strong>ΔNeff</strong> 展示为 Gain 候选（Gain candidate）：'
        'ΔNeff = N_eff(SVC parent) − N_eff(Raw Leiden)，计算范围是共同有效窗口；'
        '它是描述性证据，不能单独解释为生物学改善。</p>'
    )


def _group(relative: str) -> str:
    value = relative.lower()
    if "membership" in value or "changed_map" in value:
        return "membership"
    if any(term in value for term in ("moran", "pathway", "program")):
        return "genes"
    if any(term in value for term in ("anatomy", "window", "state", "gain")):
        return "spatial"
    return "overall"


def _observations(artifacts: list[dict]) -> list[str]:
    """State only values directly present in compact saved tables."""
    messages = []
    for artifact in artifacts:
        relative, path = artifact["relative"], artifact["path"]
        if relative.endswith("overview.csv"):
            table = _read_table(path)
            if table is not None and {"side", "n_units", "n_genes"}.issubset(table):
                detail = ", ".join(
                    f"{ {'raw': 'Raw', 'svc': 'SVC'}.get(str(row.side), row.side) }："
                    f"{int(row.n_units)} 个单位、{int(row.n_genes)} 个基因"
                    for row in table.itertuples(index=False)
                )
                messages.append(f"独立概览记录了 {detail}。")
        elif relative.endswith("partition_summary.csv"):
            table = _read_table(path)
            if table is not None and {"side", "scope", "n_clusters"}.issubset(table):
                messages.append("分群摘要记录了每个可用侧和 scope 的原生 Leiden cluster 数。")
        elif "region_extent" in relative:
            table = _read_table(path)
            if table is not None and not table.empty and {"n_region_windows", "n_valid_windows"}.issubset(table):
                row = table.iloc[0]
                messages.append(
                    f"{artifact['label']} 在 {int(row['n_valid_windows'])} 个有效窗口中选中了 "
                    f"{int(row['n_region_windows'])} 个窗口。"
                )
    return list(dict.fromkeys(messages))


def _read_table(path: Path) -> pd.DataFrame | None:
    try:
        return pd.read_csv(path)
    except (OSError, pd.errors.ParserError, UnicodeDecodeError):
        return None
