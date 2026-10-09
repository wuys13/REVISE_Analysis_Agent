# 00. 生信分析槽位总览与渐进式披露范式规范

## 一、 规范背景与核心目标

在 REVISE 空间转录组重建生态中，下游 `analysis/` 槽位不仅是生信算法运行产物的存放位置，更是**供 AI Agent 后续自主进行科学发现与证据检索的标准化知识库**。

为了避免“单样本巨量数据导致 Agent 上下文溢出”、“各分析模块命名与格式随意混乱”以及“强行跨模态不平等对比导致的统计失真”，本规范确立两大核心基石：
1. **双域架构（Dual-Domain Architecture）**：将“方法保真度评估（Benchmarks）”与“全转录组生物学挖掘（Discoveries）”在物理目录上彻底解耦。
2. **渐进式披露（Progressive Disclosure）**：以标准 4 级梯次交付结构，让 Agent 在数 KB 的摘要与底表中自顶向下完成科学假设提出、证据检验与机制发现。

---

## 二、 顶级双域目录组织架构

每个样本的分析输出槽位（`output/<group>/<sample>/analysis/`）必须严格按照以下结构组织：

```text
output/<group>/<sample>/analysis/
├── README.md                              # [Level 0 导航] 全局生信分析索引与 Agent 导航指南
├── analysis.json                          # [Level 0 机器索引] 全槽位运行状态、参数与核心指标总览
│
├── 00_benchmarks/                         # =====【域 A：保真度与方法对比】=====
│   ├── README.md                          # 评测标准、基因对齐协议与显著性阈值说明
│   ├── summary.json                       # 核心质控指标（Moran 提升率、白噪声残差、保真度打分）
│   ├── FINDINGS.md                        # 面向 Agent 的评测裁定（结论：重建质量是否达标）
│   ├── 01_spatial_continuity/             # 空间连续性与白噪声消除对比
│   └── 02_microenvironment_fidelity/      # 微环境局部物理状态保真度与 Gain 检验
│
└── 01_discoveries/                        # =====【域 B：全转录组自主科学发现】=====
    ├── README.md                          # 发现域图谱（全细胞、全空间、全基因）
    ├── summary.json                       # 科学发现全景摘要（Top 新 Marker、最活跃通路、主导互作轴）
    ├── FINDINGS.md                        # 科学发现陈述（生物学见解与机制推断）
    ├── 01_differential_expression/        # 全局无偏差异基因与特异 Marker 挖掘
    ├── 02_pathway_activity/               # 全谱信号通路与表型程序（AUCell / GSEA）
    ├── 03_cell_communication/             # 空间物理约束下的受配体网络 (CCI)
    └── 04_spatial_colocalization/         # 空间生态位 (Niche) 与细胞共定位关系
```

---

## 三、 四级渐进式披露契约（Progressive Disclosure Contract）

为保证 Agent 在消耗极低 Token 的前提下拥有完整推理链条，**每个具体的分析子槽位**（如 `01_spatial_continuity/` 或 `01_differential_expression/`）必须统一包含以下 4 级契约文件：

| 级别 | 文件对象 | 格式与大小 | 面向对象与核心职责 |
| :--- | :--- | :--- | :--- |
| **Level 0** | `README.md` | Markdown (<5 KB) | **概念说明**：该槽位的科学假说、依赖输入前提、算法公式、参数网格与质控阈值。 |
| **Level 1 (机器)** | `summary.json` | JSON (<5 KB) | **机器速读**：核心标量指标（如中位数差异、显著性 $p$ 值、通过率）、Top-5 关键生物学实体。供 Agent 代码做条件分支判定。 |
| **Level 1 (人类/语义)** | `FINDINGS.md` | Markdown (<10 KB) | **科学裁定**：结构化自然语言结论。包含“已确认事实”、“机制推断”、“已知局限性”与“建议后续探索方向”。 |
| **Level 2** | `tables/*.csv` | Tidy CSV (100 KB ~ 5 MB) | **统计底表**：标准二维数据表（强制包含标准列头），供 Agent 编写轻量 Pandas/SQL 脚本靶向过滤与排名。 |
| **Level 3** | `figures/` & `report.html` | PNG/SVG & HTML | **深层证据与交互**：多色空间切片图、火山图、小提琴图及沉浸式 HTML 报告，供深层证据追溯。 |

---

## 四、 顶层索引规范：`analysis.json` 与 `README.md`

### 1. `analysis.json` 机器格式
位于 `analysis/analysis.json`，是 Agent 首次接触样本分析时的**单点总控入口**：
```json
{
  "schema_version": 2,
  "sample_id": "P2_CRC_Xenium",
  "generated_at": "2026-10-09T22:00:00Z",
  "domains": {
    "benchmarks": {
      "status": "succeeded",
      "verdict": "high_fidelity_smoothness_validated",
      "summary_path": "00_benchmarks/summary.json",
      "aspects": ["01_spatial_continuity", "02_microenvironment_fidelity"]
    },
    "discoveries": {
      "status": "succeeded",
      "verdict": "rich_biology_recovered",
      "summary_path": "01_discoveries/summary.json",
      "aspects": [
        "01_differential_expression",
        "02_pathway_activity",
        "03_cell_communication",
        "04_spatial_colocalization"
      ]
    }
  }
}
```

### 2. `README.md` 人机自解释格式
清晰列出：
- 当前样本所属病人、癌种、测序平台（从 `sample.yaml` 继承）；
- 双域的核心结论摘要一览表；
- 建议 Agent 调用的优先顺序与快速过滤代码片段。
