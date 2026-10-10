# 01. 交付架构与目录形态规范 (Architecture & Directory Layout)

> **上级入口**：[返回交付框架总览 (README.md)](README.md)  
> **关键关联**：[05. Agent 数据契约与载荷规范](05_agent_data_contract_and_payload.md) ｜ [03. 重建五层系统解构矩阵](03_five_layer_reconstruction_impact.md)

---

## 1. 倒金字塔渐进式披露理念 (The Inverted Pyramid Paradigm)

在处理超分辨率空间转录组（如 50,000+ 细胞、13,000+ 基因）数据时，直接将数吉字节（GB）的 H5AD 表达矩阵输入下游大模型（LLM / Agent）会导致极度严重的 **Token 爆炸** 与 **幻觉漂移**。

为此，REVISE 采用**倒金字塔渐进式披露架构（Inverted Pyramid Progressive Disclosure）**：
* **越往顶层（L1 顶层入口）**：信息密度越高度提炼，突出宏观科学发现、多维交互大屏与标准化机器载荷；
* **越往深层（L2 业务槽位 $\rightarrow$ L3 统计表 $\rightarrow$ L4 全量矩阵）**：数据粒度越精细，保留完整的统计学分布与基因级表达事实。

```
[Level 1: 顶层宏观大屏与机器载荷]  report.html (交互式仪表盘), summary.json, agent_input_payload.json
             │
[Level 2: 业务槽位与解构报告]    analysis/
             ├── 01_bioinformatic_analysis/   (平权先验生信计算发现)
             ├── 02_reconstruction_impact/    (五层系统去掩盖刻画)
             └── 03_biological_discoveries/   (结构化科学假设 H1–H4)
             │
[Level 3: 统计指标汇总底表]      各槽位 tables/*.csv (135 DEGs, 40 ORA, 14 CPDB)
             │
[Level 4: 细粒度全量数据矩阵]    各槽位 detailed_data/*.csv (全谱系基因级底表)
```

---

## 2. 单样本标准交付目录拓扑树 (Sample Delivery Tree)

以结直肠癌原发灶 Xenium 重建样本 `P2_CRC_Xenium` 为例，标准化交付目录形态如下：

```text
output/CRC/P2_CRC_Xenium/
├── report.html                           # [L1] 宏观综合交互式科学大屏 (人类 PI / 学术审阅入口)
├── README.md                             # [L1] 样本宏观结论、假设摘要与导航索引
├── summary.json                          # [L1] 样本全局统一事实 JSON (33 字段核心概览)
├── agent_input_payload.json              # [L1] 面向下游 Agent 的紧凑结构化载荷 (43KB，防幻觉事实源)
└── analysis/                             # [L2] 渐进式下钻分析槽位根目录
    ├── analysis.json                     # 槽位清单清单与架构定义
    ├── README.md                         # 槽位调用指引
    │
    ├── 01_bioinformatic_analysis/        # 槽位 1：全谱系平权生信分析
    │   ├── summary.json, FINDINGS.md     # 槽位结论摘要
    │   ├── tables/                       # [L3] 统计汇总表
    │   │   ├── pathway_scores_by_subtype.csv (40 行 MSigDB Hallmark ORA)
    │   │   ├── significant_lr_pairs.csv      (14 行 CellPhoneDB 显著受配体对)
    │   │   └── colocalization_zscores.csv    (描述性空间共定位峰值记录)
    │   └── detailed_data/                # [L4] 细粒度基因矩阵
    │       ├── focal_subtypes_top_degs.csv   (9 焦点亚群 × 15 top DEGs = 135 行)
    │       └── subtype_degs_full.csv         (全量亚群差异基因底表)
    │
    ├── 02_reconstruction_impact/         # 槽位 2：重构赋能与去掩盖刻画
    │   ├── summary.json, FINDINGS.md     # 槽位结论摘要
    │   ├── tables/
    │   │   ├── raw_vs_svc_unmasking_evidence.csv (关键表型 DEG 增幅对比，Tprolif 达 18.6x)
    │   │   └── systemic_impact_5layers.csv       (系统化五层覆盖矩阵：Unit/Gene/Path/CCI/Niche)
    │   └── legacy_archive/               # 早期模拟指标归档隔离区 (退出活跃事实层)
    │
    └── 03_biological_discoveries/        # 槽位 3：生物学发现与假设推演
        ├── summary.json, FINDINGS.md     # 槽位结论摘要
        ├── tables/
        │   ├── notable_subtypes_summary.csv      (8 大焦点亚群全特征宽表)
        │   └── hypotheses_evidence_matrix.csv    (4 大科学假设多维证据矩阵)
        └── hypotheses/                   # 结构化科学假设 Markdown 详细文档
            ├── HYPOTHESIS_1_Treg_Mreg_TLS_Niche.md
            ├── HYPOTHESIS_2_SPP1_TAM_MyoCAF_Invasive_Shell.md
            ├── HYPOTHESIS_3_Tprolif_Inflammatory_TAM_Core_Recruitment.md
            └── HYPOTHESIS_4_Contractile_vs_Proteolytic_CAF_Heterogeneity.md
```

---

## 3. 分层职责、数据规模与消费对象矩阵

| 层级 | 入口文件名 | 数据规模 | 核心职责 | 消费对象 | 容错与防幻觉设计 |
|---|---|:---:|---|---|---|
| **Level 1** | `report.html` | ~21 KB | 全局可视化交互、假设卡片、图表汇总 | 人类研究员、评审专家 | 仅展示已计算事实，支持交互检索 |
| **Level 1** | `summary.json` | ~5 KB | 快速元数据核验、样本亚群概览 | 轻量级自动化工作流 | 严格受限的标量与布尔值 |
| **Level 1** | `agent_input_payload.json` | ~43 KB | 机器消费的标准结构化证据全集 | 下游科学发现主 Agent | **限制在 40–50KB 最佳认知窗口内，消除伪推理** |
| **Level 2** | `01/02/03` 下的 `FINDINGS.md` | 1–3 KB/篇 | 业务模块叙述性发现、机理解读 | 领域专职 Sub-agent | 分槽位独立下钻，避免跨槽位上下文混淆 |
| **Level 3** | `tables/*.csv` | 5–15 KB/表 | 差异基因、富集通路、受配体互作统计量 | 生信核验与二次作图脚本 | 包含源文件行号、键名与真实检验统计量 |
| **Level 4** | `detailed_data/*.csv` | 50–500 KB | 全量基因 logFC、未截断的表达谱矩阵 | 深度挖掘与定制化分析 | 备用全量事实源，不默认全量载入上下文 |

---

## 4. 防幻觉物理隔离与安全性设计

1. **计算代码与展示数据的物理解耦**：
   原始庞大的 AnnData 矩阵（如 7.2GB 的 `SVC.h5ad`）保存在集群持久存储中；分析交付包仅抽取经过严格置换检验、多重假设检验校正（FDR/padj）的显著性表格。
2. **单一事实写入源 (Single Source of Truth)**：
   所有中间指标提取由标准提取器 [`prepare_agent_payload.py`](../../scripts/prepare_agent_payload.py) 统筹完成，严禁由大模型在生成 Markdown 时临时“心算”统计量。
3. **退役旧生成器加锁**：
   历史上的原型脚本（`organize_three_tier_analysis.py` 等）已被打上安全防御守卫（`sys.exit(1)`），杜绝并发误运行覆盖生产事实。
