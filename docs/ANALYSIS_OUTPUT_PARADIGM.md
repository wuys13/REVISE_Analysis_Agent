# REVISE 分析成果组织范式与层级架构协议 (Analysis Output Paradigm)

## 1. 核心设计哲学：金字塔型渐进式披露 (Progressive Disclosure)

为支撑**多样本/批量化分析（Batch Analysis）**以及**下游 AI Agent 自主科学发现**，输出目录必须遵循**“由外向内、由宏观到细节、自顶向下递进”**的严格分层范式：

```
[顶层 / 批次级]  README.md, 批量索引, 跨样本比较
       │
[样本级入口]    report.html (可视化仪表盘), README.md (宏观结论), summary.json
       │
[分析槽位级]    01_生信分析 / 02_去掩盖刻画 / 03_生物学发现 (结构化结论与假设)
       │
[统计指标级]    tables/*.csv (统计汇总表、通路分值矩阵、假说证据矩阵)
       │
[底层细节级]    detailed_data/*.csv (全量差异基因、原始边缘权重、逐点物理噪音明细)
```

### 分层数据规则：
1. **越往内部、越往下层，内容越偏向细节数据（CSV）**：保留最全量、无过滤的矩阵与明细数据，供精细重打分或生信专家深钻。
2. **往上开始出现统计级别的 CSV 与文本总结**：对细节数据进行降维提炼（如全谱系 Hallmark 均值矩阵、配体受体互作显著对、结构化假说陈述书）。
3. **最外层为宏观说明文档、综合科学结论与可视化网页报告**：供人类审阅者与上层规划 Agent 一秒把握全局，不需在海量底层 CSV 中盲目翻找。

---

## 2. 标准输出目录层级结构 (Standard Hierarchy)

```text
output/
├── README.md                                  # 批量总览：全批次索引、项目元数据、范式说明
└── <Cancer_Type>/                             # 癌种 / Batch 维度 (例如 CRC, BRCA)
    ├── README.md                              # 批次级概览：样本列表、队列特征、批次统计
    └── <Sample_ID>/                           # 单样本根目录 (例如 P2_CRC_Xenium)
        ├── sample.yaml                        # 样本数据接入契约与元数据
        ├── report.html                        # ★ 样本综合可视化仪表盘 (宏观结论/生物学假设/去掩盖全景)
        ├── README.md                          # ★ 样本级宏观科学结论与导航指南
        ├── summary.json                       # ★ 样本级机器可读全景事实汇总字典
        │
        └── analysis/                          # 样本分析成果标准根目录
            ├── analysis.json                  # 槽位注册清单与元数据
            │
            ├── 01_bioinformatic_analysis/     # 【阶段一：平权生信分析】
            │   ├── summary.json               # 模块级统计元数据
            │   ├── FINDINGS.md                # 模块级文字总结
            │   ├── tables/                    # 统计级 CSV
            │   │   ├── pathway_scores_by_subtype.csv # 全量亚群 Hallmark 通路活性
            │   │   ├── colocalization_zscores.csv    # 全量亚群空间共定位 Z 矩阵
            │   │   └── significant_lr_pairs.csv      # 空间物理限制性配体受体通讯对
            │   └── detailed_data/             # 最底层细节 CSV
            │       └── subtype_degs_full.csv  # 24 亚群 One-vs-Rest 全量 120 DEGs 明细
            │
            ├── 02_reconstruction_impact/      # 【阶段二：重建前后变化与原始去掩盖刻画】
            │   ├── summary.json               # 模块级统计元数据
            │   ├── FINDINGS.md                # 模块级文字总结
            │   ├── tables/                    # 统计级 CSV
            │   │   ├── raw_vs_svc_unmasking_evidence.csv # 原始掩盖 vs 重建揭示定量对比 (20.7x 增益)
            │   │   ├── moran_gain_by_lineage.csv         # 谱系空间连续性增益 (+0.234)
            │   │   └── diversity_entropy_delta.csv       # 微环境邻域熵与保真度审计
            │   └── detailed_data/             # 最底层细节 CSV
            │       ├── contamination_clearance.csv       # MUC2/KRT8 假性交叉污染清除明细
            │       └── noise_reduction_diffusion.csv     # 横向物理扩散与光学 Dropout 恢复明细
            │
            └── 03_biological_discoveries/     # 【阶段三：综合分析下的生物学发现与科学假设】
                ├── summary.json               # 模块级统计元数据
                ├── FINDINGS.md                # 模块级文字总结
                ├── hypotheses/                # ★ 核心结构化假说陈述书 (供 Agent 自主推演)
                │   ├── HYPOTHESIS_1_Treg_Mreg_TLS_Niche.md
                │   ├── HYPOTHESIS_2_SPP1_TAM_MyoCAF_Invasive_Shell.md
                │   ├── HYPOTHESIS_3_Tprolif_Inflammatory_TAM_Core_Recruitment.md
                │   └── HYPOTHESIS_4_Contractile_vs_Proteolytic_CAF_Heterogeneity.md
                ├── tables/                    # 统计级 CSV
                │   ├── notable_subtypes_summary.csv      # 24 亚群多维度优先级评估大表 (Priority, SNR, Compactness)
                │   └── hypotheses_evidence_matrix.csv    # 假说量化证据支撑矩阵 (Z-score, CCI, Pathway)
                └── detailed_data/             # 最底层细节 CSV
                    └── target_lineage_programs.csv       # 核心靶标功能谱系特征程序基因集
```

---

## 3. 各层级内容职责与消费契约

| 层级 | 所在位置 | 包含内容形态 | 核心消费群体与用途 |
|---|---|---|---|
| **L0: 批次层** | `output/README.md`, `output/<Cancer>/README.md` | 队列说明、样本索引表格、全局组织架构 | 批量运行脚本、高层项目管理 Agent，快速定位感兴趣样本 |
| **L1: 样本层** | `output/<Cancer>/<Sample>/` (`report.html`, `README.md`, `summary.json`) | 宏观结论、生物学发现摘要、交互式网页仪表盘、全景事实 JSON | 人类PI审阅成果、下游自主科学探索 Agent 建立全局心智模型 |
| **L2: 槽位层** | `analysis/01_...`, `02_...`, `03_...` | 模块级 `FINDINGS.md`, `summary.json`, 结构化假说文档 | 垂直领域分析 Agent（如免疫微环境 Agent、基质硬化 Agent）获取子领域结论 |
| **L3: 统计层** | 各槽位下 `tables/*.csv` | 通路矩阵、共定位 Z 矩阵、多维排序大表、假说量化矩阵 | 统计建模、假设打分算法、交叉比对脚本 |
| **L4: 细节层** | 各槽位下 `detailed_data/*.csv` | 120+ 行全量差异基因、原始边缘、噪音污染明细 | 基因级深度检索、生信流水线二次过滤、可复现性审计 |

---

## 4. 批量运行扩展性机制 (Batch Extensibility)

当后续启动多样本（如 P1, P2, P3, P4...）批量分析时：
1. 每个样本按照本规范生成独立的 `output/<Cancer>/<Sample_ID>/` 树状目录；
2. 批量分析器自动读取每个样本的 `summary.json`，在 `output/<Cancer>/README.md` 中聚合生成**多样本交叉对比矩阵**与**跨样本假说泛化表**；
3. 任何外部消费者无需修改解析协议，即可实现从宏观批次到微观单细胞基因表达的平滑导航。
