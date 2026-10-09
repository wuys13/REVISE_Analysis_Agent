# 02. 科学发现域规范（Discoveries Specification）

## 一、 领域定位与先验哲学

`01_discoveries/` 专职回答**“全转录组空间重建带来了哪些新的生物学见解与机制假说”**。

### 核心科学原则
1. **全基因组无偏分析为主**：
   原位 Xenium 仅覆盖几百个靶向基因，绝大多数未测标志物、配体-受体对与全套信号通路在 Raw 中缺失或残缺，**不强行进行与 Raw 的无意义差值对比**。所有分析以重建后的全转录组空间矩阵为准开展无偏挖掘。
2. **亚群先验与无指向性分析的辩证结合**：
   - **锁定三大核心主群（Lineages）**：研究焦点明确锁定在微环境三大核心谱系：`Fibroblast`（成纤维细胞）、`Mono_Macro`（单核/巨噬细胞）、`T`（T 细胞）。
   - **细分亚型（sc-SVC subtypes）无指向性全量计算**：对三大谱系内部的所有 sc-SVC 细分亚型，**不做预设偏好的过滤，统一、完整地全部计算一遍**（全部跑差异基因、通路活性、细胞互作与共定位）。
   - **统计指标驱动自动筛选（Metric-Driven Prioritization）**：通过严密的统计量（特异性打分、空间方差、受配体相互作用显著性），自动筛选出“最值得深入探究的亚型与生物学实体”，辅助 Agent 进行自主发现。

---

## 二、 发现域子槽位计算契约与底表规范

### 1. `01_differential_expression/`（全亚型无偏差异基因与特异 Marker）
* **计算范围**：
  - 对锁定谱系下的所有 sc-SVC Level2 细分亚型进行一对多（One-vs-Rest）Wilcoxon 差异表达检验；
  - 覆盖全部 13,000+ 重建基因，特别标记出在原位 Panel 之外新恢复的基因（`is_expanded_gene == True`）。
* **自动化筛选评分公式**：
  $$\text{Priority Score} = |\log_2 \text{FC}| \times (-\log_{10} p_{\text{adj}}) \times \text{AUC}_{\text{specificity}}$$
* **产物契约**：
  - `summary.json`：
    ```json
    {
      "lineages_analyzed": ["Fibroblast", "Mono_Macro", "T"],
      "subtypes_analyzed_count": 12,
      "total_significant_degs": 1420,
      "top_priority_subtypes": ["Fibroblast_c2", "Mono_Macro_c1"],
      "top_novel_expanded_markers": {
        "Fibroblast_c2": ["POSTN", "COL11A1"],
        "Mono_Macro_c1": ["SPP1", "APOC1"]
      }
    }
    ```
  - `tables/subtype_degs_full.csv`：
    表头包含 `lineage`, `subtype`, `gene_symbol`, `log2FC`, `p_val_adj`, `auc`, `pct_in`, `pct_out`, `is_expanded_gene`, `priority_score`。
  - `tables/notable_subtypes_summary.csv`：按 `priority_score` 汇总出的最显著亚型及其 Top Marker 列表。

---

### 2. `02_pathway_activity/`（全谱信号通路与空间梯度）
* **计算范围**：
  - 基于 MSigDB Hallmark 50 个完整基因集，使用 AUCell 计算逐细胞表型程序活性；
  - 在组织空间网格上计算各通路的空间自相关（Moran's I）与空间变化方差。
* **自动化筛选依据**：
  筛选空间方差与空间自相关均显著偏高的通路，界定为“具有组织极性与浸润梯度的关键通路”（如 EMT、乏氧、糖酵解、干扰素响应）。
* **产物契约**：
  - `summary.json`：
    ```json
    {
      "pathways_evaluated": 50,
      "spatially_active_pathways_count": 8,
      "top_gradient_pathways": [
        {"name": "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "moran_i": 0.52, "gradient_strength": "high"},
        {"name": "HALLMARK_HYPOXIA", "moran_i": 0.44, "gradient_strength": "moderate"}
      ]
    }
    ```
  - `tables/pathway_scores_by_subtype.csv`：各细分亚型在 50 条 Hallmark 上的均值与 Z-score 矩阵。
  - `tables/spatial_gradient_pathways.csv`：各通路的空间自相关、方差及热点富集坐标统计表。

---

### 3. `03_cell_communication/`（空间物理约束下的受配体互作 CCI）
* **计算范围**：
  - 依赖细胞空间 Delaunay / kNN 临近图（设定物理距离截断，如 50μm）；
  - 全基因组配体-受体数据库（如 CellChat / LIANA 核心对）；
  - 在邻近的细分亚型配对之间，计算受配体协同表达置信度与置换检验 $p$ 值。
* **自动化筛选评分**：
  提取具有最高互作特异性且跨越不同谱系（如成纤维细胞 $\rightarrow$ T 细胞、巨噬细胞 $\rightarrow$ 肿瘤）的受配体轴。
* **产物契约**：
  - `summary.json`：
    ```json
    {
      "distance_threshold_um": 50,
      "tested_pairs_count": 2150,
      "significant_lr_pairs_count": 86,
      "dominant_cross_lineage_axes": [
        {"sender": "Fibroblast_c2", "receiver": "T_c1", "lr_pair": "CXCL12_CXCR4", "score": 3.82},
        {"sender": "Mono_Macro_c1", "receiver": "Fibroblast_c2", "lr_pair": "TGFB1_TGFBR2", "score": 3.45}
      ]
    }
    ```
  - `tables/significant_lr_pairs.csv`：
    表头包含 `sender_subtype`, `receiver_subtype`, `ligand`, `receptor`, `interaction_score`, `p_val`, `physical_edge_support`。

---

### 4. `04_spatial_colocalization/`（空间生态位 Niche 与共定位）
* **计算范围**：
  - 计算全量细分亚型间的空间邻域富集（Neighborhood Enrichment Z-score）；
  - 基于细胞局部组成聚类划分组织空间生态位（Spatial Niches，如肿瘤核心区、浸润前沿交界面、三级淋巴结构 TLS 样区）。
* **产物契约**：
  - `summary.json`：
    ```json
    {
      "niches_detected_count": 4,
      "dominant_colocalization_patterns": [
        "Fibroblast_c2 colocalizes with Mono_Macro_c1 (Z=4.12)",
        "T_c1 excludes Fibroblast_c4 (Z=-3.25)"
      ]
    }
    ```
  - `tables/neighborhood_enrichment_matrix.csv`：所有亚型间的共生/排斥 Z-score 矩阵。
  - `tables/spatial_niche_composition.csv`：各个空间生态位的细胞亚型比例与解剖区域对应表。
