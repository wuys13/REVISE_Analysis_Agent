# 02. 科学发现域规范（Discoveries Specification）

## 一、 领域定位与先验哲学

`01_discoveries/` 专职回答**“全转录组空间重建带来了哪些新的生物学见解与机制假说”**。

### 核心科学原则
1. **全基因组无偏分析为主**：
   原位 Xenium 仅覆盖几百个靶向基因，绝大多数未测标志物、配体-受体对与全套信号通路在 Raw 中缺失或残缺，**不强行进行与 Raw 的无意义差值对比**。所有分析以重建后的全转录组空间矩阵为准开展无偏挖掘。
2. **谱系范围控制开关（Lineage Scope Switch）**：
   - **默认聚焦模式（Default）**：`["Fibroblast", "Mono_Macro", "T"]`，聚焦肿瘤微环境（TME）最具代表性的三大免疫与基质轴线。
   - **全量无偏模式（All Lineages）**：一键切换覆盖单细胞参考集中的全部 14 个 Level 1 谱系（Tumor, SMC, Plasma, Intestinal Epithelial, B, Vascular EC, Pericyte 等）。
   - **自定义模式（Custom List）**：按研究课题按需指定谱系列表。
3. **未注释谱系的自适应 Leiden 聚类回退（Unsupervised Leiden Fallback）**：
   - 单细胞参考集内部分谱系（如 Tumor, SMC, Plasma, Intestinal Epithelial）缺乏人工细分 Level 2 标签。
   - **机制规范**：重建侧在单细胞端提前运行 Leiden 聚类，为缺乏 Level 2 的群体生成聚类标签并作为重建先验；分析侧自适应检测标签，若仅有粗分类则自动调用空间/转录本 Leiden 聚类簇或优先索引 `SVC_cluster`，保证全量展开时不遗漏任何谱系的内部异质性。
4. **双层粒度架构与黄金生物学基准参照系（Dual-Resolution & Biological Ground Truth）**：
   - **宏观对照层（Level 2 Macro Subtypes）**：提供跨谱系的一致性概览（19个类别）；
   - **微观高分辨率层（Fine-Grained Subclusters）**：以发表论文与复现案例中的黄金亚群作为绝对参照系：
     * `Treg (T_5)`：`FOXP3+`, `IL2RA+`, `CTLA4+`, `LTB+`，交界面富集，与 Mreg 共定位；
     * `Mreg (TAM_0)`：`MMP12+`, `ENPP2+`, `DNASE1L3+`，交界面富集，介导 ANXA1-FPR3 与 LTB-LTBR；
     * `SPP1+ TAM (TAM_4)`：`SPP1+`, `CHI3L1+`, `LPL+`，形成肿瘤外周屏障（Peritumoral Shell）；
     * `CAF 亚群分工`：`CAF_1` (ECM-remodeling, COMP+), `CAF_4` (Hypoxia/LPL), `CAF_5` (SMC-derived, MYH11+), `CAF_8` (ECM-breakdown, MMP1+)；
   - **严禁纸面评估，以真实实验结果作为检验准绳**：若未检出黄金亚群，即刻反思是否分析粒度过粗或存在过度平滑；若产生新亚群（如 Tumor 侵袭簇），结合空间前沿与互作特征科学评估其真实性。

---

## 二、 发现域子槽位计算契约与底表规范

### 1. `01_differential_expression/`（全亚型无偏差异表达与 Marker 挖掘）
* **计算范围**：
  - 采用双层粒度计算：同时输出宏观 Level 2 与微观高分辨率子群（包含无监督 Leiden 簇）的差异基因；
  - 覆盖全部 13,000+ 重建基因，特别标记出在原位 Panel 之外新恢复的基因（`is_expanded_gene == True`）。
* **自动化筛选评分公式**：
  $$\text{Priority Score} = |\log_2 \text{FC}| \times (-\log_{10} p_{\text{adj}}) \times \text{AUC}_{\text{specificity}}$$
* **产物契约**：
  - `tables/notable_subtypes_summary.csv`：宏观 Level 2 亚型优先级汇总；
  - `tables/benchmark_subclusters_gold_standard.csv`：微观高分辨率黄金亚群与 Leiden 簇表；
  - `tables/raw_vs_svc_subcluster_degs_comparison.csv`：Raw vs SVC 在亚群差异基因检出数（18.6x ~ 23.9x 放大）与特异性提升上的定量对比表；
  - `tables/subtype_degs_full.csv`：全量差异基因明细。

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
