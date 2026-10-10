# 02. 客观统计基线与焦点亚群排序规范 (Statistical Baseline & Focal Subtypes)

> **上级入口**：[返回交付框架总览 (README.md)](README.md)  
> **关键关联**：[03. 重建五层系统解构矩阵](03_five_layer_reconstruction_impact.md) ｜ [04. 生物学故事与文献对齐](04_biological_validation_and_stories.md)

---

## 1. 25 个重建亚群宇宙全景 (Subcluster Universe)

在结直肠癌原发灶 Xenium 数据（`P2_CRC_Xenium`，50,626 个高质量重构细胞）中，单细胞超分辨率重构涵盖三大核心微环境谱系，共解耦出 **25 个精细亚群**：

* **成纤维细胞大类 (Fibroblast Compartment)**：10 个亚群（`CAF_0` 至 `CAF_9`）
* **单核巨噬细胞大类 (Mono_Macro Compartment)**：8 个亚群（`TAM_0` 至 `TAM_7`）
* **T 淋巴细胞大类 (T Compartment)**：7 个亚群（`T_0` 至 `T_6`）

全量亚群差异基因底表保存在 [`detailed_data/subtype_degs_full.csv`](../../output/CRC/P2_CRC_Xenium/analysis/01_bioinformatic_analysis/detailed_data/subtype_degs_full.csv)，保证了全谱系平权先验生信基线的完整性。

---

## 2. 核心焦点亚群多维排序宽表 (Top Focal Subtypes Multi-metric Table)

为了在下游科学推理中重点攻坚核心微环境机制，系统基于真实生信特征（标记基因特异性 logFC、padj、微环境解剖位置与功能通路）聚焦出 8 大核心亚群。该排序与论文后验（Cell Fig 6/7）完全吻合，消除了任何人工公式构造的伪分数：

| 选择顺位 | 亚群 ID | 论文编号与表型 | 所属大类 | 代表性核心标记物 (Top DEGs) | 组织解剖生态位 (Spatial Niche) | 关键 MSigDB Hallmark ORA 通路 (Combined Score / padj) | 预后与生物学本质角色 |
|:---:|:---:|:---:|:---:|---|---|---|---|
| **#1** | **T_5** | T_3 (Treg) | T 淋巴 | `FOXP3, CCR4, TIGIT, RTKN2, CTLA4, LTB` | 肿瘤-正常交界区 / TLS 生态位 | **IL2_STAT5** (323.19, $p_{\text{adj}}=3.04\times 10^{-8}$)<br>**IL6_JAK_STAT3** (243.88, $p_{\text{adj}}=1.09\times 10^{-5}$) | **有利 DFS (HR < 1)**；TLS 情境依赖的保护性调节中枢，非单纯促癌抑制 |
| **#2** | **TAM_0** | TAM_3 (Mreg) | 单核巨噬 | `MMP12, RNF152, ENPP2, DNASE1L3` | 肿瘤-正常交界区 / TLS 生态位 | **KRAS_UP** (37.93, $p_{\text{adj}}=0.044$)<br>**INFLAMMATORY** (22.44, $p_{\text{adj}}=0.123$) | **有利 DFS (HR < 1)**；与 T_5 空间高亲和共存，ANXA1 自分泌维持慢性稳态 |
| **#3** | **TAM_4** | TAM_2 (SPP1+) | 单核巨噬 | `SPP1, CSTB, CHI3L1, FCRL2, LPL` | 肿瘤外周促纤维增生带 (Desmoplastic rim) | **COAGULATION** (148.24, $p_{\text{adj}}=1.36\times 10^{-5}$)<br>**ANGIOGENESIS** (119.11, $p_{\text{adj}}=0.00188$) | **不良预后**；包绕肿瘤外周侵袭边缘，协同 CAF_1 构建纤维化促血管屏障 |
| **#4** | **CAF_1** | CAF_1 (myoCAF) | 成纤维 | `GSTA1, COMP, COL11A1, TIMP3, COL1A1` | 肿瘤侵袭前沿基质 | **EMT** (4762.76, $p_{\text{adj}}=3.16\times 10^{-31}$)<br>**COAGULATION** (119.57, $p_{\text{adj}}=3.99\times 10^{-5}$) | **不良生存 / 免疫排斥**；驱动致密结缔组织增生，限制效应免疫实质浸润 |
| **#5** | **CAF_5** | CAF_5 (Contractile) | 成纤维 | `TAGLN, ACTA2, MYH11, CNN1, BTNL8` | 血管外周与基质束 | **EMT** (436.83, $p_{\text{adj}}=1.33\times 10^{-10}$)<br>**MYOGENESIS** (59.61, $p_{\text{adj}}=0.00132$) | **谱系追踪支持平滑肌起源**；提供基质机械物理张力，促进肿瘤转移 |
| **#6** | **T_1** | T_0 (Tprolif) | T 淋巴 | `MKI67, TYMS, RRM2, TOP2A, NUSAP1` | 肿瘤核心低纤维化实质 (Tumor core) | **E2F_TARGETS** (1258.46, $p_{\text{adj}}=1.09\times 10^{-31}$)<br>**G2M_CHECKPOINT** (722.92, $p_{\text{adj}}=1.49\times 10^{-21}$) | **克隆扩增**；在低纤维化肿瘤实质内部活跃扩增，与 TAM_1 互为趋化招募 |
| **#7** | **CAF_4** | CAF_4 (Proteolytic) | 成纤维 | `SCIN, MMP11, MMP14, HTRA3, LPL` | 缺氧浸润前沿 | **COAGULATION** (94.34, $p_{\text{adj}}=0.00133$)<br>**KRAS_UP** (30.82, $p_{\text{adj}}=0.0406$) | **IFA 证实 LPL 原位高表达**；高表达 LPL 与肿瘤细胞代谢耦合 (逆 Warburg 假说) |
| **#8** | **CAF_8** | CAF_8 (ECM-breakdown)| 成纤维 | `MMP1, MMP3, PLAU, ZNF280B` | 肿瘤细胞团簇交界面 | **EMT** (560.47, $p_{\text{adj}}=9.10\times 10^{-13}$)<br>**TNFA_NFKB** (83.58, $p_{\text{adj}}=0.000167$) | **局部播散**；在肿瘤细胞团交界处降解 ECM 物理屏障，利于集体侵袭 |

> 对应真实数据底表：[`tables/notable_subtypes_summary.csv`](../../output/CRC/P2_CRC_Xenium/analysis/03_biological_discoveries/tables/notable_subtypes_summary.csv)

---

## 3. 统计提取的有界覆盖范围 (Bounded Extraction Scope)

为了防止分析范围无限膨胀，本交付包的统计提取设定了明确的有界性声明：

1. **焦点差异基因 (Focal DEGs)**：
   * 包含 9 个核心亚群（`T_1, T_3, T_5`, `TAM_0, TAM_1, TAM_4`, `CAF_1, CAF_4, CAF_5, CAF_8`）；
   * 每个亚群提取 Top 15 显著标记基因，共 **135 行**（包含真实 logFC 与校正后 padj），保存在 [`focal_subtypes_top_degs.csv`](../../output/CRC/P2_CRC_Xenium/analysis/01_bioinformatic_analysis/detailed_data/focal_subtypes_top_degs.csv)。
2. **Hallmark 通路过表达分析 (ORA)**：
   * 基于 GSEApy 对 8 个焦点亚群的前 5 条关键通路执行过表达分析，共 **40 行**，包含真实 Combined Score 与 padj（其中 27 条 $p_{\text{adj}} < 0.05$），保存在 [`pathway_scores_by_subtype.csv`](../../output/CRC/P2_CRC_Xenium/analysis/01_bioinformatic_analysis/tables/pathway_scores_by_subtype.csv)。
3. **细胞通讯互作对 (CellPhoneDB)**：
   * 基于 1000 次置换检验提取的 **14 条真实显著配受体对**（$p < 0.05$），覆盖 T 淋巴与单核巨噬谱系，保存在 [`significant_lr_pairs.csv`](../../output/CRC/P2_CRC_Xenium/analysis/01_bioinformatic_analysis/tables/significant_lr_pairs.csv)。
