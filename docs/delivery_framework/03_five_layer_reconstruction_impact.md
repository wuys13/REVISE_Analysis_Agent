# 03. 重建五层系统解构矩阵规范 (Systemic 5-Layer Reconstruction Impact)

> **上级入口**：[返回交付框架总览 (README.md)](README.md)  
> **关键关联**：[02. 客观统计基线与亚群排序](02_unbiased_statistical_ranking.md) ｜ [04. 生物学故事与文献对齐](04_biological_validation_and_stories.md)

---

## 1. 重构算法的核心使命：破除原位物理局限 (The Unmasking Mission)

原始原位空间转录组（以 10x Xenium 为代表）具有极高的亚微米级光学分辨率，但受制于物理探针化学合成上限与光学衍射限制，存在两大本质瓶颈：
1. **靶标基因空间受限与严重 Poisson Dropout**：通常仅有 300~500 个预设探针（本样本为 422 个），导致绝大多数关键转录因子（如 `FOXP3`、`CTLA4`、`MAF`）和细胞周期下游调控分子近乎零检出；
2. **基质与上皮信号光学弥散**：组织致密区高丰度结构分子横向扩散，导致微弱免疫亚群被泛胶原或上皮角蛋白假性遮盖。

**REVISE 超分辨率重建**的本质，是通过引入单细胞参考先验分布，将受限的原位单细胞测序映射恢复至全转录组单细胞真实状态。为了系统化评估这一跃迁，我们建立了**“五层系统化解构覆盖矩阵”**。

---

## 2. 五层系统化覆盖矩阵 (The 5-Layer Coverage Matrix)

矩阵将空间区域（Niche）作为贯穿始终的微环境物理坐标，自底向上全面刻画重构前后的状态跃迁：

| 层面 (Layer) | 原始原位测量状态 (Raw Measurement) | REVISE 重构状态 (SVC Reconstruction) | 关键去掩盖机制与收益 (Key Biological Impact) | 对应空间微环境生态位 (Spatial Niche) | 事实底表支持 |
|---|---|---|---|---|---|
| **1. ST unit / 细胞亚型** | 粗粒度大类混合；弱表型亚群被高丰度细胞淹没 | 解耦出 T (7群)、TAM (8群)、CAF (10群) 共 25 个精细亚群 | 均质间质解耦出收缩型、代谢型、降解型 4 类成纤维细胞 | 全组织微环境覆盖 | [`subtype_degs_full.csv`](../../output/CRC/P2_CRC_Xenium/analysis/01_bioinformatic_analysis/detailed_data/subtype_degs_full.csv) |
| **2. Gene / 分子特征** | 422 目标探针；弱表达转录因子受 Poisson Dropout 零值截断 | **扩增至 13,088 基因**；`FOXP3`/`CTLA4`/`MAF` 恢复特异表达 | 高变基因空间表达恢复；Tprolif 检出 DEG 实现 **18.6 倍增幅** | 交界区 (Treg) / 核心 (Tprolif) | [`raw_vs_svc_unmasking_evidence.csv`](../../output/CRC/P2_CRC_Xenium/analysis/02_reconstruction_impact/tables/raw_vs_svc_unmasking_evidence.csv) |
| **3. Pathway / 功能程序** | 靶标探针过少，无法富集完整的 Hallmark 生物学程序 | 完整映射 Hallmark 癌症特征通路网络与功能活性 | E2F/G2M 增殖、IL2-STAT5 稳态调节及 EMT 基质重塑网络清晰涌现 | 核心 (增殖) vs 边界 (调节) vs 前沿 (基质) | [`pathway_scores_by_subtype.csv`](../../output/CRC/P2_CRC_Xenium/analysis/01_bioinformatic_analysis/tables/pathway_scores_by_subtype.csv) |
| **4. CCI / 细胞通讯网络** | 缺失绝大部分受体与配体，无法推断分子接触通讯 | 基于 CellPhoneDB 置换检验捕获候选配受体通讯流 | 识别出 LTB-LTBR (100.0)、CXCL9/10-CXCR3 (100.0) 及 ANXA1 自分泌 | 交界区/TLS 维持轴 与 肿瘤实质招募轴 | [`significant_lr_pairs.csv`](../../output/CRC/P2_CRC_Xenium/analysis/01_bioinformatic_analysis/tables/significant_lr_pairs.csv) |
| **5. 空间区域 / Niche** | 光学横向弥散模糊了亚群的物理空间分界 | 确立微米级精度的微环境空间共存与组织协同网络 | 确立交界区 TLS、外周环状纤维壳、实质浸润多尺度空间结构 | Interface / Desmoplastic rim / Tumor core | [`colocalization_zscores.csv`](../../output/CRC/P2_CRC_Xenium/analysis/01_bioinformatic_analysis/tables/colocalization_zscores.csv) |

> 对应真实数据底表：[`tables/systemic_impact_5layers.csv`](../../output/CRC/P2_CRC_Xenium/analysis/02_reconstruction_impact/tables/systemic_impact_5layers.csv)

---

## 3. 核心表型层面的去掩盖定量对比 (Phenotype-level Unmasking Facts)

根据真实接受性笔记本（`Xenium_sc_SVC_T.ipynb`、`Xenium_sc_SVC_Monocyte.ipynb`、`Xenium_sc_SVC_Fibroblast.ipynb`）中的配对与表型计算，重构带来的关键生物学去掩盖事实如下：

```
[原始原位 Xenium (422 探针)]                [REVISE 超分辨率重建 (13,088 基因)]
───────────────────────────────────         ────────────────────────────────────────
• Tprolif vs TRM: 仅检出 123 个 DEG    ───>  • Tprolif vs TRM: 检出 2287 个 DEG (18.6x 提升)
  (MKI67 下游周期分子严重截断)                    (恢复 TYMS, RRM2, TOP2A, NUSAP1 等完整网络)

• Treg vs TRM: 仅检出 64 个 DEG        ───>  • Treg vs TRM: 检出 90 个 DEG (1.41x 提升)
  (弱表达调控因子原位难以分离)                     (恢复 FOXP3, TIGIT, CTLA4, LTB 特异性)

• 均质促癌间质 (单一纤维化大类)          ───>  • 解耦出 4 类功能正交成纤维细胞亚群:
  (泛胶原探针无法区分精细亚态)                    - CAF_5: ACTA2/MYH11+ (平滑肌收缩)
                                                  - CAF_4: LPL/MMP11+ (代谢耦合)
                                                  - CAF_8: MMP1/MMP3+ (基质降解)
                                                  - CAF_1: COMP/COL11A1+ (ECM 硬化)
```

### 证据溯源索引

* **T 淋巴细胞表型对比**：`REVISE/reproduce/case/Xenium_sc_SVC_T.ipynb` (Cell 28) 对应 Cell Fig 6D–G；
* **单核巨噬表型解离**：`REVISE/reproduce/case/Xenium_sc_SVC_Monocyte.ipynb` (Cell 30) 对应 Cell Fig 6H–I；
* **成纤维细胞异质性解构**：`REVISE/reproduce/case/Xenium_sc_SVC_Fibroblast.ipynb` 对应 Cell Fig 7A–G。
