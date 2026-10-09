# 01. 评测对比域规范（Benchmarks Specification）

## 一、 领域职责与科学边界

`00_benchmarks/` 目录专职回答**“重建算法的质量、抗噪性与保真度是否达标”**。

### 科学原则
1. **严格物理/统计量配对**：仅对在 Raw 与 SVC 之间具有对等物理定义或共同度量空间的指标开展定量对比；
2. **多方法横向基准（Multi-method Benchmark）**：必须将 `Raw`、`Random donor`、`Mean prior` 与 `Within_cluster OT` 纳入同一比较框架，用实测数据证明当前优化配置是否最优；
3. **消除白噪声是核心门禁**：通过平滑先验压制离散随机抽样引入的高频白噪声，是判定重建成功的关键技术门禁。

---

## 二、 子槽位规范与计算契约

### 1. `01_spatial_continuity/`（空间连续性与白噪声消除）
* **分析对象**：
  - 全基因组空间自相关统计量（全局 Moran's I 与 Geary's C）；
  - 局部信噪比与空间梯度衰减率。
* **对齐口径**：
  - **共同基因轴（Common Genes）**：锁定 Raw 与各 SVC 方法完全共有的 343 个基因（严格相同大小写、顺序与过滤条件）；
  - **全量基因拓展轴（All Expanded Genes）**：评估 13,000+ 恢复基因的空间连续性是否与生物学先验自洽。
* **产物契约**：
  - `summary.json`：
    ```json
    {
      "scope": "common_343_genes",
      "methods_compared": ["raw", "random", "mean", "within_cluster"],
      "raw_median_moran": 0.312,
      "mean_median_moran": 0.460,
      "mean_moran_delta": 0.148,
      "positive_gain_fraction": 0.872,
      "white_noise_suppression": "validated",
      "verdict": "pass"
    }
    ```
  - `tables/moran_raw_vs_methods.csv`：
    表头包含 `gene_symbol`, `moran_raw`, `moran_random`, `moran_mean`, `moran_within_cluster`, `delta_mean_vs_raw`, `p_val_raw`, `p_val_mean`。

---

### 2. `02_microenvironment_fidelity/`（微环境物理状态与 Gain 检验）
* **分析对象**：
  - 空间物理多尺度窗口（50μm, 100μm, 150μm）下的局部细胞组成多样性（Shannon Entropy、Simpson Index）；
  - 经过严格 K-control 控制网格校正后的微环境状态（State）与增益（Gain）候选。
* **对齐口径**：
  - 使用相同的空间原点（Spatial Origin）和统一网格切片；
  - 严格剔除由聚类群数（Cluster K）人为不平衡引入的假阳性熵增。
* **产物契约**：
  - `summary.json`：
    ```json
    {
      "window_scale_um": 100,
      "raw_mean_entropy": 1.25,
      "svc_mean_entropy": 1.48,
      "entropy_delta": 0.23,
      "k_matched_control_status": "consistent",
      "true_gain_candidate_count": 3
    }
    ```
  - `tables/diversity_entropy_delta.csv`：各物理窗口的熵值配对对比表。
  - `tables/k_control_gain_audit.csv`：微环境候选 Gain 区域的 K 值敏感度与 bootstrap 检验表。

---

## 三、 评测判定（Verdict）机制

在 `00_benchmarks/FINDINGS.md` 中，系统自动依据以下阈值输出结构化裁定：
1. **Pass（通过）**：
   - 共同基因轴中 Moran's I 正向提升比例 $\ge 75\%$；
   - 且不存在显著的空间几何变形或坐标错配。
2. **Conditional / Warning（带限制通过）**：
   - 共同基因轴正向提升在 $50\% \sim 75\%$ 之间；
   - 或局部微环境增益在 K-control 网格边缘出现波动。
3. **Fail（未通过）**：
   - 类似早期 random 模式出现的 Moran's I 全局暴跌（正向率 $< 10\%$）；
   - 提示单细胞离散白噪声严重，阻断下游精细生物学发现模块的信任级别。
