# 按科学问题找到能力

## 正式分析入口

| 名称 | 回答什么 | 需要什么 | 结果与说明 |
| --- | --- | --- | --- |
| `reconstruction_impact` | 表征、成员、局部状态与空间特征发生什么变化 | SVC 主标签、scope 标签及空间坐标/尺度支持标签空间链；表达仅为相应 baseline/分子分支前提 | [Impact 主线](reconstruction-impact.md)，分群、区域、Moran、通路及诊断图 |
| `spatial_autocorrelation` | 基因在各自原生空间单位上的自相关 | 各侧表达与空间坐标 | [Moran](spatial-autocorrelation.md)，每侧基因表和比较图 |
| `pathway_activity` | 基因集覆盖和通路评分分布有什么差异 | 显式基因集、AUCell 可选依赖 | [通路活性](pathway-activity.md)，覆盖、原生单位分数和分布图 |

调用见[输入输出协议](../input-output.md)。默认分析 scope 是逻辑标签，不是重建文件夹。Impact Notebook 支持 sample/project/override 参数层和可失效阶段；参数变化后只重跑实际受影响的阶段。正式入口不会自动运行新的 CCI、TLS 或其他跟进分析。

## 可复用计算工具

这些是 Python 函数，不代表已经建立相应批量工作流或完成真实样本验证。

| 任务 | 模块 | 主要入口／资源 |
| --- | --- | --- |
| 独立分群与局部成员比较 | `methods.partition` | `compute_partition`、`compare_membership` |
| 原生空间 Moran | `methods.spatial` | `compute_moran`、`native_knn_graph` |
| 空间窗口、区域与背景 | `methods.regions` | `compute_window_diversity`、`select_window_scale`、`select_region_threshold` |
| AUCell | `methods.programs` | `compute_pathway_scores`；安装 `[pathway]` |
| 差异表达 | `methods.differential_expression` | `get_degs` |
| 基因集评分与 GMT | `methods.gene_set_scoring` | `score_genes`、`read_gmt` |
| 富集 | `methods.enrichment` | `get_enrichment_local`；显式背景／基因集及 gseapy |
| CellPhoneDB | `methods.cell_communication` | `run_cellphonedb_v5`；安装 `[cci]`，另提供数据库 |
| Palantir | `methods.trajectory` | `infer_palantir`；安装 `[trajectory]` |
| 生物学／聚类指标 | `methods.biological_metrics`、`methods.metrics` | 具体函数文档说明输入标签、图和意义 |
| 保留的旧通用分析 | `methods.unsupervised`、`methods.spatial_autocorrelation` | 原有单对象工具；后者需 `[spatial]`，不等同新正式 Moran 入口 |

## 从独立探索到稳定分析

先选择现有正式流程或通用函数。新科学问题可用独立 Notebook 探索，说明输入来源、表达前提、预处理、单位/基因范围、资源、随机参数、输出与局限；不必追加到 Impact 主 Notebook。通用函数不自动执行正式流程的全部前提检查或结果登记。

成熟计算复用 `methods/`；需要重复运行的问题组合进入 `analyses/`，接入 `runner.ANALYSES` 及现有参数、状态、保存和阅读协议，并补充相应行为测试。临时实验无需提前注册。历史 Recoverability/TLS 等仍是扩展方向，已有 CCI/轨迹工具也不等于真实流程已完成。上游四方法比较和跨库职责见[审阅入口](../cross-repo-review/README.md)。无需 Agent registry 或通用证据框架。
