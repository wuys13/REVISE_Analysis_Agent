# 按科学问题找到能力

## 正式分析入口

| 名称 | 回答什么 | 需要什么 | 结果与说明 |
| --- | --- | --- | --- |
| `reconstruction_impact` | 表征、成员、局部状态与空间特征发生什么变化 | Raw/SVC 表达；分群用 broad 标签；空间部分用坐标与尺度 | [Impact 主线](reconstruction-impact.md)，分群、区域、Moran、通路及诊断图 |
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

新增分析：先用函数与主 Notebook 探索；需要跨样本复用时下沉稳定计算，再注册一个明确分析入口。说明它的问题、输入、参数、结果与解释限制即可，无需 Agent registry 或通用证据框架。
