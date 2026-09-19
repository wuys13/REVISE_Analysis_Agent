# 架构

## 边界

REVISE／上游交付 analysis-ready Raw 与统一 SVC；本仓库从这两个对象开始。正式 `.X` 契约固定为 finite、nonnegative、unlogged linear；每侧分别声明 identity，允许 unknown。loader 不从数值猜测；Leiden、Moran、AUCell 在消费表达时检查已确认身份。主标签和空间分析不依赖表达。Raw Level2 由上游提供，缺失只影响依赖它的分析。

## 分层

| 位置 | 负责什么 |
| --- | --- |
| `revise_analysis/io.py` | 读取样本配置与两个独立对象，不配对、不归一化 |
| `methods/` | 可复用计算及必要中间结果 |
| `analyses/` | Impact、Moran、通路活性的稳定组合 |
| `plotting/` | 消费计算结果的诊断与结果图 |
| `runner.py`、`batch.py` | 单任务、sample × analysis、结果发布 |
| `reporting/` | 读取已保存表、图和结果索引，不读 H5AD |
| `notebooks/` | 连续的科学开发工作台，调用相同参数解析、阶段对象和保存协议 |

## 简化原则

默认不配对；只有 ARI、Hungarian 匹配和成员变化定位使用同单位对应。Moran 在两侧原生坐标分别建 6 邻居图；通路分别评分；Anatomy 用完整 Raw 独立网格，parent State/Gain 用各自配置的网格。跨网格将每个观察点落入 Anatomy 再汇总组成，不直接按 window_id 连接。State 是主 SVC 标签的局部多样性；ΔNeff vs Raw Leiden 是不同标签系统的描述性差值。

不建立通用 pairing、preflight、Region class、Agent planner、缓存指纹或旧 batch 兼容框架。Agent 通过能力说明与 `index.json → result.json → tables/figures/report.html` 工作。

## 计算顺序与阅读顺序

科学阶段保留真实依赖；报告采用固定的 Impact 问题树：重建后特征与两侧差异、空间位置、空间上的分子与成员关联。阶段状态用于解释执行情况，不直接决定网页章节。Moran 和固定单位级 AUCell 评分是并行证据，空间关联消费已有评分和区域；integration 表是关系汇总，不是新的独立验证。

关系图只读取本次 outputs 登记的保存表，由 Notebook/batch 共用展示函数生成。每张图登记到负责它的阅读 section；科学阶段产物与图登记分开，因此重绘不会伪装成重新计算。HTML 渲染只读取登记结果、图和底表；不能为了调整阅读结构重新调用科学方法或输入 loader。完整文件索引保留为审计入口。

Impact 的交互状态是显式的：完整 `STAGE_ORDER` 始终出现在结果中，未运行或参数变化后失效的阶段是 `pending`。参数影响图定义每个阶段真正依赖的配置；失效沿该图传播，`run_stage` 同时检查前置状态。这里是该 workflow 的明确依赖，不扩展成通用调度框架。Notebook 运行前重新解析 sample/project/override 并要求先应用不一致的参数，避免用新参数解释旧产物。
