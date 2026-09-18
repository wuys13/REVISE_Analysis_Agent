# 架构

## 边界

REVISE／上游交付 analysis-ready Raw 与统一 SVC；本仓库从这两个对象开始。输入 `.X` 为未归一化、未取对数的非负表达，重建值可为浮点数。Raw Level2 由上游提供，缺失只影响依赖它的分析。

## 分层

| 位置 | 负责什么 |
| --- | --- |
| `revise_analysis/io.py` | 读取样本配置与两个独立对象，不配对、不归一化 |
| `methods/` | 可复用计算及必要中间结果 |
| `analyses/` | Impact、Moran、通路活性的稳定组合 |
| `plotting/` | 消费计算结果的诊断与结果图 |
| `runner.py`、`batch.py` | 单任务、sample × analysis、结果发布 |
| `reporting/` | 读取已保存表、图和结果索引，不读 H5AD |
| `notebooks/` | 连续的科学开发工作台，调用相同计算模块 |

## 简化原则

默认不配对；只有 ARI、Hungarian 匹配和成员变化定位使用同单位对应。Moran 在两侧原生坐标分别建 6 邻居图；通路分别评分；空间窗口可共用而单位与抽样独立。State 与 Gain 不是同一概念。

不建立通用 pairing、preflight、Region class、Agent planner、缓存指纹或旧 batch 兼容框架。Agent 通过能力说明与 `index.json → result.json → tables/figures/report.html` 工作。
