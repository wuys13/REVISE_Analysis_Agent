# REVISE 分析交付框架

本框架说明如何阅读和组织已经发布的分析事实。它是人类与 Agent 的阅读视图，不创建第二套科学计算流程，也不替代正式结果索引。

## 先看哪里

正式结果的唯一导航链是：

```text
output/<sample_id>/index.json
    -> <analysis>/result.json
    -> result.outputs 中登记的 tables、figures、report
```

正式分析能力、输入前提和参数见[分析能力索引](../analyses/README.md)与[输入输出协议](../input-output.md)。连续科学主线见[Impact Notebook](../../notebooks/3_standard_pipeline/01_reconstruction_impact.ipynb)。历史来源与验收限制见[迁移与验收](../migration.md)。

## 三个阅读槽位

下列槽位用于把同一批正式结果组织成可读主线。它们是阅读视图，不要求在 `output/` 中生成同名目录，也不改变现有 `index.json` / `result.json` 责任。

| 阅读槽位 | 回答的问题 | 主要事实来源 |
| --- | --- | --- |
| `01_bioinformatic_analysis` | 当前样本和各侧有哪些标签、基因、程序、通讯或区域事实？ | 对应正式分析的 `result.json`、完整 tables 和 figures；缺少的能力保留不可用原因 |
| `02_reconstruction_impact` | Raw/SVC 在 unit、gene、pathway、CCI 和区域层面发生了什么变化？ | `reconstruction_impact`、`spatial_autocorrelation`、`pathway_activity` 及其明确可比范围 |
| `03_biological_discoveries` | 哪些观察可以形成可检查的生物学解释或假说？ | 已保存事实、空间上下文、外部来源和科学审阅；候选机制不升级为因果结论 |

一个正式分析可以被多个槽位引用；槽位之间不是强制串行依赖。上游重建方法比较和 benchmark 继续归属并列的 `REVISE` 仓库，本仓库只链接其真实结果或 contract，不复制静态模拟 Notebook。

## 渐进式披露

```text
L0  工作区/批次入口       项目、样本和结果索引
L1  样本入口               index.json、结果状态、报告入口
L2  分析入口               result.json、参数、阶段状态、限制和输出 manifest
L3  统计底表               tables/*.csv、可复核的元数据和来源
L4  细节与图               figures/、report.html、按需读取的明细
```

越深入，数据越接近具体输入和计算；越靠外，内容越适合导航。摘要、载荷和报告必须引用同一批已登记结果，不用手写数字或静态示意值补足缺口。

## 六个专题页面

| 页面 | 用途 |
| --- | --- |
| [01. 目录架构与组织形态](01_architecture_and_directory_layout.md) | 说明阅读视图与正式结果树的关系，以及发布边界 |
| [02. 统计基线与范围](02_unbiased_statistical_ranking.md) | 说明排序、全量覆盖和可评估范围如何从实际结果读取 |
| [03. 重建影响五层视图](03_five_layer_reconstruction_impact.md) | 组织 unit、gene、pathway、CCI、区域五个问题，不预设结果 |
| [04. 生物学故事与证据边界](04_biological_validation_and_stories.md) | 区分观察、解释、候选假说和外部验证 |
| [05. Agent 载荷与事实契约](05_agent_data_contract_and_payload.md) | 说明来源、状态、缺失和 Agent 消费方式 |
| [06. 假说生命周期与队列扩展](06_hypothesis_lifecycle_and_cohort_expansion.md) | 记录从计算观察到后续验证的边界和未来扩展 |

## 稳定原则

1. **事实同源**：数字、标签、程序和结论必须能回到本次结果的输入、参数、方法和输出文件；没有来源就写不可评估或未知。
2. **状态可见**：`succeeded`、`partial`、`skipped`、`failed`、`pending` 和 `unavailable` 保留其实际含义，不用成功退出码替代科学验收。
3. **证据分层**：技术计算、历史运行、科学解释和外部文献分别标识；外部临床、IFA、谱系或湿实验材料不成为当前生信管线的隐式运行依赖。
4. **空间有上下文**：unit、gene、pathway 和 CCI 的解释尽量连接实际坐标、区域和单位语义；表达候选与物理接触、相关性与因果性分开。
5. **不默认扩大范围**：没有统一的样本数、亚群数、基因数、载荷大小或 token 窗口。每次交付的范围由样本配置和结果 manifest 声明。
