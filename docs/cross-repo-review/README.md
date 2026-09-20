# REVISE × Analysis Agent：协同审阅入口

本仓库的初衷是让 Agent **直接完成下游分析，并能围绕新科学问题扩展分析**。从 REVISE 分支选择性迁入方法，是为了保留科学计算能力和经验，同时解除对重建 backend、旧载体和旧批处理的依赖。Impact 是目前最完整的阅读工作台，不代表全部能力。

2026-09-20 双仓 mini 真实联测已完成，当前结果和限制见[现状与缺口](evidence-and-gaps.md)。正式项目配置为 `configs/p2_project.yaml`、`p1_hd_project.yaml`、`p2_visium_project.yaml`，直接引用 REVISE 发布的 sample.yaml。

本说明从消费端组织审阅；[上游审阅入口](../../../REVISE/docs/development/reconstruction-analysis/cross-repo-review.md)从生产端组织。两者引用各自实现，不定义第二套协议。适用版本、证据和未完成工作见[现状与缺口](evidence-and-gaps.md)。

## 一条端到端主线

**REVISE 准备与重建 → 发布 `sample.yaml + raw.h5ad + SVC.h5ad` → 分析库加载并判断能力 → 正式分析或独立探索 → 保存结果与科学审阅。**

| 责任方 | 负责什么 | 边界 |
| --- | --- | --- |
| REVISE | reference 准备、重建、原 Raw 保护、SVC 表达与标签、来源声明、同次交付及失败回滚 | 生产代码能说明当前生成路径，不能反推历史 H5AD 的处理历史 |
| Analysis Agent | 原生对象加载、方法与流程、参数解析、Notebook 探索、批处理、结果与静态报告 | 不导入 `revise`，不改写输入、不补 Raw Level2、不引入旧 spatial/expr 兼容层 |
| Agent／研究者 | 选择问题与范围、检查前提、调用或扩展分析、阅读诊断并提出下一问题 | 不以执行成功代替科学解释，不把 unavailable 当无差异 |

Raw/SVC 默认保留独立单位与基因轴；仅特定问题采用局部明确对应。State 描述 SVC 主标签的局部多样性，Gain 是相对 Raw baseline 的描述性差值；二者都不自动证明重建带来生物学改善。

## Agent 的三种使用方式

1. **直接执行已有分析。** 从[能力索引](../analyses/README.md)选择 Impact、Moran 或 pathway，通过 API、CLI、batch 或连续 Impact Notebook 调用。先检查结果状态与分母，再读表图。
2. **独立探索新问题。** 使用独立 Notebook，优先复用 `methods/`，就近展示参数、中间对象和图。独立科学问题不强塞进 Impact。通用函数不自动提供正式流程的全部前提检查、状态记录与产物登记，调用方需显式负责这些事项。
3. **将成熟分析稳定化。** 可复用计算进入 `methods/`；稳定的问题流程进入 `analyses/` 并接入现有 runner、参数与结果协议。补充行为测试、能力说明及适合该问题的结果阅读，不要求每个临时实验立即成为正式 workflow。

独立 Notebook 的最小约定：说明问题、输入来源及表达前提、预处理、单位/基因/比较范围、随机参数、资源依赖、输出位置与解释限制；保护输入，保留缺失与失败。具体入口见[架构](../architecture.md)、[调用协议](../input-output.md)和[Notebook 说明](../../notebooks/README.md)。不新增 Agent planner、通用缓存或插件注册框架。

## 上游额外比较怎样协同

上游已有[四方法比较 Notebook](../../../REVISE/reproduce/case/assembly_comparison.ipynb)和[比较协议](../../../REVISE/docs/development/reconstruction-analysis/assembly-comparison-contract.md)。它读取 mean、random、within_cluster、outside_cluster 四种已生成表达，各自重新分群，与指定的历史 `SVC_cluster` 标签比较。该标签是比较基准，不是独立生物学真值；ARI/NMI 不自动决定赢家或默认方法。

该实验采用共同 ID/基因范围，只适用于该比较问题，不改变本仓库默认独立分析规则。本轮保持其现有归属并链接复用；将来若迁入分析库，应单独处理依赖拆分、入口和输出协议，避免双重维护。2026-09-20 已完成上游四方法 mini 真实运行，T/Mono_Macro/Fibroblast 共同 ID 为160/157/68、共同基因13088。见[执行版比较 Notebook](../../../REVISE/output/mini-acceptance/20260920/assembly/comparison/assembly_comparison.executed.ipynb)；全量与科学解释尚未验收。

历史讨论中的 Recoverability、TLS 等属于未来方向；CCI、轨迹、富集已有部分可调用工具，但不是全部建立了正式流程或真实验收。旧审阅意见是设计背景，当前能力以源码和[证据边界](evidence-and-gaps.md)为准。

## 阅读顺序与给 GPT 的任务

| 问题 | 下钻位置 |
| --- | --- |
| 双方字段和行为是否真正对应？ | [协议与调用链索引](contract-and-code.md) |
| 哪些完成、哪些未验证、谁负责？ | [证据与缺口](evidence-and-gaps.md) |
| 某项方法的科学定义是什么？ | [能力索引](../analyses/README.md)及专题 |
| 最初迁入了什么、后来为何改变？ | [迁移记录](../migration.md)，历史段落不作为当前状态表 |

可将以下任务直接交给能访问两个仓库的 GPT：

> 同时审阅 REVISE 和 REVISE_Analysis_Agent 当前工作树。先读本入口、协议索引和证据表，再核对双方分支、HEAD、未提交修改与证据版本。沿“生产配置与来源 → 实际交付 H5AD/YAML → loader 与能力判断 → 方法/Notebook/batch → 当前结果与报告”检查，不只比较文档字段名。重点关注 Raw 保护、SVC 表达、主标签与注释的区别、未知来源、独立轴与局部对应、物理坐标、失败与部分可用状态，以及 Agent 扩展是否绕过必要前提。也核对上游独立 assembly 比较的输入、范围及解释边界。
>
> 每个发现给出严重度、实际触发条件、双方文件和符号、影响、可复核证据与建议归属（生产端/消费端/共同约定）。分别列出已确认缺陷、未验证风险和后置科学问题。明确哪些检查是本次执行、哪些引用历史记录。不要把暂缓研究扩展变成当前阻塞，不以测试通过宣称生物学改善。只审阅时不改代码或输入；需要真实运行时先说明具体输入与范围。

包内文档使用相对链接。跨库链接假设两仓库并列；若只提供一个仓库，须同时提供被引用的上游文件或另一仓库访问权限。`output/` 中证据可能被 Git 忽略，缺文件时应重取或明确未核验，不能当作已通过。
