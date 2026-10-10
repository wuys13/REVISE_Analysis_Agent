# 03. 重建影响的五层阅读视图

> [返回交付框架总览](README.md) ｜ [Impact 分析](../analyses/reconstruction-impact.md) ｜ [输入输出协议](../input-output.md)

五层视图用于系统回答“重建改变了什么”，但不是五个自动成功的分析模块。每一层都必须从本次结果读取实际事实，或明确记录不可评估原因。

## 1. 五层问题与当前能力

| 层面 | 要回答的问题 | 当前可复用入口 | 当前边界 |
| --- | --- | --- | --- |
| **ST unit / 细胞亚型** | 单位语义、标签覆盖、组成和局部状态是否改变？ | `reconstruction_impact`、`methods.partition`、`methods.regions` | Raw/SVC 默认独立；成员变化只在明确共同 ID 时比较。unit uncertainty 仍需单独定义和校准 |
| **Gene / 分子** | 哪些基因被观测、恢复、改变或具有特异性？ | `impact_molecular`、`spatial_autocorrelation`、`methods.differential_expression` | 表达 identity、共同基因轴、背景和检验条件必须明确；更多检出量不自动等于灵敏度提高 |
| **Pathway / program** | 基因集覆盖、程序评分和空间分布怎样变化？ | `pathway_activity`、`impact_molecular`、`methods.programs` | provider、基因集和有效基因覆盖是前提；单侧分数与跨侧差值分开解释 |
| **CCI / 通讯** | sender/receiver、ligand/receptor 及其空间支持到什么程度？ | `methods.cell_communication` 与已登记的 producer/case 结果 | wrapper 不是完整 spatial CCI workflow；表达候选、邻近支持和功能因果不可混写 |
| **区域 / niche** | 变化位于哪些坐标、区域、边界或邻域？ | `methods.regions`、Anatomy/State/Gain 结果、空间图 | 区域来源和尺度需登记；坐标和共定位不单独证明 TLS 成熟度、力学、血流或因果机制 |

## 2. 证据链

建议按以下顺序阅读，而不是用一组正向 KPI 代替全层覆盖：

```text
输入身份与单位语义
  -> 标签/表达/坐标的可评估范围
  -> 单层结果与分母
  -> 跨层连接（基因 -> 程序 -> 通讯/邻域 -> 区域）
  -> 支持、反证、替代解释与不可评估项
```

`State` 描述主 SVC 标签的局部多样性；`Gain` 是在共同有效窗口上定义的描述性差值；Moran 在两侧原生空间上独立计算；AUCell/program 先固定单位级评分再做空间聚合；membership 只比较明确的共同单位。它们是平行证据视图，不能单独证明整体生物学改善。

## 3. ST unit 与 gene uncertainty

用户已确认 unit 和 gene 两维 uncertainty 是后续目标。当前文档只规定消费边界，不把实验性稳定性分数改名为正式后验：

- 先声明量针对测量、表达预测、类型映射还是位置；不同对象分别解释。
- 全局 gene 或 unit 分数不足以说明某个局部 gene×unit 组合可靠；需要时保留条件查询能力。
- 生物学异质性、低检测和技术/推断不确定性分开；不把真实罕见状态自动压成噪声。
- Gene、Pathway、CCI 和区域解释都应传播适用的覆盖/可靠性信息，并保留校准和缺失原因。

具体算法、单位、校准和保存形式尚未在本仓库形成正式 contract，不能由当前结果声称已完成。

## 4. 故事与比较边界

重建影响可以显示可观测特征、差异、空间位置和候选机制，但不能把“恢复的基因”“显著的程序”“共定位”直接改写成因果、生物学改善或外部实验已证实。论文 case、临床队列、IFA、谱系追踪和湿实验是外部或历史证据，需保持来源和支持范围，不能作为本次分析的隐式运行前提。

任何层面缺少数据或方法时，结果应保留 `partial`、`unavailable`、`skipped` 或 `failed` 及原因。五层矩阵是导航和审阅工具，正式事实仍来自 `result.json` 及其 outputs。
