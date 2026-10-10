# 01. 交付架构与目录形态

> [返回交付框架总览](README.md) ｜ [输入输出协议](../input-output.md) ｜ [分析能力索引](../analyses/README.md)

## 1. 阅读视图与正式结果

交付框架采用由外向内的渐进式披露，但只定义阅读顺序：

```text
L0  工作区/批次导航       项目、样本和当前结果入口
L1  样本导航               output/<sample_id>/index.json
L2  分析状态               <analysis>/result.json、参数、阶段状态和限制
L3  统计事实               result.outputs 登记的 tables/*.csv
L4  细节与展示             figures/、report.html、按需明细
```

`01_bioinformatic_analysis`、`02_reconstruction_impact` 和
`03_biological_discoveries` 是把同一组事实组织给读者的三个槽位。它们不要求写成正式输出目录，也不允许绕过 `index.json`、`result.json` 和输出 manifest 另建事实登记。

## 2. 正式输出树

```text
output/
└── <sample_id>/
    ├── index.json                         # 样本级结果导航
    └── <analysis>/
        ├── result.json                    # 参数、状态、限制、输出 manifest
        ├── tables/                         # 本次结果登记的统计表
        ├── figures/                        # 本次结果登记的图
        └── report.html                    # 只读已保存结果的报告（如有）
```

实际发布的分析名称由 `runner.ANALYSES`、项目配置和样本执行决定。目前正式入口见[分析能力索引](../analyses/README.md)；缺少的 CCI、TLS 或其他后续能力不能因为阅读槽位存在而被假定为已经运行。

结果目录可能包含受保护的历史 `.previous/` 或独立 review 证据。它们保留原始运行身份，不自动进入当前 `index.json`，也不能被当前报告当作本次产物。

## 3. 内容职责

| 层级 | 内容 | 使用方式 |
| --- | --- | --- |
| L0/L1 | 项目、样本、`index.json`、报告入口 | 先确认样本、提交、来源和当前状态 |
| L2 | `result.json`、参数、阶段状态、错误或不可用原因 | 判断本次到底运行了什么 |
| L3 | `tables/*.csv`、来源和比较范围 | 做定量复核和二次筛选 |
| L4 | `figures/`、报告和按需明细 | 解释空间上下文、局部现象和完整证据 |

报告只能读取本次 `result.outputs` 声明并实际存在的文件。目录中残留的旧表、截图或 Notebook 不因路径相近而成为当前事实。

## 4. 载荷与批量扩展边界

`agent_input_payload.json` 是可选的机器消费摘要，不是独立事实源。它应携带样本/输入身份、来源登记、实际范围、证据文件、状态和缺失原因；下游需要时再回到 `result.json` 与底表核查。载荷大小、token 窗口和字段数量随样本与消费任务变化，不在这里设置固定目标。

批量运行沿用现有 `batch.py` 的样本 × 分析执行和失败隔离。跨样本汇总必须登记独立样本、可比口径、结果版本和缺失/失败状态；不能只把多个 `summary.json` 拼成新的总结果。新的科学问题先在独立 Notebook 或方法函数中探索，成熟后才进入 `analyses/` 和现有发布协议。

## 5. 事实与证据规则

- 代码负责确定性计算和输出登记；阅读文档和 Agent 可以组织解释，但不能构造统计量。
- 输入缺失、身份未知或方法不支持时，保留 `partial`、`unavailable` 或 `failed` 及原因，不补零、不用随机值代替。
- Raw/SVC 默认独立；共同基因、共同单位和成员变化仅在实际比较条件满足时使用。
- 技术结果、科学审阅和生物学接受分开记录。报告存在、测试通过或进程退出 0 都不自动等于科学验收。
