# REVISE 生信分析槽位体系与科学发现规范索引

本套规范定义了 REVISE 与 REVISE Analysis Agent 在全流程重建后，下游生信分析结果的文件组织、渐进式披露以及面向 AI Agent 自主科学发现的统一标准。

## 规范文档列表

1. [00. 生信分析槽位总览与渐进式披露范式规范](00_overview_and_progressive_disclosure.md)
   - 双域架构（`00_benchmarks` 与 `01_discoveries`）的职责分工
   - 4 级渐进式披露契约（Level 0 `README.md` $\rightarrow$ Level 1 `summary.json` / `FINDINGS.md` $\rightarrow$ Level 2 `tables/` $\rightarrow$ Level 3 `figures/` / `report.html`）
   - 单点入口 `analysis.json` 与顶层机器交互协议

2. [01. 评测对比域规范（Benchmarks Specification）](01_benchmarks_specification.md)
   - 空间连续性（Moran's I / Geary's C）与白噪声消除
   - 微环境局部物理状态保真度与 Gain 检验
   - 严格的 343 共同基因轴与多方法横向基准（Raw vs Random vs Mean vs Within_cluster）
   - 评测门禁判定机制（Pass / Conditional / Fail）

3. [02. 科学发现域规范（Discoveries Specification）](02_discoveries_specification.md)
   - 全转录组无偏分析核心原则（不强行与 Raw 进行无意义差值对比）
   - 锁定先验三大主群（Fibroblast, Mono_Macro, T），内部细分亚型（sc-SVC subtype）无指向性全量计算
   - 统计显著性指标驱动自动筛选（`priority_score`, 空间梯度, 互作置信度）
   - 4 大子槽位（差异表达 DEGs、功能通路、受配体通讯 CCI、空间生态位与共定位）的数据模型与底表表头标准

4. [03. Agent 自主科学发现工作流规范（Agent Discovery Workflow）](03_agent_scientific_discovery_workflow.md)
   - 四阶段自主发现状态机（质控门禁 $\rightarrow$ 假说种子 $\rightarrow$ 靶向底表检索 $\rightarrow$ 证据链综合）
   - 轻量 Python/Pandas 代码调用范式
   - 科学发现备忘录（Discovery Memo）自动产出标准

---

## 快速导航规则

```text
analysis/
├── analysis.json                          # 第一步：Agent 检视全局状态
├── 00_benchmarks/summary.json             # 第二步：Agent 校验重建质量与抗噪门禁
├── 01_discoveries/summary.json            # 第三步：Agent 提取高亮假说种子
└── 01_discoveries/<aspect>/tables/*.csv   # 第四步：Agent 靶向抓取多维证据链
```
