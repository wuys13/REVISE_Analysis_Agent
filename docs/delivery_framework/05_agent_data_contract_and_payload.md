# 05. Agent 载荷与事实契约

> [返回交付框架总览](README.md) ｜ [目录架构](01_architecture_and_directory_layout.md) ｜ [输入输出协议](../input-output.md)

本页把确定性计算与 Agent 文字综合分开。代码读取实际结果并登记来源、范围和状态；Agent 可以据此组织 `FINDINGS.md`、故事或假说草稿，但不能补写缺失数字、改变分母或把候选机制改成已证实结论。

## 1. 载荷的最小内容

载荷不是新的结果数据库，也没有固定大小或 token 目标。适合机器消费的载荷至少应能回答：样本和输入是谁、本次覆盖了什么、哪些证据可读、哪些缺失或不可用、每条事实回到哪里。

```json
{
  "sample_metadata": {
    "sample_id": "<from registered sample>",
    "source_status": "<observed status>",
    "acceptance_level": "<technical or scientific level>",
    "scientific_acceptance": false,
    "evaluation_scope": "<counts and ranges derived from current tables>",
    "provenance_registry": "<optional source registry>",
    "evidence_files": {
      "<table role>": {"path": "<relative path>", "sha256": "<current hash>", "rows": "<current rows>"}
    },
    "pipeline_status": "complete | partial",
    "missing_inputs": []
  },
  "quantitative_evidence": {
    "historical_count_comparisons": "<records copied from the saved comparison table>",
    "pathways_ora_hallmark": "<records copied from the saved ORA table>",
    "preselected_lr_rows": "<records copied from selected LR candidates>"
  },
  "contextual_evidence": {
    "focal_subtype_descriptions": "<saved focal subtype descriptions>",
    "systemic_impact_5layers": "<saved five-layer context>",
    "literature_spatial_colocalizations": "<literature-derived spatial context>",
    "hypothesis_candidates": "<saved candidate hypotheses>"
  },
  "generation_targets": ["<relative paths that a writer may update>"]
}
```

字段值必须来自当前样本和结果；示例中的占位符不构成任何样本事实。载荷可以引用五层阅读视图，但正式状态和输出身份仍以 `index.json`、`result.json` 和 manifest 为准。

## 2. 生成与消费边界

[`scripts/prepare_agent_payload.py`](../../scripts/prepare_agent_payload.py) 负责从已经准备的分析材料读取 metadata、表格和可选 source registry，并登记当前行数、路径和 hash。它不是科学分析 runner，也不触发 Agent 写作；脚本退出成功只表示载荷生成过程完成，不能证明表格内容正确或科学结论已接受。

载荷把保存的数值比较与上下文材料分开：`quantitative_evidence` 只包含历史计数比较、Hallmark ORA 行和预选 LR 行；焦点亚型描述、五层系统影响、论文描述性空间共定位和假说候选放在 `contextual_evidence`。当前 LR 与空间表的文件名分别是 `selected_lr_candidates.csv` 和 `literature_spatial_context.csv`。`evaluation_scope.selected_cci_rows` 是预选记录行数，不能解释为显著配对数量。`scientific_acceptance` 直接复制登记的 `summary.json`；登记缺失时保留 JSON `null`，不能按状态或文件存在推断为 `true` 或 `false`。

生成前应检查：

1. 样本身份来自登记的 `summary.json` 或等价事实来源，而不是脚本常数。
2. 缺少必需材料时状态为 `partial` 并列出 `missing_inputs`；缺失项不能静默变成空数组或零。
3. 表格中的缺失值保持 JSON `null`，不要改变缺失和零的含义。
4. evidence file 的相对路径、当前 hash 和行数与本次读取对象一致；历史路径只作为历史来源标识。
5. `generation_targets` 只列实际存在且允许被写入的阅读材料，不能暗示已经自动生成。

## 3. Agent 写作纪律

- 每条数字、排名、分母、程序或通讯候选都能回到载荷记录和底表；需要更多细节时下钻到 `result.json` 和原表。
- 明确区分观察、派生量、候选解释、历史 case 和外部证据。
- 报告 `partial`、`unavailable`、`failed`、`pending` 与 `scientific_review_pending`，不把它们压成“无差异”或“已完成”。
- 空间邻近、群级表达、类别 posterior 或未校准分数不能直接称为物理接触、功能因果或正确概率。
- 目标文档应保留证据入口、支持范围、替代解释、限制和下一步检查；文字总结不能成为第二份事实数据库。

## 4. 推荐工作顺序

```text
确定性分析/结果发布
  -> index.json / result.json / output manifest
  -> prepare_agent_payload.py（可选摘要）
  -> Agent 读取载荷并按需下钻底表
  -> 生成阅读视图或假说草稿
  -> 人类科学审阅与单独验收
```

三槽位是阅读组织，不是 Agent 自动编排协议。外部临床、IFA、谱系追踪和湿实验材料可以作为关联证据登记，但不自动成为当前载荷或生信运行的隐式前提。
