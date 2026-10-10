# 05. 面向 Agent 消费的数据契约与载荷规范 (Agent Data Contract & Payload Specification)

> **上级入口**：[返回交付框架总览 (README.md)](README.md)  
> **关键关联**：[01. 交付架构与目录形态规范](01_architecture_and_directory_layout.md) ｜ [04. 生物学故事与文献对齐](04_biological_validation_and_stories.md)

---

## 1. 机器消费载荷的设计哲学 (Design Philosophy)

为赋能下游 AI Agent 进行自主科学发现，我们必须为大模型提供一个**既高度结构化、又具备严格事实边界（Factually Bounded）**的数据契约：

1. **紧凑尺寸与 Token 预算控制**：
   全量提取产物 [`agent_input_payload.json`](../../output/CRC/P2_CRC_Xenium/agent_input_payload.json) 大小约为 **43 KB**（约 9,000 Token）。这完美适配主流大模型的输入窗口，避免了由于输入几十万行原始矩阵导致的性能衰减与上下文截断。
2. **纯粹确定性的事实边界 (Zero-hallucination Contract)**：
   载荷中包含的所有基因名、差异倍数、富集得分、受配体对与空间区域，均由提取脚本从底层经过统计检验的 CSV 中无损抽取，杜绝大模型在下游推理时虚构数值。
3. **元数据与溯源登记 (Provenance Registry)**：
   载荷头部显式登记了数据来源路径、亚群宇宙总数、评估范围与管线状态，确保每一项结论均具备完整的计算可追溯性。

---

## 2. Payload 核心 Schema 字典与数据模型

标准化载荷 `agent_input_payload.json` 的顶层结构定义如下：

```json
{
  "sample_metadata": { ... },          // 样本身份、亚群宇宙与溯源登记
  "quantitative_evidence": {           // 定量事实证据库
    "unmasking_impact": { ... },        // 重构去掩盖对比与五层覆盖矩阵
    "top_subtypes_multimetric": [ ... ],// 8 大焦点核心亚群多维特征
    "pathways_ora_hallmark": [ ... ],   // 40 条 Hallmark ORA 真实富集结果
    "cell_cell_interactions": [ ... ],  // 14 条 CellPhoneDB 显著受配体对
    "spatial_colocalizations": [ ... ], // 空间共定位描述性记录
    "hypotheses_evidence": [ ... ]      // 4 大科学假说的多维支撑矩阵
  },
  "generation_targets": [ ... ]        // 规范约定的下游产物目标清单
}
```

### 字段级数据契约对照表

| 顶级键 | 二级键 / 字段路径 | 数据类型 | 物理与统计含义 | 下游 Agent 使用范式 (Prompting Directive) |
|---|---|:---:|---|---|
| `sample_metadata` | `sample_id`, `cancer_type`, `modality` | string | 样本基础生物学与技术身份 | 设定分析报告的全局上下文与标题 |
| `sample_metadata` | `subclusters_universe` | int (25) | 重构解耦出的全谱系亚群总规模 | 约束全量宇宙，防止声称仅有局部亚群存在 |
| `sample_metadata` | `evaluation_scope` | dict | 包含的焦点亚群数 (8)、DEG行数 (135) 等 | 明确当前分析的有界范围，说明未覆盖全量的原因 |
| `sample_metadata` | `pipeline_status` | string | 生成状态 (`complete` / `partial`) | 门禁检验：若为 `partial` 则提示缺口，不强行归纳 |
| `quantitative_evidence` | `unmasking_impact.key_comparisons` | list[dict] | Tprolif (18.6x) 等关键表型的去掩盖证据 | 撰写“重构揭示了什么”章节，论证算法核心价值 |
| `quantitative_evidence` | `unmasking_impact.systemic_impact_5layers` | list[dict] | Unit/Gene/Path/CCI/Niche 五层系统覆盖大表 | 渲染系统化全景矩阵，展示微环境信息解脱 |
| `quantitative_evidence` | `top_subtypes_multimetric` | list[dict] | 8 大焦点亚群标记物、生态位与预后角色 | 生成亚群核心特征卡片，驱动靶标分子识别 |
| `quantitative_evidence` | `pathways_ora_hallmark` | list[dict] | ORA 富集通路名、Combined Score 与 padj | 绘制功能通路气泡图，推演亚群生物学状态 |
| `quantitative_evidence` | `cell_cell_interactions` | list[dict] | CellPhoneDB 置换检验显著受配体对 | 构建微环境分子通讯流与受配体互作网络 |
| `quantitative_evidence` | `hypotheses_evidence` | list[dict] | H1–H4 假说结构化对应关系与文献锚点 | 核心科学假说推演、撰写假说报告与设计实验 |

---

## 3. 标准载荷生成器与健壮性保障 (Generator Script)

载荷由专用无偏提取脚本 [`scripts/prepare_agent_payload.py`](../../scripts/prepare_agent_payload.py) 生成，内置以下安全机制：

1. **依赖底表完整性校验**：
   自动检测必需的 7 张核心底表；若任意底表缺失，自动标记 `pipeline_status: partial` 并列出 `missing_inputs`，拒绝生成伪装成功的载荷。
2. **跨平台确定性输出**：
   统一浮点数解析与 JSON 序列化编码（`ensure_ascii=False`, `indent=2`），保证在 macOS（arm64）与 Linux（x86_64）环境下生成的 SHA-256 绝对一致。

### 命令行调用方式

```bash
# 从指定样本目录生成标准化 Agent 载荷
python3 REVISE_Analysis_Agent/scripts/prepare_agent_payload.py --sample-dir output/CRC/P2_CRC_Xenium
```

---

## 4. 下游 Agent 提示词调用范式 (Downstream Prompting Guide)

当下游 Agent 消费本载荷以生成科学假说或实验设计备忘录时，建议使用以下系统约束提示词：

> **下游 Agent 系统提示词模板 (System Prompt Snippet)**：  
> “你是一个肿瘤空间转录组与微环境科学发现专家 Agent。你将接收一个经过严格统计检验的结构化载荷 `agent_input_payload.json`。  
> 你的推理必须严格遵守以下事实纪律：  
> 1. **严禁捏造未经载荷声明的数值**：报告中引用的每一个差异倍数、通路富集分、受配体互作分必须直接取自载荷中的字段；  
> 2. **尊重证据层级**：明确区分体系内生信计算证据（如 LTB-LTBR 互作分 100.0）与外部文献佐证（如 Cell 论文 16 队列 DFS 预后）；  
> 3. **空间情境锚定**：所有提出的生物学机制必须结合载荷中的 `spatial_niche` 进行解释，说明其发生在哪一个微环境解剖区域。”
