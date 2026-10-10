# REVISE 架构与科学发现基建：审阅事实与讨论稿

> [!NOTE]
> **已归档历史审阅记录（Archived Review: 2026-10-10）**
> 本目录记录 2026-10-10 架构审阅与脚本整改过程。文中提及的硬编码历史脚本（如 `run_unbiased_all_subcluster_analysis.py`）已在提交 `2759f93` 中清退；当前维护记录见[内容整理跟进](../../../output/maintenance/20261010-content-followup/README.md)，本目录仅供审计溯源与历史依据参考，不代表当前活跃工作流。

该时点的判断见[第三轮只读复核与收口反馈](THIRD_PASS_REVIEW.md)：执行方当时更新至 `6b158ed`。下方及第二轮文件保留原时点事实；旧 payload 截断/常数描述是该轮观察，不作为当前代码、部署或科学验收状态。当前维护记录见[内容整理跟进](../../../output/maintenance/20261010-content-followup/README.md)。

后续入口：[第二轮审阅](SECOND_PASS_REVIEW.md)、[第二轮执行反馈](SECOND_PASS_EXECUTION_FEEDBACK.md)、[合作方开发说明](COLLABORATOR_BRIEF.md)。用户已确认并行开发，并明确当前体系止于生信分析与生物学假设；临床/IFA/谱系只与假设关联，属于体系外后期验证。下文为第一轮审阅记录，其较宽的验证整合建议已由上述新反馈修正。

该时点状态为：`scope_clarified / feedback_draft_ready`。用户当时已明确系统化层级、uncertainty 方向与 case 验收范围；Fig6/7 对应材料已定位。具体 uncertainty 方法和实施细节尚未定稿。本文件保留该时点核查事实，并按当时共识整理建议；可交给执行 Agent 的组织化反馈见 [EXECUTION_FEEDBACK.md](EXECUTION_FEEDBACK.md)。

本轮目标：核查本地与 qz 实际落实，判断交接材料是否服务于重建信息评估与可检验科学发现，形成与用户讨论后的执行反馈。只读核查业务实现；本目录仅保存审阅说明。未授权删除旧内容、同步部署或重跑正式科学计算。

## 1. 目标理解与初步判断

两个目标应相互支持：

1. **重建的信息评估**：系统覆盖 gene、Pathway、CCI、细胞/亚型及具体空间区域，分层组织各侧状态与前后变化。保留、扩展、削弱或改变是每个层面的比较内容，不能替代这些层面的完整覆盖。
2. **自主科学发现**：近期让 Agent 基于 Fig6/7 的既有具体 case，基本完成几类亚型及关键基因的生物学故事；组织真实计算、空间定位与机制解释，保留可检验预测与下一问题。外部队列验证作为额外证明，不作为当前交付的必备门槛。

用户补充的 uncertainty 方向：从 ST unit 和 gene 两个维度刻画不同 ST 技术，并融入重建及后续置信度。方向已明确，具体估计对象、公式、校准与传递方式尚不擅自定为既有实现。

三层阅读视图与渐进式披露值得保留。当前最高优先级是事实可靠性：现有新增脚本将若干直接赋值、人工选择及随机生成的数字写为统计证据。代码确定性、文件齐备、网页可读和远端同步均不能消除这个问题。

## 2. 已核实的核心问题

以下行号对应本次读取的工作树；发现由负责各范围的子代理与主代理整合复核。没有重跑完整生信分析。

| 问题 | 直接依据 | 对最终目标的影响 | 建议 |
| --- | --- | --- | --- |
| 排名并非统一指标自然产生 | `scripts/run_unbiased_all_subcluster_analysis.py:98` 起给 TAM 指定分数，`:150` 起给 CAF 指定分数，`:219` 起给 T 指定分数；声明的 `compute_priority_score` 没有用于这套排序 | 不能以“前八名命中已知亚群”证明无偏或后验验证；排名已被先验选择影响 | 清除亚群专属调分；真实计算所有可评估亚群并保存组件、参数、缺失及排序稳健性 |
| T DEG 数值与 p 值直接指定 | 同脚本 `:195-255`，预定义细胞数、marker、FC 和 p 值；写入 `evidence_source=T/acceptance_results.json` 不等于实际读取该来源 | 不能把这部分表格当本样本检验结果 | 建立输入与实际检验链；演示数据显式隔离标为 synthetic |
| 通路、CCI 与 Z 值没有对应计算 | 同脚本 `:310-355` 直接给通路及 LR 数字；`:364-386` 共定位仅覆盖十个指定亚群，关键 Z 直接赋值，其他值随机生成 | 不能声称全 24 亚群通路、空间限制性 CCI 与共定位已全量算完 | 使用有科学前提的真实方法，保存基因集/LR 库版本、图、尺度、零模型、检验和多重校正 |
| 重建改善数字由整理脚本生成 | `scripts/organize_three_tier_analysis.py:85-124` 直接构造 Moran、降噪、DEG 增益、污染与熵/模块度表 | 输出组织器承担了科学结论的制造，不能用于评价重建 | 整理器只引用已登记结果，不产生科学数值；无结果就明确不可用 |
| 多指标并非独立计算 | 同脚本 `:164-167`，SNR 是 `priority_score * 0.42 + 2.5`，compactness 按亚群直接指定 | 无法把三个指标的一致视为三条独立证据 | 对指标给出科学定义与真实输入；不把同一指标变换重复计为证据 |
| 顶层与底表口径冲突 | 当前底表 CAF_1 排名 1、分数 48.5；样本 `summary.json` 则 T_5 排名 1、分数 14.85 | 同一载荷与网页消费不同事实口径 | 顶层摘要从同次登记表生成；检查数值、单位、排名、分母一致 |
| 载荷包含固定 P2 内容并省略关键限制 | `scripts/prepare_agent_payload.py:25-45` 固定癌种、模态、50,626、24 亚群及 +0.234/20.7x/62.1%；`:63-70` 固定四份假说目标 | 新样本可能继承 P2 数字；缺文件仅跳过，不携带不可评估原因 | 从当前样本身份和正式结果提取；载荷携带执行状态、来源性质、证据定位和限制 |

20.7x（若是多比较均值）与 23.7x（单一比较）可以有不同口径，本身不构成矛盾；本轮问题是这些数字的真实计算链没有由新增实现建立。

科学子代理补核：120 行为各亚群前五行组成，现有 CSV 只有 31 行满足 `p_adj ≤ 0.05`，其中 T 的 30 行 p 值由代码直接指定；不能把全部 120 行称为显著 DEGs，更不能称全量底层差异表达结果。五个预设比较的倍数平均值为 `20.66≈20.7`，不代表 24 亚群总体灵敏度。

上述 31 行的准确分解为：从旧 Fib/Mono CSV 读取的 90 行只有 1 行（CAF_7，p_adj≈0.009198）达到 0.05；T 的 30/30 行是脚本指定的极小 p 值。即 `1/90 来源表行 + 30/30 指定值行`。这里“来源表”不授予当前 run 的统计有效性，仍需核查原检验与来源身份。

输入范围另核：`REVISE_Analysis_Agent/data/P2CRC_Xenium/sample.yaml` 所属旧交付的本地 `SVC.h5ad` 为 `17,455×343`、仅 Fibroblast、10 个 SVC_cluster（SHA256 `d20c8c7bc72c93146e332674f8241c2d550a3b0e2b0c0b817dca3a9f923edc15`）。这是独立存在的另一份输入，并非新 31 文件展示树已绑定的输入。不能将它、2026-09 的 50,626 完整交付、2026-08 的 Fib/Mono/T CSV 与新树的 24 focal 标签未经 manifest 直接拼为同次结果。

### 本地与 qz 实际落实

| 对象 | 本地 | qz 本次观测 | 结论 |
| --- | --- | --- | --- |
| REVISE 代码 | `revise-2.0 @ b4ea718` | `revise-2.0 @ 95c4cf2429bad6374a49ac9b0d8770e5dd0e8bbe`，子仓 clean | 最新配置改动未对齐；不能由本地提交推定远端运行采用它 |
| Analysis Agent 代码 | `sl @ 6d41babb84ffc14d7c4024981d5966bb1db150ce` | `sl @ e33d12217ec1ebed55ee607ed365dae1535e883a`，子仓 clean | 缺 `56df57f/6d41bab` 对应的两份新规范及 prepare/build 脚本；已存在的 unbiased/organize 两脚本 hash 与本地相同 |
| 新展示产物 | 父级 output 与 REVISE/output 两树各 31 文件，全部相同 | REVISE/output 中同名 31 个展示文件逐文件 hash 匹配，missing/extra/diff 均为 0 | 只证明这份展示包已同步 |
| 新样本整体目录 | 两份本地展示树没有三文件输入 | 远端 REVISE/output 实际有 52 文件，约 8.05GB，额外 21 文件为输入、方法和运行元数据 | 两端整体样本目录不完全相同；不能称全目录 100% 对齐 |
| 远端父级 output | 本地已有 31 文件 | qz `REVISE_code/output/CRC/P2_CRC_Xenium` 仅空 analysis 目录 | 同名根路径不是同一交付来源；需选择唯一事实入口与展示索引 |
| 远端运行占用 | 本轮无正式计算 | 当次确认 RUNNING 的 cpu-download；未观察到目标重建或新分析脚本进程 | 这是有限观测时点的进程检查，不代表全部 notebook 永久无任务 |
| 远端重建记录 | 本轮不重跑 | `reconstruction.json` 指向 `code_commit=0abf8515...`、REVISE `.venv/bin/python` 3.10.12，`technical_verified / engineering / scientific_acceptance=false` | 属于历史工程记录；不是最新配置或新科学计算的验证 |
| 远端旧记录证据 | 旧审阅文档可读 | 两个源输入 hash 现场重新匹配，但记录所指 receipt、delivery-audit 和 donor assignments 未找到 | 输入身份可依赖到已核范围，完整回执链尚缺；不得把元数据存在当完整 formal receipt |

远端父工作区 `main @ 6d40a295` 的 `workspace.yaml` 存在未提交修改；本地生产仓另有既存未跟踪文件。后续版本对齐须先比较提交关系、保护这些修改，不直接覆盖或重置。旧目录删除、远端部署操作的历史过程未在本轮复现，不能从当前目录状态推定它们当时的全部步骤。

### 生产侧配置：已落实与尚需修正的部分

| 项目 | 本轮证据与判断 | 建议 |
| --- | --- | --- |
| 字段归属迁移 | `revise/application/config.py:672-712,825-876` 与 `inputs.py:80-88` 支持新 local 字段，旧 global 字段明确拒绝。当前 request 配置测试 37 passed | 承认配置归属迁移；真实执行尚未解除全局锚定对谱系过滤的依赖 |
| 过滤执行时序 | `reconstruct.py:133-137` 先 preprocess；`application/inputs.py:80-88` 传 include/exclude，`application/preprocess.py:125-149` 过滤 reference 与 spatial；`recon/pipeline.py:41-42` 随后 global anchoring → local refinement | 如果共识目标是全谱系 global、local 才筛选，必须改变过滤作用阶段并保护 Raw/全局输入；只改 YAML 归属不能验收这个目标 |
| 90/90 声明 | 没找到对应当前提交的日志；找到的 759 passed 属于旧 `a26d36f...`。当前本机测试经 capture/plugin/cache workaround 验证配置模块 37 项 | 更新为明确提交、命令、解释器与覆盖范围的实际测试结果；不称全工程或科学通过 |
| 旧 batch 字段 | `revise/batch/runner.py:127-132,579-583` 仍读 global include/exclude，Application 已拒绝这些字段 | 入口统一迁移政策：显式迁移并记录来源，或在任务选择前给明确错误；避免运行到一半才失败 |
| provenance 的有效配置身份 | `revise/application/publication.py:40-46,93-105` 没有登记新 include/exclude 与 local mapping 来源；子代理隔离验证：仅 include 列表不同，source SHA 不同，但 effective request/hash 相同 | 把所有影响参考范围和方法的字段纳入 effective request；否则结果身份不能区别不同实验 |
| mapping 多处冲突 | `config.py:709-712` 优先 mapping 而非 ist_mapping；`:768-775` output.ist_mapping 又可覆盖 local；根 batch 两处目前值一致 | 建立唯一归属，重复字段相同可明确归一，冲突必须显式处理并登记，不静默改变含义 |
| 生成器和模板 | `scripts/setup_batch_structure.py:30-40,85-103` 未同步 root 新语义；Application/Batch 维护模板没有演示新字段 | 同步所有仍用的生成入口，避免再次生成旧结构；字段省略合法，不把缺示例当运行缺陷 |
| 参考文档 | `docs/source/application-reference.rst:99-113` 未反映新字段与完整 mapping 方法，`:299-307` 又部分记载 within/outside | 更新当前文档；保留冻结旧配置与运行记录的原始身份 |

上述迁移测试属于当前配置解析的有界技术检查。宿主 Python 3.12 不是 README 声明的正式 Python 3.10/3.11 环境；普通 pytest 的 capture 崩溃不能单独归因于项目变更，也不能以 workaround 的 37 项通过代替受支持运行环境的完整回归验证。

子代理本次测试命令（在 REVISE 子仓）：`PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONDONTWRITEBYTECODE=1 PYTHONFAULTHANDLER=1 pytest -q tests/application/test_request.py -p no:capture -p no:cacheprovider`。环境解释器 `/Users/stephen/miniconda3/bin/python`，Python 3.12.4，输出 `37 passed in 0.34s`。未运行正式重建。

## 3. 科学表述需要调整

| 交接材料表述 | 建议的科学边界 |
| --- | --- |
| 恢复真实表达、突破仪器物理极限 | 核对 Fig6/7 既有 case 的实际支持范围，并区分原位测量与参考推断。强定量恢复主张需对应 case 的有效对照/真值证据；遮蔽、独立样本与正交观测是可能补证方式，不把全部新实验或外部队列设为当前完成门槛 |
| 更多 DEGs 等于检出灵敏度提高 | 先限定共同检测背景、比较对象、阈值与模型；灵敏度需要真阳性分母，同时报告假阳性。预测误差也应进入推断 |
| MUC2/KRT8 下降证明横向扩散污染已清除 | 需本样本转录本位置、分割归属或正交证据；重建后典型谱系 marker 下降也可能是参考压制真实异质性 |
| 共定位还原 TLS 因果网络 | Fig6 有既有共定位、TLS score 和候选 LR case；真实结果可复用，但共定位及 RNA 通讯推断不单独等同因果。TLS 成熟度若成为具体问题，再补所需结构观测，不作为当前故事交付的统一门槛 |
| CAF 收缩、溶骨、逆 Warburg 与侵袭已被揭示 | Fig7 既有程序/空间 case，论文另描述 LPL/MYH11 IFA 与 SMC→CAF 谱系追踪，不能忽略这些既有证据；它们分别支持定位、临床关联与谱系关系，不自动证明全部力学、代谢通量或侵袭机制。“ECM breakdown”也不应直接替换成骨吸收 |
| 代码硬指标 + 小载荷杜绝幻觉 | 代码和载荷也会传播错误事实；需源证据可追溯、缺失与限制传播、数值一致性和独立复核 |

外部依据：表达预测误差会影响下游假发现，[TISSUE](https://pubmed.ncbi.nlm.nih.gov/38347138/) 给出不确定性感知的推断框架；[LIANA+](https://www.nature.com/articles/s41556-024-01469-w) 将通讯推断定位于假说生成，RNA 与空间关联不直接等同功能互作；[TLS 原始研究](https://www.nature.com/articles/s41467-025-59341-w) 使用病理及多重染色确认结构和成熟度。这些是判断方法边界的依据，不能替代本样本验证。

### Fig6/7 既有 case：本次补核与范围修正

论文实际路径为 [Cell/main.tex](/Users/stephen/Documents/wuyushuai_project/wuyushuai_paper/REVISE_NBT/product/Latex/Cell/main.tex:419)。子代理初次回报曾多写 `wuyushuai_research_project` 路径段，主代理读取失败后重新定位并实际读取 415–520 行，已解决路径问题。

| case | 既有材料与主要 panel | 本次核查得到的范围 |
| --- | --- | --- |
| T | `REVISE/reproduce/case/Xenium_sc_SVC_T.ipynb`；Fig6 C–G，正文 419–422 | 论文 7 个 T 亚群，重点 T_0 增殖、T_2 TRM、T_3 Treg；marker/DEG/特异性与具体区域分布的故事已有基础 |
| TAM / T–TAM | `Xenium_sc_SVC_Monocyte.ipynb`；Fig6 H–J、M–N，正文 423、458–465 | 8 个 TAM 亚群，重点 TAM_1 招募、TAM_2 SPP1+、TAM_3 Mreg；既有 CellPhoneDB case 与空间共定位/TLS score，可以串联分子、通讯和区域 |
| CAF | `Xenium_sc_SVC_Fibroblast.ipynb`；Fig7 A–G，正文 472–486 | 10 个 CAF 亚群，重点 CAF_1/4/5/8；Hallmark/GO BP、CAF_4/LPL、CAF_8/ECM breakdown、CAF_5/SMC-like 的程序与空间故事已有基础 |
| 已有验证材料 | Fig6 K–L 及补图；Fig7 F、H–M，caption 427–455、492–519 | 论文描述 TCGA/独立队列、LPL 与 MYH11 IFA、Myh11 谱系追踪；外部队列是额外支持，Fig7 正交/谱系证据属于对应 case 的验证；本轮未独立复核临床/湿实验原始数据 |

因此，用户指出 Fig6/7 已有具体 case 验证是成立的。之前把近期闭环收窄成“单假说另建验证”不符合本轮目标；现有 case 可作为主线，新队列与新实验不成为统一完成门槛。当前审阅发现主要针对新增展示包的计算与来源，不应外推为已有 case 全部没有科学证据。

另有实质叙事冲突：原论文 T_3–TAM_3 强调空间情境下的 TLS-associated program 与有利 DFS 关联，候选 ANXA1–FPR3/LTB–LTBR；新增 H1 的“TLS 协同抑制轴”不能替代这一核心发现。调节性名称不能自动推出该空间生态位抑制或不良预后。应组织观察、情境依赖与候选机制，并保留不同证据支持范围。

论文 case 的 T 7 / TAM 8 / CAF 10 与新展示包 T 6 / TAM 8 / CAF 10 存在范围差异，标签编号也不能直接互换。复用时登记 cluster→论文标签→当前生物学实体的映射依据与输入身份；不能把旧 AA Fib-only 交付的 expression 限制套给全部 case。`spatial` 与 `expr` 两侧按实际表型/群体关系比较，不默认为同一批 unit 的逐一配对观测。

## 4. 三层架构与两个战略议题

### 保留三层阅读视图，明确计算与知识责任

| 层 | 建议职责 | 必要补充 |
| --- | --- | --- |
| 01 生信分析 | 当前样本/侧的表达、标签、程序和空间特征 | 全覆盖范围与不可评估部分；完整结果而非只保留 top5/top8；不要预设某个亚群必为焦点 |
| 02 重建影响 | 在 gene、Pathway、CCI、亚型与具体区域各层系统记录前后状态、变化和可信度，并组织空间解释 | 可比基因/群体/尺度、reference 来源、混杂和负面结果；改变与收益分别解释，不能以几个正向指标代替完整覆盖 |
| 03 科学发现 | 围绕既有 Fig6/7 case 串联亚型/基因的生物学故事，形成证据支持的解释与假说 | 所有关键陈述可追溯；支持/反证、替代解释、可检验预测、下一步及假说更新记录；外部队列追加证明 |

三个槽位是阅读入口，不应成为所有科学计算强制串行的依赖关系。沿用 `methods → analyses → runner/batch → index/result` 的已有责任划分，不为目录美观建立第二套科学流程。

已有方法复用入口：Gene 的 `methods/differential_expression.py`、`analyses/impact_molecular.py` 与 `analyses/spatial_autocorrelation.py`；Pathway 的 `methods/programs.py`、`analyses/pathway_activity.py`；CCI 的 `methods/cell_communication.py`；亚型/区域的 `methods/partition.py`、`methods/regions.py`、`analyses/reconstruction_impact.py`。CCI wrapper 本身不接收空间半径/区域，不等于完整空间通讯流程；区域模块已有 Anatomy/State/Gain 与物理窗口，不等于 TLS 成熟度或力学拓扑。逐 case 检查前提并补正式结果衔接，避免另建重复框架。

### 02 槽位：多层覆盖与组织矩阵

| 层面 | 需要覆盖的前后内容 | 空间连接与比较条件 |
| --- | --- | --- |
| ST unit / 细胞与亚型 | 覆盖范围、状态/亚型分辨、标签与组成变化、unit confidence；仅在需要时做明确的局部成员对应 | 单位代表 cell/bin/spot/virtual cell 的语义；亚型出现在哪里；不将所有 Raw/SVC 分析强制配对 |
| Gene | 检测与表达、亚型特异性/差异表达、关键 marker、共同实测基因与重建新增基因、gene confidence | 表达位于哪些 unit/区域；区分变化、新增预测与实际验证；保留完整检验表与可评估分母 |
| Pathway / program | 活性与富集、涉及基因及覆盖、亚型程序、Raw/SVC 差异及可信度 | 程序在哪些区域集中或渐变；评分背景、覆盖与尺度可比时才做定量差值，单侧证据仍保留 |
| CCI | sender/receiver、ligand/receptor、候选强度/检验、空间支持及前后变化 | 哪些区域有表达与邻域支持；表达、细胞身份与邻近关系的可信度分别说明；缺基因不补零 |
| 空间区域 / niche | 区域身份、细胞组成、局部程序/通讯网络、边界/邻域及可支持的连续变化 | 统一物理单位与实际坐标范围；明确区域来源；联系以上各层，不能只加一张总览空间图 |
| 综合故事 | 某亚型/基因对应什么程序、通讯与区域，重建增添或改变了哪些证据 | 对应 Fig6/7 case 的图/底表/命题；分清观察、解释、机制假说及仍待追加的证明 |

阅读入口可以同时支持“亚型/故事 → 分子与通讯 → 空间证据”和“区域 → 亚型组成 → Gene/Pathway/CCI”；两条路径引用同一组事实，避免为展示重复生成不同数字。完整覆盖指每层有明确结果或不可评估原因，不要求每个亚型机械地产出所有显著结果。

### Uncertainty：已明确方向与进一步建议

ST unit 与 gene 两轴继续作为主要刻画与展示维度。建议在此基础上区分测量支持和重建输出可靠性，并明确每个置信度究竟针对表达、类型/映射还是位置。两轴汇总不足以说明某 gene 在某个局部区域一定可靠；在科学问题需要时可查询条件性 `unit × gene` 可信度，不预先要求保存一张庞大全矩阵。

| 消费位置 | 建议消费的可靠性信息 |
| --- | --- |
| 重建 | 使用与算法估计对象相符的 unit/gene 信息来约束映射或表达推断；如何赋权/正则由后续方法方案决定 |
| Gene 与 DEG | 展示当前 gene/群体的测量和预测支持、比较的不确定性；避免只显示一个全局 gene 分数 |
| Pathway | 结合有效基因覆盖、贡献及其不确定性；不假定相关基因的误差彼此独立 |
| CCI | ligand/receptor 表达、sender/receiver 状态及邻域支持共同进入解释；表达可信不等于功能互作已验证 |
| 空间区域 | 区分实测/继承坐标与推断位置、区域归属与状态定位的可靠性；坐标未重建的路线不伪造位置重建置信度 |

校准可优先复用既有 case 中适当的对照/遮蔽/已知结构：检查声明更高可信度的结果是否有更小实际误差，置信区间是否符合预定覆盖。原始未测、低检测、参考缺失与预测不稳定的含义分别保留；不能简单把任意 uncertainty 分数变为 `1-u` 后称为正确概率。[TISSUE](https://pubmed.ncbi.nlm.nih.gov/38347138/) 提供预测区间与下游不确定性感知推断的参考，不作为必须迁入的算法。

本次子代理对已有实现的窄核查：

- `revise/backend/kernels/gene_uncertainty.py:89-106,161-199,218-237` 的生产调用见 `sc_svc_impute_benchmark.py:102-108`，输入是单细胞 reference，图来自表达 PCA/kNN 并按 cell type 分块，未使用物理 spatial 坐标。输出为 gene 的 Moran/Geary 与均值/方差；`test` 仅表示是否在 overlap_genes，不是显著性检验。它是表达流形结构稳定性的实验线索，尚不是 ST 测量误差、重建预测区间或 unit×gene uncertainty。
- `OTKernel._publish_annotation()`（`ot.py:454-474`）的 `Confidence` 为类别 posterior 的行最大值；熵与 top2 margin 也描述类别区分。当前没有校准过程，POT/TACCO 对 posterior 行归一化的处理还有差异，不直接解释为跨方法可比的正确概率。
- 当前坐标检查验证维度、有限性和对应一致，不提供配准、segmentation 或真实位置误差估计。现有 Confidence 不能代替位置置信度。

因此最新“两维 uncertainty”方向不能仅把旧指标改名接入。建议定义估计对象后再确定算法：生物学表达变异与技术/推断不确定性分开，避免把真实异质性当噪声；unit 轴上的类型、表达与位置可靠性分别解释。此处是方法建议，不要求用户放弃两轴，也不授权本轮实现新 uncertainty 方法。

### 空间连续性与高阶结构

战略议题成立，但应由具体假说决定模块。最小候选包括：有独立区域依据的肿瘤边界/血管距离场；多尺度、丰度与组织区域控制的邻域关联；具备所需细胞类型与病理依据的局部结构分析。距离场可放 01，Raw/SVC 差异放 02，结构机制假说放 03。

RNA 坐标不能单独给出真实血流灌注、缺氧物理场、胶原纤维方向和力学张力。缺少影像/标记/区域身份时列为不可评估，不通过增加图指标制造结构真实性。空间连续变化也不自动等于时间轨迹或因果过程。

### 队列扩展

本节作为未来/额外验证的能力约定保留，不作为近期完成 Fig6/7 case 与 Agent 生物学故事的验收门槛。工程上的多样本执行与失败隔离仍应复用既有 batch。

| 必要约定 | 最小内容 |
| --- | --- |
| 统计实验单位 | patient/donor、sample/section、region/cell 的关系；独立患者数单独报告 |
| 跨样本实体 | 稳定假说/细胞状态定义及注释依据；不能直接把不同样本的 `T_5` 当同一实体 |
| 可评估性 | panel、细胞/区域覆盖、QC、所需输入与执行状态；缺失、失败和不显著分开 |
| 每样本证据 | 同定义效应量及不确定性、方向、方法和比较口径；不能只收 p 值或显著票数 |
| 聚合 | 先分层核对可比性，保留异质性；P2/P3 显著、P1 不显著本身不证明亚型特异性 |
| 增量更新 | 固定纳入样本/run/结果身份的快照；新结果形成新版本，相关队列结论失效后再生成 |

已有逐 sample × analysis 的 batch 错误隔离应复用。动态资源调度、断点恢复、队列统计和 Agent 自动综述是不同能力，需分别声明已实现、计划和未验证，不能由 `summary.json` 文件存在推出全部具备。

## 5. 旧内容与当前内容的处理规则

用户补充：可能有以前内容过时而未更新或删除；本轮核对后列举即可。

只列清单并说明引用关系，不执行清理。冻结 receipt、provenance、旧运行与科学审阅保留其当时身份；当前入口与规范需要避免继续引用过时字段、旧目录或不存在的脚本。新代码中的直接赋值/随机证据属于当前实现问题，不能统一归为旧文档未更新。

| 入口或内容 | 本轮核查结果 | 后续建议（本轮不执行） |
| --- | --- | --- |
| `docs/specification/README.md:7-39` 与双域 overview | 仍定义 `00_benchmarks/01_discoveries`，`docs/analyses/README.md:3-14` 仍引用它作为 Agent 当前入口 | 确立当前规范与 legacy 说明，修正仍有效的阅读链接 |
| `scripts/build_dual_domain_slots.py` | 仍能生成旧双域目录 | 明确属于历史构建器，或在统一契约下保留过渡入口 |
| `scripts/run_full_subcluster_analysis.py`、`integrate_gold_standard_analysis.py` | 仍写旧双域路径 | 核对当前调用关系后标记用途，不能只删输出目录 |
| `scripts/run_unbiased_all_subcluster_analysis.py:274-386` | 新脚本依然写 `01_discoveries` | 当前执行路径未迁移完成；“旧目录已删”不等于旧架构已退出 |
| `scripts/organize_three_tier_analysis.py:39-48,158-162` | 依赖旧 `01_discoveries` 的 CSV；旧输入和目标表都不存在时无法重建 | 输入改为明确的登记结果；构建器应从既定事实来源重建阅读视图 |
| `docs/specification/03_agent_scientific_discovery_workflow.md:41-85` | 使用旧路径/字段，包含 `neighborhood_enrichment_matrix.csv` 等旧表入口 | 列出旧新字段对应，明确不可直接替换不同统计定义 |
| 新 `output/CRC/.../README.md` 与 `analysis/README.md` | 其中 `../../REVISE_Analysis_Agent/...` 从样本目录解析不到真实仓库规范 | 修复具体相对链接；重复展示树分别核验 |
| CRC 顶层 README 的 `aggregate_cohort_analysis.py` | 当前代码库没有该脚本，规范中自动队列聚合也没有 batch 接入 | 改为明确规划项，不呈现为现有可调用能力 |
| 新三槽位树与正式 API 协议 | 样本缺 `index.json`、各分析缺 `result.json`，没有正式每项状态；现有正式 runner/batch 仍独立有效 | 让阅读视图引用正式结果，避免两套互不一致的事实登记 |

`build_sample_analysis_delivery.py` 当前只建目录并写 `analysis.json`/`analysis/README.md`，不生成顶层报告与 summary，也不核验 CSV 或 hash；打印验证完成不构成验证。`AGENT_GENERATION_SPEC.md` 中自动触发 Agent 的流程尚为规范示意，未找到实际生成调用。以上是当前能力差距，应与旧内容过时分列。

此前正式记录 `p2-xenium-impact-20260922` 具有独立任务身份，说明当时 random SVC 参考表达转移、Moran 负向变化及 EMT 跨侧可比性限制；它不自动决定另一次新运行的科学结论。新增成果若要改变这些判断，必须提供新的输入、方法、运行身份与对照证据。

## 6. 推进建议与待讨论项

最新推进范围：以 Fig6/7 的具体 case 为近期验证与故事主线，复用可追溯真实历史成果，按 gene、Pathway、CCI、亚型、空间区域系统组织；补齐当前展示包中无真实依据的计算与来源缺口，将摘要和载荷接入同一事实链。只补影响这些 case 和目标故事的必要部分，不默认重跑全部历史结果或增加外部队列验证。

此前“先用一个假说验证再扩展”的建议已被用户的系统化 case 目标修正。Fig6/7 的 panel 与 notebook 已定位；后续需绑定可复用底表、输入与标签身份，以及细化 uncertainty 的估计对象与传播方式。事实修正服务于这组 case 和故事，队列作为额外证明。

### 给执行 Agent 的反馈候选（尚非共识任务书）

| 优先级 / 归属 | 反馈内容 | 可审阅验收方式 |
| --- | --- | --- |
| P0 / 共同输入身份 | 明确当前 focal analysis 来自哪些实际 H5AD/历史 CSV、哪些只是演示值；选定正式事实入口与派生展示入口。保护旧证据，不以最新目录名字替代 run 身份 | 输入与来源 manifest 能定位每个分析 scope、基因轴、标签、来源侧、reference、版本/参数和当前 hash；两个副本只作为副本登记 |
| P0 / REVISE | 修正 effective request/hash、旧 batch 入口与 mapping 冲突；核验字段迁移是否真正改变筛选作用阶段 | 影响计算的字段变化使有效身份变化；冲突行为显式；相同合法请求的 direct/batch 含义一致；给受支持解释器的有界检查记录 |
| P0 / Analysis 科学计算 | 取消正式事实里的手动分数、指定检验量和随机 Z；真实计算覆盖声明的可评估亚群。无法算的内容写不可评估，不补数字 | 每条关键数字能回到底表与真实方法/输入；声明测试数量、显著数量、完整结果与限制；known biology 不进入调分后再用于声称独立验证 |
| P0 / 报告与载荷 | 三层保留为阅读视图，引用正式 index/result；样本元数据与数值均从同次结果生成；演示、计算、科学审阅分别声明 | 摘要、payload、底表与 HTML 数值/分母一致；缺文件/partial/failed 可见；证据引用能解析；不存在固定四假说或固定 P2 数字被套给新样本 |
| P1 / 故事与假说 | 不固定 H1–H4；修正 H1 与原 case 的方向冲突，补问题、预测、反证、替代解释、下一检验；复用实际存在的正交/谱系证据 | 每项故事区分当前计算、已有 case 验证与候选解释；不把所有机制同时升级为因果验证，也不否定已有明确支持 |
| P1 / 既有 case 与完整故事 | 建立 Fig6/7 case→图/底表→验证命题→亚型/基因故事对应；按 gene、Pathway、CCI 和区域补齐必需证据 | 优先复用有效结果；区别旧 case 与新 run；每个层面有真实结果或不可评估原因，核心故事能读通、查证并解释重建带来的变化 |
| 方法讨论 / uncertainty | 保留 ST unit 与 gene 两维；明确测量支持、映射/状态、表达与位置的估计对象及下游消费 | 不能用旧实验指标冒充最新已确认算法；按适用 case 校准，并将方法定义、输出轴和含义交接给下游 |
| 后续 / 队列追加证明 | 队列实现可评估性、同定义效应与快照增量索引；外部验证作为额外证明 | 不阻塞近期 case 故事交付；未来 patient 独立性、失败/缺失/不显著与快照版本均明确 |

不应要求执行 Agent 一次实现全息表征、自动资源调度和所有 meta-analysis。先复用已有方法、runner/batch、发布保护和结果协议；扩展范围按本轮讨论确定。该顺序保留两大目标，同时避免扩大一个未验证展示包的影响范围。

## 7. 本轮核查登记

- task_id：`revise-architecture-review-20261010`
- owner：本轮审阅主代理；子代理分别核查生产端、分析科学计算、架构、qz 和外部科学边界。
- upstream/downstream：REVISE 三文件交付 → Analysis Agent 计算/报告 → 人与 Agent 科学审阅。
- execution_location：本地只读文件与代码审阅、隔离配置检查；qz 只读核查。正式业务结果与代码均未改写。
- 本地两份新样本目录：`REVISE_code/output/CRC/P2_CRC_Xenium` 与 `REVISE/output/CRC/P2_CRC_Xenium` 均为 31 个文件，相同相对路径的全部文件 SHA256 一致。
- Analysis Agent 当前正式能力：`reconstruction_impact`、`spatial_autocorrelation`、`pathway_activity`；新三槽位脚本未接入正式 runner。已有 batch 在 `revise_analysis/batch.py:38-54` 按 sample × analysis 隔离错误并写 `batch_result.json`；runner 已有临时计算与 `.previous/` 发布保护，应优先复用。
- 版本、远端事实、旧入口清单、运行身份、隔离测试命令与验证范围见本报告。
- 本轮主代理记录时间：`2026-10-10T11:28:01+08:00`；远端 curated hash 观测为 `2026-10-10T03:17:51Z`（上海 11:17:51），重建记录与回执路径核查为 `03:20:47Z`（11:20:47），目录计数为 `03:24:30Z`（11:24:30）。
- qz 实例：`cpu-download / 5cd61cb5-b0ec-40a6-ab66-874b1b941573`，本次从 RUNNING 清单解析；不固化为后续默认实例。
- 输入身份验证范围：qz 两个历史源输入当前 hash 与重建记录匹配；新三槽位科学表尚无已闭合的新 run 身份与计算链。本轮没有复算真实生信全流程。
- 验证命令类型：本地只读 `rg/nl/git` 与 Python SHA256 inventory；qz 通过 `qz list --json`、`qz run` 的只读 git/文件/输入 hash/进程核查。没有部署、重算或自动调用 Agent。
- 核心产物摘要：`summary.json` SHA256=`660a8e672b8d52957ff683a2a749bb932214745b042c226aced906a7d977d82e`；`agent_input_payload.json` SHA256=`3f596614c00059f359f4bdade419686aa7f512a82c3dbc9d9d39f293e2cf6cbb`。
- acceptance_level：仅本轮审阅发现与证据核查；不授予当前样本生物学接受。

本地取证：[31 文件当前清单](local-curated-manifest.json)。远端取证：[子代理根据 qz 返回值整理的观测索引](remote-observation-index.txt)，该索引不是原始 stdout，原始输出留在子代理工具调用记录中。主代理复核索引时发现 contamination_clearance 的 hash 有一处转录截短，原索引保留不改；已通过[单文件原始 JSON 补证](remote-single-file-supplement.json)解决，观测时间为 `2026-10-10T03:35:57.177823Z`（上海 11:35:57）。补证 size/hash 与当前本地一致，其余 30 项与本地完整 hash 相符。因此新展示包 31 文件的字节同步结论成立，完整目录、代码同步与科学验证分别判断。
