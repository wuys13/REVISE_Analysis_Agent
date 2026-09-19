# 当前证据、未完成工作与审阅边界

[审阅入口](README.md) · [协议与代码](contract-and-code.md)

## 核对对象与时效

本次文档核对日期：**2026-09-19**。本轮只读检查代码和保存证据，没有重跑科学分析或全量测试。以下通过数均引用已有运行记录。

| 仓库 | 分支与 HEAD | 本轮编辑前工作树 | 审阅含义 |
| --- | --- | --- | --- |
| REVISE | `revise-2.0`，`a26d36fe3b11d652f2cd9eeafb92e990518eb3ad` | 多处 tracked 修改及 untracked 交付/比较文件 | 本页指当前工作树，不能仅 checkout HEAD 复现 |
| REVISE_Analysis_Agent | `sl`，`291cf4f06aad187af0ccbf5bb465bf3b03f4ab88` | 仅 `other_reviews/` 未跟踪 | 本轮新增/修订文档；运行代码与数据不变 |

上游[联合证据索引](../../../REVISE/docs/development/reconstruction-analysis/verification-review-handoff.json)保存双方文件摘要与执行信息。本轮写文档前比对其中 **119 个 producer 文件、41 个 consumer 文件，均无摘要差异**。文档编辑会有意改变其中 consumer README 等文档摘要，不表示计算源码发生变化；后续审阅应区分文档差异与计算差异，并重新检查当时工作树。这里不复制或另建全文件摘要清单。

分析库交互验收记录的 source HEAD 仍是提交前的 `e028806` 加未提交实现，不应误读成当前 HEAD 的自动验收。本轮比对该记录所列的 12 个运行代码/Notebook 文件摘要，与当前文件一致；这仅确认被记录文件的对应关系，不扩展为所有文件或新一轮测试通过。

计划阶段看到的 `in_progress` 已被上游更新为 **`engineering_verified`**，关联最终 75 项检查及产物；因此“联合索引未闭合”不再列为本次当前缺口。这不扩大该验收的科学范围。

## 已有证据分别证明什么

| 对象与层次 | 已保存结果 | 能支持的结论／不能支持的结论 |
| --- | --- | --- |
| 分析库软件行为 | [交互验收记录](../../output/verification/interactive_acceptance.json)：103 项测试；最后局部绘图调整另有 14 项聚焦检查 | 支持输入、参数、重跑和展示行为；两批数量不累加，不等于全部通用方法真实验证 |
| P2 标签空间链 | 同一记录：17 个 Notebook 代码单元执行且源码一致，Notebook/batch 23 份表/JSON 一致，原 19 份科学产物未变 | 支持当前临时 SVC 主标签、坐标、State/Anatomy 及派生比例链；不证明其表达可用于生信对照 |
| 页面检查 | [浏览器记录](../../output/verification/browser_reading_qa.json)、[迁移记录](../migration.md)：桌面与窄屏目录、图片、锚点及键盘检查 | 界面技术证据；不是用户对阅读效果或区域定义的认可 |
| 上游最终联合 fixture | [最终运行摘要](../../../REVISE/output/review-handoff/20260919-final/summary.json)：75 tests、0 failures/errors/skipped，Notebook passed，4 PNG；[联合索引](../../../REVISE/docs/development/reconstruction-analysis/verification-review-handoff.json)已关联 | 标签空间消费者为 partial，已知线性 fixture 原生 Moran succeeded；覆盖生产交付与正式消费者入口，不是真实全量重建 |
| 四方法比较 fixture | [执行 Notebook](../../../REVISE/output/review-handoff/20260919-final/comparison.executed.ipynb)、[覆盖表](../../../REVISE/output/review-handoff/20260919-final/coverage.csv)、[指标](../../../REVISE/output/review-handoff/20260919-final/metrics.csv) | 验证合成 T 范围及实际运行的 resolution 0.6/0.7；不能宣称真实 T/Macro/CAF 或所有配置已比较 |
| 通用方法迁移 | [迁移记录](../migration.md)与[通用回归测试](../../tests/test_generic_parity.py) | 源码/行为迁移证据；CCI、trajectory 等可选资源与真实应用仍需各自核验 |

### 历次记录不能混用

上游 `20260919-a` 为 58 passed、1 次 kernel socket 权限失败；`20260919-b` 为 59 passed，但执行版未嵌入 PNG；`20260919-c` 为 59 passed、4 PNG，发生在最终 sST 修正之前。当前引用 `20260919-final` 的 75 项结果，不将旧批次数量相加，也不继续把已解决的环境问题称为当前阻塞。详细过程由上游联合索引维护。

分析库迁移记录中的 39/78 项测试、旧 Notebook 单元数和曾未完成的浏览器检查属于各自历史批次。最近交互验收见本页所链记录；文档整理没有产生新的软件或科学通过项。

## 自我审视与剩余工作

| 分类／事项 | 责任方 | 当前状态与证据 | 下一步完成条件 |
| --- | --- | --- | --- |
| 文档定位与扩展说明 | 消费端 | 本轮补充三个正式流程之外的独立探索路径；修正 Impact 总体“必须表达”的误导表述 | 入口、能力索引与 Notebook 说明一致；不把通用工具写成正式已验收流程 |
| 上游 scale 迁移旧描述 | 消费端文档 | 当前 `consumer_declaration` 已知线性去 scale，旧迁移段落已标历史 | 引用当前实现，保留 unknown 限制；不反推历史数据 |
| 真正生产到消费的真实数据链 | 生产端＋消费端 | fixture 已验证；P2 全类型真实交付及 sST 真实联合验收仍未完成，见上游 `remaining` | 有明确来源的真实三文件交付，核对覆盖、标签、坐标与表达，再按实际能力运行并保存结果 |
| P2 表达分析 | 生产端确认来源；消费端执行 | [样本声明](../../data/P2CRC_Xenium/sample.yaml)仍为 unknown；临时空间载体仅验标签空间链 | 正式交付或足够来源证据后更新配置；Raw 确认可推进 baseline/Gain/K-control/membership，Moran/AUCell 按侧验收 |
| P1 历史输入 | 生产端＋研究者 | [迁移记录](../migration.md)记录处理历史未确认，未宣称已分析 | 确认表达来源和主标签，或重新交付，再确定实际运行范围 |
| 四方法真实比较 | 当前由上游实验负责；用户审阅 | 三个历史标签载体来源已确认；整合基准和生成四种真实表达仍待完成 | 按[比较协议](../../../REVISE/docs/development/reconstruction-analysis/assembly-comparison-contract.md)准备输入，运行实际 scope/参数，检查覆盖与结果；不自动决定默认方法 |
| 扩展方法成熟度 | 消费端＋具体分析调用方 | 富集、CCI、轨迹等已有工具，但不是全部具有正式流程/真实验收 | 新任务按实际资源、前提与科学问题验证；成熟后再接入稳定流程，不为目录完整而开发 |
| old P2 parity | 消费端＋历史结果提供方 | 同口径历史窗口/结果不足，仍未完成 | 找到可比输入、定义、随机过程与结果，逐项区分保持、适配、方法修正；缺证据不填通过 |
| 科学与阅读判断 | 用户／研究者 | State/Anatomy、区域选择、assembly 优劣与阅读认可均不是软件测试结论 | 基于真实结果审阅；State/Gain 不以更大区域、更大差值作为成功标准 |

没有从本轮文档核对推导出必须新增的科学方法或通用框架。更深的代码缺陷检查留给后续联合审阅；本表不是“代码已无缺陷”的证明。

## 证据携带与交接

本次只保存轻量说明，未复制 H5AD、运行结果或日志。两库的 `output/` 可能不随 Git 分发；审阅者应同时获取相应产物目录，或使用[上游复现脚本](../../../REVISE/scripts/verify_review_handoff.py)在独立新目录运行。依赖、权限与资源不满足时记录未核验，不覆盖现有结果。

本轮不改上游。上游后续负责真实交付及其验收记录；分析端负责真实消费、独立扩展与报告解释。两侧后续变更应重新核对受影响协议与证据，不要求无关全量重跑。

## 本轮文档交付检查

两项独立只读复核分别检查职责/扩展意图，以及协议/代码索引/证据范围，均未发现需要返工的实质问题。主代理检查本次涉及文档的 100 个本地链接，全部可定位；`git diff --check` 通过。运行代码、Notebook 源码、输入数据及上游源码未由本轮修改；没有新增科学通过项。
