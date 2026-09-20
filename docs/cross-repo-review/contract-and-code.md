# 交付协议与关键代码索引

[审阅入口](README.md) · [证据与缺口](evidence-and-gaps.md)

本页是双向实现导航。消费者的正式调用定义仍在[输入输出协议](../input-output.md)，生产侧见[交付说明](../../../REVISE/docs/development/reconstruction-analysis/plans/analysis-handoff.md)。核对日期为 2026-09-20；行号只作该工作树定位提示，函数名优先。

## 生产与消费对应关系

| 交接项 | REVISE 生产行为 | Analysis Agent 消费行为 | 证据与限制 |
| --- | --- | --- | --- |
| 文件、路径与身份 | `sample_document` 生成 schema v1、同目录 raw/SVC 文件名；sample ID 作可逆 percent encoding，保留原 ID | `read_sample_config`、`resolve_path` 按 YAML 位置解析；`safe_segment` 检查目录段 | 上游 `test_analysis_delivery.py` 覆盖真实消费入口与相对路径；不是接受所有任意路径身份 |
| Raw 与 SVC | `reconstruct` 保留/重读原 Raw；`prepare_raw` 按实际 ID 回填推断别名，保留原矩阵、原标签及冲突信息；发布 pipeline 的 SVC | `load_sample` 读取两个对象，不用 Raw 表达覆盖 SVC，不写回 H5AD | 原始侧矩阵不一定等于整数测序 counts；回填目标为上游保留的推断别名，不保证重复回填任意历史别名无冲突 |
| 表达接口 | 已确认线性来源经数值检查后输出 matrix=X、identity；未知来源保留限制 | 固定 finite/nonnegative/unlogged linear；允许小数及小值；旧 log 声明拒绝，identity unknown 或旧 scale unknown 阻止表达消费者 | `expression.consumer_declaration` 对已知来源不写旧 scale；本轮P2输入由用户确认counts；P1 HD仍unknown。loader 不猜历史，数值检查在表达方法执行 |
| 标签 | 交付映射指向推断 broad/subtype；SVC 实际有 `SVC_cluster` 才声明 reconstruction | 按映射读列，cell type仅统一 /→_、保留NA；State使用原主标签；Raw Level2按当前scope有效ID计算 | 主重建分群不等于已注释 biological subtype；可读取不代表每条 route 都具备 State 所需标签 |
| 单位及基因轴 | 保留真实 ID、坐标及 route 对应关系，不强造通用 Raw/SVC 一一配对 | ID/基因名须唯一，两侧数量、顺序、基因集合可不同；membership 使用明确共同 ID | sST parent 信息不是共享同单位 ID；不能按行号或坐标推断成员身份 |
| 空间坐标 | 检查两侧声明坐标形状及有限性，um 转 micron；未知比例不猜 | 物理窗口按声明单位与转换比例判断；Anatomy 与 State 独立尺度 | sST 当前发布虚拟单位坐标；本轮sST实际parent/坐标通过；0.73尺度暂定 |
| whole-sample | iST 按实际 broad type 分析；不符合条件的类型带原因跳过；eligible 计算错误使样本失败，不发布半样本成功 | 分析按实际交付对象与 scope 判断，不假定 SVC 覆盖所有 Raw 单位 | “whole-sample 交付”不等于所有类型都成功重建；查看上游覆盖与跳过记录 |
| 参数 | 生产参数与来源保留在上游交付信息 | 默认值 → sample → project → override；资源路径按声明文件解析；Notebook/batch 共用解析 | 重建参数与分析参数是两套职责，不要求值或随机过程相同 |
| 失败与发布 | 同次三文件发布、备份与回滚 | runner 发布当前分析结果；缺前提与执行错误分开；Impact 失效阶段撤销旧登记 | 上游样本失败与下游分析 partial 语义不同；旧磁盘文件不自动成为本次证据 |
| 报告与扩展 | 提供生产与比较证据 | `index.json → result.json → tables/figures/report.html`；HTML 仅读保存产物 | 通用方法调用不自动拥有正式流程的 manifest、状态和报告保障 |

Raw 表达确认可推进 Raw baseline、Gain、K-control 和 membership，但还需各自的标签/坐标/共同 ID。Moran/AUCell 按侧判断，不能因另一侧缺失而删除已有原生结果或补零。完整科学解释见相应方法页。

## 按调用链找代码

链接是实际文件；符号名用于搜索，括号中的行号为本次核对位置。`REVISE/` 指相邻生产仓库，`revise_analysis/` 指本仓库包。

| 阶段 | 生产端入口 | 消费端入口及验证 |
| --- | --- | --- |
| 配置与来源 | [application/config.py](../../../REVISE/revise/application/config.py)；[expression.py](../../../REVISE/revise/application/expression.py)：`parse_expression_declaration`、`bind_expression_sources`、`consumer_declaration`（15/48/63） | [io.py](../../revise_analysis/io.py)：`_expression_declaration`（13）、`Sample.expression_unavailable`；[线性接口测试](../../tests/test_linear_contract.py) |
| 重建与交付 | [reconstruct.py](../../../REVISE/reconstruct.py)：`reconstruct`（28）、`run_application`（99）；[delivery.py](../../../REVISE/revise/application/delivery.py)：`prepare_raw`（56）、`sample_document`（91）；[publication.py](../../../REVISE/revise/application/publication.py)：`publish_outputs`（195） | [io.py](../../revise_analysis/io.py)：`read_sample_config`（129）、`load_sample`（145）；上游[联合消费测试](../../../REVISE/tests/integration/test_analysis_delivery.py) |
| whole-sample 与 sST | [adapters.py](../../../REVISE/revise/backend/adapters.py)：whole-sample iST 分支；[sST runner](../../../REVISE/revise/backend/runners/sc_svc_super_resolution_application.py)：最终 SVC 表达与坐标；[生产交付测试](../../../REVISE/tests/application/test_sample_delivery.py) | [Impact](../../revise_analysis/analyses/reconstruction_impact.py)：`ImpactWorkflow.stage_input` 及各 scope；缺标签/表达前提按能力记录 |
| 方法与正式流程 | 输出事实是分析输入，不调用消费者选择算法 | [partition.py](../../revise_analysis/methods/partition.py)、[regions.py](../../revise_analysis/methods/regions.py)、[spatial.py](../../revise_analysis/methods/spatial.py)、[programs.py](../../revise_analysis/methods/programs.py)；三个正式流程及通用工具见[能力索引](../analyses/README.md) |
| 参数与交互 | 重建配置仍由生产侧负责 | [runner.py](../../revise_analysis/runner.py)：`ANALYSES`、`resolve_analysis_parameters`（93）、`run_analysis`（169）；[batch.py](../../revise_analysis/batch.py)：`run_batch`（8）；[Impact Notebook](../../notebooks/01_reconstruction_impact.ipynb)；[参数测试](../../tests/test_project_parameters.py)、[重跑测试](../../tests/test_interactive_workflow.py) |
| 保存与阅读 | [上游验收脚本](../../../REVISE/scripts/verify_review_handoff.py)收集 fixture 联测和比较产物 | [ImpactWorkflow](../../revise_analysis/analyses/reconstruction_impact.py)：`run_stage`、`apply_parameters`、`result`；[保存表绘图](../../revise_analysis/plotting/impact_figures.py)：`render_impact_figures`；[报告](../../revise_analysis/reporting/report.py)：`render_report`；[报告测试](../../tests/test_reporting.py) |

## 四方法比较是独立的研究问题

入口：[Notebook](../../../REVISE/reproduce/case/assembly_comparison.ipynb) → [assembly_comparison.py](../../../REVISE/revise/analysis/assembly_comparison.py) 的 `load_assembly_inputs`、`compare_assembly_methods` → [对应测试](../../../REVISE/tests/analysis/test_assembly_comparison.py)。真实输入与历史路径以[比较协议](../../../REVISE/docs/development/reconstruction-analysis/assembly-comparison-contract.md)为准。

四种方法在共同范围内各自进行表达预处理及 Leiden；旧空间标签不用于生成新表达分群。保留多 resolution、覆盖表、列联表及 ARI/NMI，不自动挑赢家。三个历史carrier已汇集为49279×0 baseline，四种真实mini产物及比较已完成，坐标按共同ID严格校验。本轮不复制该实验，也不将其局部配对规则写入普通 loader。
