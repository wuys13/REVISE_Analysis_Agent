# 当前证据与剩余事项：2026-09-20

[协同入口](README.md) · [协议与代码](contract-and-code.md) · [上游完整验收及服务器命令](../../../REVISE/docs/development/reconstruction-analysis/acceptance.md) · [联合机器记录](../../../REVISE/docs/development/reconstruction-analysis/verification-mini-2026-09-20.json)

本轮实际修改两个工作树，并执行约1%连续空间ROI真实联测。生产端定义SVC，消费端按正式sample.yaml加载；没有用手写临时YAML替代正式发布。初始HEAD：REVISE `5b21cec38106b2d864b7d6d9ec8fb5cb87cf893a`；本仓库 `4118061983f15d0a51b6fdac8c9e4f534cd2c814`。验收包含未提交修改，用户other_reviews材料保留。

## 修改与复用

| 本轮任务 | 实际落实 | 证据 |
|---|---|---|
| 1–3 名称与scope | 删除label_aliases；cell type统一 /→_并保留NA；scope去重；parent-window/Anatomy参数同规则；主cluster不改 | Sample.labels、effective_parameters及聚焦回归 |
| 4 Raw Level2 | scope内有效ID及同ID坐标；保留原label-count summary，新增coverage并显示有效/排除数 | P2 All有效501/3409；Fibroblast112/112、Mono_Macro191/191、T183/183；部分/全缺失回归 |
| 5 K-control | 缺主标签、目标scope或窗口unavailable；真实计算异常保留stage_errors | hST/sST mini开启K-control，无阶段异常；sST原因是缺SVC_cluster，hST先受unknown表达限制 |
| 6 项目配置 | mini/full六配置直接引用上游正式sample.yaml；P2窗口/抽样/EMT参数迁入项目 | configs/p2_project.yaml、p1_hd_project.yaml、p2_visium_project.yaml及对应full文件 |
| 7–8 真实执行与修复 | 复用参数解析、17代码单元连续Notebook、三个workflow batch、保存结果报告 | 三份Notebook与batch有效参数一致；科学表/JSON逐字节相同；报告仅读保存结果 |

## 本轮真实状态

| 样本 | Raw → SVC | Impact | 独立Moran / pathway | Notebook/batch一致科学文件 |
|---|---|---|---|---:|
| P2CRC_Xenium_mini | 3409 → 501 | partial，stage_errors=0 | partial / succeeded | 170 |
| P1CRC_HD_mini | 5077 → 4128 | partial，stage_errors=0 | skipped / skipped | 11 |
| P2CRC_Visium_mini | 43 spot → 567 | partial，stage_errors=0 | partial / succeeded | 44 |

[iST报告](../../output/mini-acceptance/20260920/iST/notebook/report.html) · [hST报告](../../output/mini-acceptance/20260920/hST/notebook/report.html) · [sST报告](../../output/mini-acceptance/20260920/sST/notebook/report.html)。对应目录保留执行版Notebook、有效参数、result.json和表图。[逐文件parity记录](../../../REVISE/output/mini-acceptance/20260920/evidence/notebook-batch-parity.json)和[12份报告只读刷新](../../../REVISE/output/mini-acceptance/20260920/evidence/report-refresh.json)位于上游证据目录。

P2小ROI缺稳定State/Gain阈值；QC后Raw推断broad存在NA使Anatomy不可用；常量基因Moran明确排除。hST输入表达身份未确认，因此Raw/SVC均unknown，表达分析未验收。hST/sST正常无SVC_cluster，仅相关能力受限；sST不开启同单位membership。Visium0.73 µm/coordinate暂定。

sST生产校正实际有821个parent–gene有效零项，占总目标表达量约0.774%，单parent最大差约19.53%；其余138636项在容差内守恒。独立分子分析执行成功不等于此校正问题已解决。原算法保持，策略讨论归生产端。

消费者全套回归在NUMBA_DISABLE_JIT=1、单线程隔离环境为109 passed；最后report聚焦31 passed，彼此重叠，不累加。上游正式发布器到本仓loader/表达分析的联合检查及mini工具检查12 passed，另有sST校验3 passed。普通环境出现过numba/OpenMP启动错误，成功隔离运行才是有效证据。历史103/78/75等记录仅作旧快照。

## 上游四方法比较

历史三carrier已汇集为49279×0；四方法来自同一mini ROI/reference/QC/GA/LR/seed，仅assembly/output改变。比较T、Mono_Macro、Fibroblast共同ID160/157/68及共同13088基因，逐ID原坐标严格校验。保存36指标行、36列联表和12图；[执行版Notebook](../../../REVISE/output/mini-acceptance/20260920/assembly/comparison/assembly_comparison.executed.ipynb)及[协议](../../../REVISE/docs/development/reconstruction-analysis/assembly-comparison-contract.md)。本仓只链接阅读，不复制算法，也不将其共同ID规则扩展为普通loader要求。

## 剩余事项

- 确认P1 HD源表达历史后再决定声明和重跑表达分支，未知不靠数值外观升级。
- 由生产侧讨论sST零/极小支持的校正策略；当前结果保持partial数值验收。
- 服务器执行全量、评估现有OT容量限制；mini成功不是全量容量通过。
- 用户审阅科学解释与默认assembly；State/Gain和ARI/NMI均不自动证明生物学改善。

运行 `python -m revise_analysis.cli batch --config configs/p2_project.yaml` 消费mini；服务器换 `configs/p2_full_project.yaml`（其他路线同理）。CLI对partial/skipped返回1，必须读取result而不能当作统一计算失败。完整命令在上游验收入口。output内产物被Git忽略，需随运行数据保留；缺文件时明确未核验。
