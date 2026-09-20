# 当前证据与剩余事项：2026-09-20 收口

[协同入口](README.md) · [协议与代码](contract-and-code.md) · [上游四表验收及服务器命令](../../../REVISE/docs/development/reconstruction-analysis/acceptance.md) · [联合机器记录（当前见 closeout）](../../../REVISE/docs/development/reconstruction-analysis/verification-mini-2026-09-20.json)

正式主链保持 REVISE 发布三文件 → 本仓库独立加载 → Notebook/batch/report 共用流程。三类路线依据字段、表达身份和坐标判断能力，没有新增平台分支。当前基线：REVISE `464ee26ac22978c3994e40af50325413e70cef13`，本仓库 `3e70ab364f89f9e24e981e3e9239e7fe563138a6`；本轮修改已提交并推送：REVISE `467af8b9a104328fba8db0e786471620eab4150f`，本仓库 `4c4278b3adde1907fd66474ff97957a9c7dab790`。

## 本轮实际修改与简化

| 对象 | 修改 | 边界 |
|---|---|---|
| Anatomy | 使用完整交付 SVC broad 与坐标，产物改为svc_anatomy_*；共享Raw原点保持 | Other只描述SVC窗口未观察到Tumor/Normal；无对应格为Unknown；不把Raw原始标签重新注释 |
| State/Gain和Raw对比 | 复用原窗口、scope、阈值及独立对象关系 | 本轮iST的44份State/Gain、diversity、support表与旧结果逐字节一致 |
| Notebook | 默认正式P2 project；保留显式sample/project/override；按字段描述能力 | 不再默认旧临时carrier；单独sample不套用默认project |
| 项目配置与报告 | 三项目输出写入closeout新目录；sST指向新正式交付；报告读取新命名和准确来源 | 旧产物保留，报告不重算科学结果；原始Raw历史来源翻译不改成SVC |
| 上游sST | stable share-then-target；正式H5AD保存校正前支持诊断 | 消费者不补表达，生产侧真零支持仍unresolved |

旧阶段已完成的cell type /→_、NA保留、label_aliases删除、Raw Level2 scope有效ID、K-control局部unavailable及正式项目配置继续复用。没有新增注册／调度／兼容框架。

## 真实运行与验收

本轮仅重跑sST生产重建；iST/hST使用上一阶段正式交付，三路分析重新执行。

| 样本 | Raw → SVC | SVC Anatomy／Raw Unknown | Impact | 独立Moran／pathway | Notebook/batch一致科学文件 |
|---|---|---|---|---|---:|
| P2CRC_Xenium_mini | 3409 → 501 | 96窗均Other；401 Raw点Unknown | partial | partial／succeeded | 172 |
| P1CRC_HD_mini | 5077 → 4128 | 214窗；275 Raw点Unknown | partial | skipped／skipped | 14 |
| P2CRC_Visium_mini | 43 spot → 567 | 1个Interface窗；0 Raw点Unknown | partial | partial／succeeded | 45 |

三份Notebook各17个代码单元、无cell error；全部结果无stage_errors。三路SVC Anatomy窗口/context均生成；sST的Anatomy stage因部分Raw parent scope缺失为partial。hST/sST无SVC_cluster不再阻断Anatomy，仍限制需要主cluster的分支。hST表达身份unknown，表达计算继续跳过。sST单窗结果依赖小ROI和暂定尺度，不作全量组织结论。

[iST Notebook](../../output/mini-acceptance/20260920-closeout/iST/notebook/01_reconstruction_impact.executed.ipynb) · [iST报告](../../output/mini-acceptance/20260920-closeout/iST/notebook/report.html) · [hST报告](../../output/mini-acceptance/20260920-closeout/hST/notebook/report.html) · [sST报告](../../output/mini-acceptance/20260920-closeout/sST/notebook/report.html)。[逐文件parity](../../../REVISE/output/mini-acceptance/20260920-closeout/evidence/notebook-batch-parity.json)与[12份报告只读刷新](../../../REVISE/output/mini-acceptance/20260920-closeout/evidence/report-refresh.json)记录有效参数一致、表JSON一致和保存结果／表／图哈希不变。

| 检查 | 结果 | 证明范围 |
|---|---|---|
| 消费者四文件聚焦回归 | 53 passed | SVC Anatomy、无cluster、Unknown、网格、报告和交互失效；[日志](../../../REVISE/output/mini-acceptance/20260920-closeout/logs/tests-consumer-focused-53passed.log) |
| Notebook入口实际5场景 | passed | 默认project、显式sample/project、override及错误归属拒绝；[记录](../../../REVISE/output/mini-acceptance/20260920-closeout/evidence/notebook-entry.json) |
| 上游正式发布到真实消费者 | 2 passed | 三文件加载、SVC背景与Raw未覆盖点、表达及输入保护 |
| 上游数值与verifier | 28 passed | 正支持／真零、诊断序列化、旧无诊断unobserved |
| 输入与旧证据 | 25个输入／旧交付和71份旧产物哈希不变 | 不覆盖旧文件和用户审阅；不代表科学或全量验收 |

本轮工程环境使用固定线程与必要隔离；初期一次pytest capture导入崩溃不计为通过。旧109/31/12/3等数值为上一阶段记录，聚焦子集与完整套件不累加。

## sST 数值结果与来源边界

新产物中139207个正支持parent–gene项全部在容差内守恒，最大相对误差约1.07e-7。旧821个有效零项中571个恢复；本轮诊断直接识别剩余250个校正前真零，涉及2个parent，目标缺口约0.253%，最严重parent约9.90%。[正式交付验收](../../../REVISE/output/mini-acceptance/20260920-closeout/evidence/delivery-P2CRC_Visium-default.json)保持partial。

这些目标是共同基因轴上的内部parent归一化表达，不是原始counts。Moran/pathway执行成功不证明真零支持已经解决，也不证明生物学改善。旧输出没有校正前诊断，不能将其最终effective-zero直接命名为true zero。

## 上游比较与剩余事项

H2历史baseline 49279×0及四assembly、H4共同ID/基因/坐标检查的代码与产物未变，本轮复用。[比较Notebook](../../../REVISE/output/mini-acceptance/20260920/assembly/comparison/assembly_comparison.executed.ipynb)保存36指标、36列联表和12图。比较类型T／Mono_Macro／Fibroblast共同ID160／157／68、共同13088基因；不把该局部规则变成普通loader要求，不自动选择assembly。

| 剩余问题 | 当前状态与后续证据 |
|---|---|
| sST真零支持 | 不自动补值；需要生产侧明确新分配语义及跨样本证据 |
| hST表达来源 | 保持unknown，待可追溯预处理记录 |
| Visium物理尺度 | 0.73 µm/coordinate仍暂定，待样本标定来源 |
| 全量容量与科学解释 | 服务器运行、实际峰值、State/Gain与跨平台解释均未验收 |
| 共坐标与覆盖差异 | 保持真实parent和单位关系；mini不证明空间统计无偏 |

默认Notebook与`python -m revise_analysis.cli batch --config configs/p2_project.yaml`使用同一正式mini项目；服务器用对应full配置。CLI对partial/skipped返回1，应读取result的具体状态与stage_errors。当前产物在`output/mini-acceptance/20260920-closeout`，旧`20260920`目录为前一阶段快照；二者均须随运行数据保留。
