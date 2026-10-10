# 第三轮只读复核：6b158ed 与当前 P2 展示包

结论：本轮有实质整改，载荷完整提取、隔离再生成、本地双目录和启智展示资产同步可以有界验收；“全部整改完成、所有虚构指标退出、展示事实100%一致、远端代码部署完成”仍不成立。当前整体为 **partial / review_complete_with_limits**，不是整包事实验收通过。

继续执行用户已确认的并行交付：合作方开发正式接口与已核查生信 case；执行方完成下面的收口。体系止于系统化生信分析与基本确立的生物学假设，临床、IFA、谱系和外部队列为体系外后期关联验证。本轮不要求整合这些验证，也不要求重算全部历史 case 或先完成 uncertainty 算法。

## 1. 对本轮汇总的逐项判断

| 声明 | 当前核查结果 | 验收范围 |
| --- | --- | --- |
| prepare_agent_payload 去除旧常数与 CCI 截断 | 成立：旧 Moran/20.7x/62.1% 字段撤出，完整读取当前表格 | 完整提取通过；样本 metadata 仍固定、来源和缺失状态未完备 |
| 43.4KB payload 与底表一致 | 43,423 bytes；隔离再生成逐字节一致 | 证明忠实复制当前表，不能证明表中叙述正确 |
| 14条 CellPhoneDB 统计显著条目 | 当前14条均 p<0.05，含2条 score=0；原来不显著的 LTB 条目已撤出 | 表达侧通讯候选；score、方向与显著性分别解释，同群候选不等于已证自分泌 |
| 40条 ORA | 8群×5条，27条 padj<0.05，13条不显著 | 不是40条均显著，也不是25群全覆盖；ORA名称已修正 |
| summary / report / README 已全部清理 | 局部成立，数字与科学叙述仍冲突；槽位文件和假说仍旧 | 整包同源和事实一致性未通过 |
| 所有人工构造指标退出 | 不成立：当前下钻入口及可执行旧脚本仍保留固定/随机科学量 | 需隔离旧入口并撤出当前事实层；不删除冻结历史证据 |
| 本地双目录完全同步 | 33/33文件 size+SHA256 一致，无路径差异 | 展示资产技术同步通过 |
| 启智资产部署完成 | 33/33 curated文件与本地两目录一致 | 展示资产技术同步通过，非新代码部署 |
| origin/sl 已推送6b158ed | 本轮启智 live ls-remote 确认成立 | 本地 tracking 与 live origin 都指向该提交 |
| 启智代码部署完成 | 不成立：AA checkout仍 e33d122，缺新 payload 脚本；对象库无6b158ed | 不可在该 checkout 声称能再生成本轮包 |

当前导出范围：DEG135行=9焦点群×top15；ORA40行=8群×5；CCI14条为预选 T/TAM case；空间5条为论文描述性记录；焦点故事8群。25是声明的亚群 universe，不能等同于每层实际评估范围。新增 focal DEG 文件与旧 full 命名文件同时存在。

本轮独立数值核对：DEG135/135与真实历史源行一致（logFC按导出三位小数容差）；ORA40/40与原始enrich输出一致；CCI14/14的pair、cell_pair、p、mean、score一致。CCI包含13个p=0及1个p=0.005；这是1000次置换结果的原始表示，p=0不解释为真实概率严格为零。真实替换应明确承认，不再重复第一、二轮已解决的数值问题。

但CSV中的 `Monocyte/20260814T130000Z-display/...` 等source别名不能从当前包直接解析；实际源文件位于并列 `REVISE/output/notebook_analysis_acceptance/{T,Monocyte,Fibroblast}/...`。登记实际base path/映射、源文件hash和选择规则即可，不需要重跑科学计算来补来源。

## 2. 仍直接影响 Agent 使用的证据

以下路径相对当前 `output/CRC/P2_CRC_Xenium/`；两份本地样本目录及启智 curated 文件相同。

| 问题 | 当前具体位置与事实 | 使用影响 / 最小收口 |
| --- | --- | --- |
| 假 Z 仍在活动阅读链 | H1正文9/17行4.85；H2正文4/7行4.35；H3正文旧Z；H4正文19行3.92；根 report 仍链接这些文档 | 不是仅有历史文件残留。撤出当前计算证据；描述性空间记录保持描述性，不从图号赋统计量 |
| 摘要通路数字与对象错误 | pathway表2/3行为 T_5 IL2=323.19、IL6=243.88；根summary30行、README26行、H1正文21/22行仍15.88/17.70，且H1错归IL6为巨噬细胞协同 | 按亚群、方法、score、padj绑定来源；两项Treg结果不自动成为跨群协同 |
| CCI表修复，故事未同步 | 当前 ANXA1_FPR3 / TAM_0|TAM_0 /31.281；H1正文26行仍 TAM_0→T_5。H2正文9/10行仍宣称TGFB/CSF1固定互作分，当前14条没有CAF输入 | 更正对象；无对应结果的LR仅作候选问题，不能放在已计算支持证据中 |
| 下钻仍恢复旧改善指标 | 02/summary.json 5/6行0.234/20.7x；02/FINDINGS.md含62.1%、23.7x；01/summary仍24群/120DEG/T6 | 更新当前入口或显式退役旧视图，不能只改顶层卡片 |
| 五层表仍传播未计算结论 | systemic_impact_5layers.csv的CCI层称Spatial-constrained，区域层称因果生态位/溶骨，原样复制入payload | CPDB未使用空间距离；分别说明表达侧候选和描述性空间关联。五层须有真实覆盖、定位与缺口，不把文字总览称全层计算比较 |
| 机制与后期验证仍混入事实 | report399行“谱系追踪证实促转移”、419行“IFA验证逆Warburg”；H4正文25/30行同类扩写；payload带后期验证字段 | 当前故事按生信观察与候选机制组织；后期材料只外置关联，不能把定位或来源验证升级为功能因果 |
| 缺失底表仍返回成功包 | payload脚本42–70行仅exists后跳过；隔离空analysis测试仍退出0，无missing/status | 必需材料缺失要拒绝或明确partial/unavailable与原因；可选缺失不得被解释为已分析无结果 |
| 身份仍固定、来源不易检索 | payload脚本26–35行直接CRC/Xenium/50626/25；CSV有source hint，未形成当前提取/生成登记 | 从已登记样本或本轮小型源表清单读取，给可解析路径/键和hash；不要求每个CSV重复整份receipt |

“ANXA1同群候选”允许提出自分泌假说，但群级表达不证明同一个细胞发出并响应。score=100也不是置信度概率。T细胞周期程序支持增殖，不能直接证明克隆扩增或持续招募闭环。CAF4的当前ORA主要为COAGULATION/KRAS_UP/IFNα，不能把其他细胞的hypoxia图当CAF4已测得程序。CAF8统一为ECM-breakdown/基质重塑。

Treg90/64=1.40625、Tprolif2287/123≈18.5935可以按具体历史比较报告；当前表已注明343 vs13,088基因空间。顶层422为panel规模，不应与实际343分析轴混同。数量比支持DEG数量变化，不能直接升级为已校准灵敏度或污染清除；uncertainty仍为后续方法方向。

## 3. 旧内容清单与再生成风险

仅列举，未删除或修改业务内容。

| 文件 / 类别 | 仍存在什么 | 当前处理要求 |
| --- | --- | --- |
| analysis/README、三个槽位README/summary/FINDINGS | 旧24/120、假Z、旧Moran/降噪/灵敏度及旧H1/H4 | 当前交付入口更新；若冻结保留则明确历史且退出默认导航 |
| hypotheses/H1–H4 | 旧统计、错误方向、过度机制 | 同步正文，不只更新标题与Pill |
| moran_gain_by_lineage / contamination_clearance / noise_reduction_diffusion / diversity_entropy_delta 等 | 先前构造结果仍置于当前结果目录 | 未获真实来源的条目标不可用并退出当前事实层；不要求为清理重算 |
| organize_three_tier_analysis.py | 85–131行固定改善量，164–167行公式SNR/固定compactness，190–208行指定假说统计；默认写当前P2 | 禁止旧生成器写当前交付；整理器不得制造科学量 |
| run_unbiased_all_subcluster_analysis.py | 364–386行指定/随机Z、写旧槽位 | 明确legacy/demo边界，防止成为正式事实入口 |
| build_sample_analysis_delivery.py | 旧24群/改善文案生成 | 更新为引用登记结果或退出当前入口 |
| AGENT_GENERATION_SPEC.md示例 | 旧字段/常数，与本轮payload不同 | 更新或标历史示例 |
| 第二轮审阅文档 | 描述整改前payload/head(6)，不是当前状态 | 保留冻结观测，以本第三轮报告更新判断；已修复项不重复当未完成 |

REVISE生产侧过滤时序、effective-request/hash与mapping问题未被本轮commit改动，继续作为已知新生产运行限制；不重复设置为历史case/合作方开发的前置门槛。

## 4. 给执行 Agent 的有组织反馈

请保留本轮真实进展，将任务声明改为“载荷提取与展示资产同步通过，内容与部署收口partial”。下一步限定为以下六项：

1. **统一当前事实源。** 同步summary/README/HTML、三个槽位和H1–H4正文；修正通路数值、亚群归属与LR对象，所有关键数字可回到源行/键。当前证据不支持的TGFB/CSF1数字撤出，必要时作为按需分析问题。
2. **撤出当前虚构量。** 处理活动导航中的假Z与旧改善表/文字；隔离会重新写入这些内容的脚本。历史审阅与冻结材料保留原观测，不全盘删文件。
3. **让五层真实可读。** 声明实际focal覆盖与缺口，并给Gene/Pathway/CCI/区域的结果定位。区域叙述与程序连接，但不要称CPDB已空间约束或生态位因果已确立。无需全25群机械产出所有显著结果。
4. **补足消费最小信息。** 登记真实源路径/键、输入或源表hash、历史方法身份与本轮提取/生成命令、解释器、时间、产物hash；metadata取自登记，缺失显式可见。历史结果提取不冒充本轮科学重跑。
5. **完成真实运行部署说明。** 若要让合作方在启智再生成包，保护remote dirty与提交关系后部署对应代码，验证脚本实际可用；否则明确远端仅展示快照，并提供可用的本地生成位置。
6. **收口验收。** 做一次隔离再生成和逐层关键声明核对；提供当前33文件或修订后清单，而非复用旧32-file快照。临床/IFA/谱系保持体系外关联，机制仍待检验的部分按假说表述。

可继续给合作方：正式分析API、可追溯的历史T/TAM/CAF case，以及已核查底表的有界快照。当前整包不得作为已验收事实库；合作方不要自动摄入旧summary/FINDINGS/H1–H4。上述六项通过后即可接入证据包，阶段交付无需等待后期验证或新研究模块。

## 5. 本轮检查身份与边界

- task_id：`TASK-20261010-P2-THIRD-PASS-READONLY-REVIEW`；owner：审阅主Agent及五个独立子代理；upstream：REVISE历史case与执行方展示包；downstream：执行整改与合作方Agent开发。
- 本地REVISE：`revise-2.0 @ b4ea718bbc9c7171f654961c4d0fa51e2193d3da`，保留未跟踪迁移脚本；AA：`sl @ 6b158ed67117a580a1818cab957480537baeaa25`，审阅写文档前clean。本commit业务变更仅payload脚本，另12项为审阅资料。
- 启智cpu-2：`58bfcc6c-e848-44a7-a227-2ee82d25da01`，当次list为RUNNING；观测`2026-10-10T05:54:31.745292+00:00`（上海13:54:31）。REVISE checkout95c4cf2、AA checkoute33d122，tracking refs陈旧，不据ahead数字推断快进关系。无匹配重构/分析进程只是该次观测。
- 隔离验证解释器：`/Users/stephen/miniconda3/bin/python` 3.12.4；执行payload脚本指向临时样本、空analysis测试、py_compile。有效再生exit0且SHA相同；空样本也exit0是缺失处理问题。测试发生本轮，未单独记录秒级时间；临时副本已清理，不伪造精确时点。
- 输入：当前33个展示文件及指向的历史源表；H5AD未重新计算或全量重哈希。命令范围：只读git/status/show/ls-remote、qz list/run小型现场盘点、本地CSV/JSON读取、SHA比较与隔离payload生成。
- 输出：本报告和[本轮证据](third-pass-evidence/review-check.json)。未改业务代码、输入、现有结果，未提交/推送本次审阅，未向其他聊天或合作方发送。
- acceptance：载荷提取/同步局部`technical_verified`；整体`partial / review_complete_with_limits`，不授予生物学接受。
