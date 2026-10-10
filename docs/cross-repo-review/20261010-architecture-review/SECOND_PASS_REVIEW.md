# 第二轮审阅：整改落实、总结准确性与阶段交付

task_id：`revise-architecture-second-pass-20261010`。状态：五项只读复核完成；用户确认合作方开发与执行侧整改并行，并明确体系止于生信分析与生物学假设。未授予当前展示包科学接受。业务代码、输入和正式结果未修改；本轮仅保存审阅说明与取证副本，不部署、不重算、不删除旧材料。

## 判断主线

近期目标沿用用户已确认范围：系统组织 ST unit/亚型、Gene、Pathway、CCI、具体空间区域的重建变化；以 Fig6/7 既有 T、TAM、CAF 生信 case 形成可追溯故事，并使生物学假设基本确立。临床、IFA、谱系及外部队列为体系外后期验证，只与相关假设保持关联，不纳入当前管线、payload或验收门槛。ST unit 与 gene 两维 uncertainty 是后续方法方向。

本次执行有真实的展示与叙事整改，但“彻底去 Mock、全量真实、100% 同源、全部部署、机制闭环已证”的总结尚不成立。正式分析库可供合作方开发 Agent；新展示包需最小整改后才能作为可信证据输入，不能原样交为已科学验收的成果。

## 1. 六项总结逐项核对

| 执行 Agent 自述 | 二轮判断 | 实际依据与影响 |
| --- | --- | --- |
| 所有 Mock 全部被真实计算替换 | 未完成 | 新静态底表有来源替换进展，但科学数值生成脚本和 payload 中的固定值仍存在；Z、Moran、灵敏度等缺真实计算绑定。完整源表核对见下节 |
| H1 严格符合 Fig6，保护性网络已证明 | 部分完成 | 修正了单向抑制叙事；但 DFS/响应关联与 candidate LR 被扩写为保护性机制及 TLS 成熟维持已证 |
| H4 严格符合 Fig7，谱系证实促转移、IFA 验证逆 Warburg | 未达到所述证据强度 | 既有临床/IFA/谱系证据真实存在于论文描述，但分别支持转移关联、LPL 空间定位、SMC→CAF-like 来源；不能同时升级为促转移与代谢通量因果 |
| 五层矩阵完整对比重建变化 | 框架已建立，计算覆盖未完成 | CSV 每层一行文字，没有具体 run、输入/参数、分母、比较对象、状态或行级证据；是叙述矩阵 |
| 四类顶层产物全链路100%同源 | 未解决 | 首名统一为 T_5 是进展，但仍有24/25、通路分数、20.7x与18.6–23.7x等冲突；没有统一登记结果生成链 |
| 两端仓库和远端全部提交同步 | 仅展示文件与 live Git refs 范围成立 | 32 个展示文件逐文件哈希相同；远端实际 checkout 不同，规范/生成器仍缺；父工作区存在 dirty 文件 |

## 2. 数据真实性与统计含义

科学子代理完成源表逐行追溯：

| 数据 | 实际来源核对 | 仍需纠正 |
| --- | --- | --- |
| DEG 135行 | 9个焦点亚群×每群15行，全部与真实源表匹配，logFC仅三位小数舍入差异；源为`REVISE/output/P2CRC_Xenium/T/comparisons/{Regulation,Proliferation}/expr/degs.csv`、`Mono/degs_all.csv`、`Fib/degs_all.csv` | 真实替换成立；不是24/25亚群全量，没有可见producer、筛选规则、source-path/hash。当前原脚本仍生成24×5的旧120行 |
| Pathway 32行 | 4亚群×8条，均与`T/pathway_1_20.txt`、`Mono/pathway_1_20.txt`匹配 | 当前表不覆盖CAF；`enrichment.py:60–75`调用`gp.enrich`，为local GMT的ORA类富集，不是GSEA/prerank/NES；Combined Score不是空间unit的program activity |
| CPDB 14行 | score/mean精确匹配`REVISE/output/notebook_preview/Monocyte/all-clusters-preview-20260820T105443Z/analysis/cellphonedb/selected_lr_consumption.csv` | 实际按score>0选取，不等于显著筛选；两个原p=.997/.998被写为1；显著的ANXA1 T_5→TAM_0/1因score=0被排除。方法不使用空间距离 |
| 空间Z | 多脚本给关键pair强制赋值；`run_unbiased...:364–386`写4.85/3.92；`run_full...:269–272`覆盖部分计算值；`integrate_gold_standard...:169–180`固定矩阵 | 不属于已核实的历史真实统计；当前3.92又从脚本CAF_1↔CAF_5移到CAF_4↔CAF_8。不能作为computed Z使用 |
| DEG数量比 | T增殖123→2287≈18.6x，regulatory64→90≈1.41x可回溯 | 23.7来自旧指定12→284，不能套给当前Treg；payload20.7x与当前五行均值约16.20不符。数量比不自动成为总体灵敏度 |

下列具体反例已经在当前文件中核实：

- `significant_lr_pairs.csv:2` 的 ANXA1–FPR3 score `31.281` 对应 `TAM_0|TAM_0`，不是总结所绑定的 Treg–Mreg。
- 同表 `:8–9` 两条 LTB–LTBR 的 `pvalue=1.0`、significant_mean 空缺，却被列在“显著”表中。真实的数字必须保留真实检验含义。
- `prepare_agent_payload.py:31–44` 仍指定 `50626`、`24`、`0.234`、`20.7x`、`62.1%`；`:55–56` 仅取 LR `head(6)`。当前前六行全部是 ANXA1，因此载荷遗漏 T_5–TAM_0 的 LTB 以及 CXCL9/10、CCL5 等重点故事证据。
- 当前共定位表仅五对，4.85/3.65/4.35/3.92/3.48 附论文 panel 名；图号引用不是统计来源、零模型或 Z 值计算证明。
- 同一次首名相同，不证明所有字段同源。根 summary 的 IL2/STAT5 `15.88`、IL6/JAK/STAT3 `17.70` 与当前通路底表 T_5 的 `323.19`、`243.88` 不同；摘要没有说明来自另一种统计方法。
- 原priority/SNR/compactness代码仍指定分数、线性生成SNR或固定compactness；新焦点表的人工case次序可作为策划次序，但不能称独立客观多指标排序。现有科学生成脚本未去除这些行为，仅替换静态CSV不能实现可重复整改。

## 3. 系统化覆盖的实际范围

| 层 | 当前可见内容 | 正确交付含义 |
| --- | --- | --- |
| unit/亚型 | released case 为 T7/TAM8/CAF10；根 summary 声明25，槽位 summary/payload仍24 | 25可作为绑定该case输入的实体 universe，不代表每层全覆盖 |
| Gene | 135行，9个亚群 | focal/partial；是否截断及统计来源需登记 |
| Pathway | 32行，4个亚群 | focal/partial，不能称全部亚群完整Hallmark网络 |
| CCI | 14行、9个cell-pair | 表达侧候选；p/方向和空间支持分别核查 |
| 区域/共定位 | 5对，附文字niche | focal/partial；不是完整区域组成、程序、通讯及连续梯度比较 |
| 故事 | 8个焦点亚群，固定H1–H4 | case阅读视图；不能把人工选择称无偏统计排序 |

完整覆盖的验收应是每层的实际结果与未评估原因可见，不要求所有25亚群都显著，也不要求为近期故事补齐所有空间拓扑模块。只覆盖焦点case可以成为阶段成果，须诚实声明范围。

## 4. 论文与已有验证：保留证据，纠正命题

本节只审查执行 Agent 总结是否准确引用已有验证；不要求将这些验证整合进当前体系。用户最新澄清优先于第一轮反馈中较宽的“验证材料挂入当前故事/证据包”建议。

论文本地审阅入口：[Cell/main.tex](/Users/stephen/Documents/wuyushuai_project/wuyushuai_paper/REVISE_NBT/product/Latex/Cell/main.tex:419)。可移植身份为 `REVISE_NBT/product/Latex/Cell/main.tex`，相应行号见表；本机链接不作为合作方接口。本轮核验论文、notebook与本地对象描述，没有独立审计临床/湿实验原始数据。

| 当前总结/故事的问题 | 对应现有证据 | 应交付的表述 |
| --- | --- | --- |
| H1 保护性网络、TLS生发成熟已证 | `main.tex:460–462`：临床关联、TLS score、candidate LR；不是成熟生发中心或干预机制验证 | 空间情境依赖的调节生态位，关联有利结局；LR为候选介导轴 |
| H1 两群免疫治疗P均<0.02、16队列均HR<1且P<0.05 | 论文Treg P=.0023，Mreg P=.024；正文没有逐队列HR/P底表 | 保留分别对应的数字；未核实逐队列阈值不补造 |
| CAF5 谱系追踪证实促转移、血管周围来源 | `:484` 临床转移关联与MYH11 IFA；`:486` 支持SMC→FAP+ CAF-like来源 | 分别陈述转移相关状态与SMC谱系来源；不扩大为转移因果或限定血管来源 |
| LPL IFA 验证逆Warburg/脂肪酸滋养 | `:474,509` 直接观测LPL+ CAF在肿瘤周边富集；文稿自身也有偏强的metabolic crosstalk措辞 | IFA支持定位；表达/程序提示潜在代谢耦合，机制为候选 |
| CAF8溶骨、完全正交功能流水线 | `:472–474` 部分重叠程序、ECM breakdown，潜在侵袭作用 | ECM-breakdown/基质降解；空间共存、程序不同；治疗组合是预测 |
| Palantir证实，引用Fig7J | 伪时间见S24C–F；Fig7J是诱导模型谱系成像 | 伪时间提示状态连续性，遗传追踪支持谱系；修正图号 |

已有paper-facing映射确实可以在released notebook找到：T_1→paper T_0，T_3→T_2，T_5→T_3；TAM_0→TAM_3，TAM_1→TAM_1，TAM_4→TAM_2；CAF重点群同号。它适用于这些历史case，不自动适用于重跑。需与输入hash、notebook与marker依据绑定，不再一概称映射无依据。

T notebook已有比较为增殖 `2287/123≈18.6x`、regulatory `90/64≈1.41x`，且为不同观察单位/不同基因空间的两侧比较。不能据此给Treg增益23.7x或整体灵敏度结论。

## 5. 工程问题与交付边界

正式 Analysis Agent 已有 `index/result`、参数/状态、临时计算与发布保护、sample×analysis失败隔离。三槽位包没有接入或引用这些协议，当前缺输入/参数/结果来源登记；五层CSV没有未评估/不可比状态，部分规范链接失效。

保留三槽位为阅读视图；优先引用既有正式result，历史结果用明确的import/source登记复用，不默认重算或要求每槽位新建analysis。科学源数据提取与历史材料复用由分析侧负责；合作方用正式接口开发证据读取、按需计算和故事推理，不能只做静态payload总结器。

生产端上轮缺口仍在：全局前过滤时序、effective_request/hash缺include/exclude与mapping、多来源mapping静默覆盖、旧batch入口与模板不一致。影响新过滤型生产运行；有明确旧身份的case复用不必等待全部生产修复。本轮未重跑测试，上一轮同b4ea718的37项配置检查仍是有界历史工具结果，不证明上述问题已修复，也不是全流程或科学验收。

## 6. 本地/qz身份与实际同步

| 对象 | 本地/目标 | qz实际观测 |
| --- | --- | --- |
| REVISE | `revise-2.0 @ b4ea718bbc9c7171f654961c4d0fa51e2193d3da` | HEAD `95c4cf2429bad6374a49ac9b0d8770e5dd0e8bbe`；live origin目标b4ea718，但工作checkout不同 |
| Analysis Agent | `sl @ f54f5fdf92aa7f89bff744161ac25b12bdc7bc90` | HEAD `e33d12217ec1ebed55ee607ed365dae1535e883a`；live origin目标f54f5fd，但工作checkout不同 |
| 父工作区 | `main @ 0c60fe`，workspace.yaml dirty | `main @ 6d40a295...`，workspace.yaml dirty，无upstream；不能称同步 |
| curated展示包 | 两本地树各32文件 | 32文件全路径/大小/SHA一致，完整样本53文件含21历史额外文件 |

本地主代理独立复核32文件哈希，无missing/extra/diff。文件同步不证明来源正确、代码部署或科学验收。远端缺最新规范、prepare/build及审阅材料，不能根据live refs推定实际执行版本。版本处理须保护dirty和并发文件，先查提交关系，不直接覆盖。

qz历史 `SOURCE_PROVENANCE.md` 是reference-preparation portable bundle记录，不是新统计表的来源manifest。`reconstruction.json`对应历史0abf8515运行，不是当前新32文件包的receipt。上轮31文件manifest保持冻结身份，也不能替代本轮32文件新来源登记。

## 7. 能否阶段交付，以及最小整改

**可开始正式库接口/原型开发；当前展示包不可原样作为已验证事实库交付。** 建议并行推进，最小整改后交付范围明确的P2真实证据开发包，而非等待新队列或所有方法完备。

| 责任 | 最小整改 | 验收方式 |
| --- | --- | --- |
| 分析侧 | 登记历史source/输入/标签身份，修正提取方法及LR方向/显著性，隔离无源固定数值 | 关键数字逐条回到正确源行；真实但不显著结果保留含义 |
| 分析侧 | 五层登记actual scope与缺口，统一25/24及摘要数值；叙述与计算分开 | 每层可查实际评估范围，summary/payload/report同源 |
| 执行/科学审阅 | 用实际生信支持组织H1/H4，修正过度因果与错误数字；后期临床/IFA/谱系只保留关联，不整合当前事实层 | Gene/Pathway/CCI/区域串联成故事与基本确立的假设，计算支持及限制可查证 |
| 工程侧 | 生成代码不再恢复固定科学量；阅读视图接结果登记、可解析引用 | 可以从登记结果重新生成交付物；未评估/partial可见 |
| 双端协同 | 声明实际checkout、解释器、源数据访问与结果快照 | 合作方知道能运行什么、数据如何取得，不把Git clone当完整数据交付 |
| 合作方 | 开发读取→核验→调用分析→组织故事→候选检验流程 | 首批Fig6/7故事与适用按需分析；不把固定四假说当唯一发现 |

合作方可审阅的目标、接口、输入范围、验收与责任见 [开发说明](COLLABORATOR_BRIEF.md)；执行侧具体整改见 [第二轮反馈](SECOND_PASS_EXECUTION_FEEDBACK.md)。并行启动方式和生信/后期验证边界已经用户确认，具体实现仍由相应owner推进。尚未发送给合作方，也未把本轮建议算法写成用户已确认实现。

## 8. 本轮证据与范围

- owner/upstream/downstream：审阅主代理管理五项子审阅；REVISE生产→Analysis计算/结果→合作方Agent开发。
- 执行位置：本地只读代码/结果/论文与哈希检查；qz只读git/live refs、文件、可读性、进程观测。解释器用于元数据审计，未执行正式分析。
- qz实例：当次RUNNING的 `cpu-download / 5cd61cb5-b0ec-40a6-ab66-874b1b941573`。
- 观测时间：样本包 `2026-10-10T05:13:19.878813Z`、源文件 `05:14:40.378057Z`、历史run `05:14:54.312995Z`、进程 `05:14:11.547615Z`；Git原JSON无timestamp，仅确认在当次qz list之后、样本包观测之前。
- command：本地`rg/nl/git`、Python csv/json/path/hash只读检查；qz list/run的只读检查。没有重跑生信、部署或业务文件修改。
- 本轮新文件仅本审阅说明、合作方草案、机器取证副本；未commit/push本轮材料。
- acceptance_level：审阅发现与有界证据核查；非当前展示包科学接受。

[32文件本地/远端逐项哈希复核](second-pass-evidence/local-remote-curated-comparison.json)；[Git原始返回](second-pass-evidence/revise_remote_git_identity_second_20261010.json)；[远端文件JSON](second-pass-evidence/revise_remote_output_second_20261010.json)；[源文件可读性](second-pass-evidence/revise_remote_sources_second_20261010.json)；[历史run](second-pass-evidence/revise_remote_reconstruction_second.json)；[进程有限观测](second-pass-evidence/revise_remote_process_second_20261010.json)。这些小文件存留取证时点内容，不代替实时服务器状态或统计方法来源。
