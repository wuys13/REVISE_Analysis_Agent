# 迁移与验收

## 来源与边界

通用分析来源：REVISE `main@c83dc97d25b6d513b59cc301255e5bdc7e9c7cd9`。实施时复核 `revise-2.0@a26d36fe3b11d652f2cd9eeafb92e990518eb3ad` 的通用分析内容相同。

AUCell 增强与 Impact 来源：`reconstruction-impact@e82dd13ce013f3f120475e79892906c97cb114d3`。2026-09-18 再次核对分支指向。选择性复制，源仓库保持不动，保留 MIT 许可。

不迁入 `revise.svc` 服务、重建 backend、Raw Level2 mapping、旧 batch 或 Notebook builder。新包没有 REVISE runtime 依赖。

## 首次迁移记录（修订前）

| 类别 | 实现 | 验收依据 |
| --- | --- | --- |
| 架构与协议 | 独立 Python 包、样本与项目 YAML、公开 API、CLI、分层文档 | 输入原生轴和标签别名检查，输入 H5AD 不改写 |
| 计算工具 | 通用分析函数、Leiden、局部成员比较、窗口多样性、稳定阈值、原生 Moran、OmicVerse AUCell | 原实现回归与新非配对语义检查；可选提供者延迟导入 |
| 正式流程 | Impact、Moran、通路活性；顺序 sample × analysis 批量运行 | 单任务失败隔离；失败重跑发布当前失败记录；已有非生成目录受保护 |
| 结果审阅 | 结果索引、逐分析表图与静态报告 | 报告从已保存结果重绘；产物链接检查与代表图视检 |
| 科学工作台 | 连续 Impact Notebook，保留 `anatomy`、`partition`、`windows`、`state` 等中间变量 | 早期版本曾用独立 Jupyter 内核执行合成样本；当前源码的阅读展示后来已修订，需重新执行后再确认 |
| 独立环境 | 目标目录 `.venv`，Python 3.11 与 `[dev,pathway]` | 不安装 REVISE；Scanpy 分群及 OmicVerse AUCell 实际执行 |

首次迁移时目标仓库独立环境自动化测试 **39 项通过**；当时的旧版 Notebook（12 个代码单元）曾在目标路径实际执行，119 个本地报告链接与产物路径检查通过。随后 Notebook 的阶段展示和阅读文本发生了修改，因此这条历史执行证据不再证明当前源码与执行版 exact match；本轮不沿用它声称当前 Notebook 已执行，待根代理完成当前源码的执行与核对后再补充新的数字。其余测试包含 ID 不相交/重排、成员按 ID 对齐、无 Level2、常量基因、无通路交集/缺少 provider、共同窗口缺失不补零、State/Gain 成功和不可计算分支。通用迁入函数的签名与来源 AST 核对一致，计算函数体除可选 Squidpy 延迟导入外与来源一致；11 项通用回归通过。部分流程测试替换昂贵提供者或受控阈值；它们用于检查组合语义，不作为真实算法或生物学验证。

## 比较口径改变

这次不是仅改包名。Raw/SVC 分群、Moran、通路评分独立运行，允许单位数与 ID 集不同。Moran 每侧原生 kNN 建图，默认 6 邻居；AUCell 先固定单位级评分，再做空间聚合，不按 Region 重复评分；局部多样性共享物理窗口但独立单位与等量抽样。确认 Raw 表达后，Raw baseline、Gain、Raw K-control 和 Raw-defined membership 可按各自前提运行；Moran/AUCell 一侧表达满足条件即可保留该侧结果，跨侧比较需要两侧条件同时满足。成员变化独立开关，使用 Raw 定义范围内的共同 ID；SVC 的独立抽样不缩小这一局部对应范围。与旧强制配对结果不要求数值相同。

修订后 State 是输入 SVC 主重建标签在父群内的高多样性区域。Gain 只比较共同有效窗口的 SVC–Raw Leiden 多样性，正差异参与阈值选择。Raw Level2 是可选的独立局部 baseline，与主 SVC 标签的窗口差异只作描述，不混入 Gain。具体参数与解释边界见[方法索引](analyses/README.md)。旧 P2 数值或图形 parity 仍待找到同口径历史窗口并完成比较；技术结果刷新与用户对 State/Anatomy/组成的阅读确认是两个状态。

## 合成验收与真实样本的区别

合成样本由 `scripts/create_example.py` 生成：Raw 240、SVC 210 个单位，600 个基因，部分共同 ID、两侧无 Level2。正式示例配置各侧独立抽样 180 个单位，显式 DEMO_PROGRAM；Moran 与通路流程成功，Impact 因缺少 Level2、窗口不足以稳定确定区域阈值等原因保持 `partial`。实际观察和参数记录在 `output/example/`。Notebook 执行版位于 `output/notebook/example/`，源码 Notebook 保持无输出。

这些结果验证调用、计算和产物组织，不能支持真实生物学结论，也不能称为 P1 或全量验收。浏览器自动化未成功启动；报告完成静态链接检查与图像视检，未宣称完整浏览器交互验收。

## P1 表达输入仍待确认

P1CRC VisiumHD 历史 Raw 与 SVC 合计约 8.8 GB：Raw 507,684 × 18,085，SVC 424,433 × 12,926。SVC 全部 ID 可对应 Raw，顺序不同；按 ID 对齐坐标一致；两侧无 Level2。

现有元信息和可追溯生成记录**不能证明这些历史文件的 `.X` 是未归一化、未取对数的线性表达**。当前生产代码的行为不能反向证明历史产物。未将它们复制为已核验输入，未对它们执行分析；不猜测、不静默逆变换。

P1 后续需要上游确认表达与主重建标签，或交付合规 Raw/SVC，然后：

1. 复制至 `data/P1CRC_VisiumHD/`（H5AD 已忽略），从样本模板明确填写标签、坐标和尺度。
2. 使用项目模板的明确抽样规模运行 Impact（Moran/program 已集成，standalone 可按需单独运行），检查中间图与局部不可计算项。
3. 将实际结果与限制补到本页；抽样运行保持标注，只有实际运行过才声称全量验证。

P2 Xenium 标签空间链已按普通单对象接口接入；其表达验收仍需正式来源确认。见下方本轮记录。

## 本轮修订与 old P2 parity

P2 Raw 来源 `REVISE/raw_data/Real_application/P2CRC_Xenium.h5ad`；临时 SVC 来源 `REVISE/results/sc_SVC_case/P2CRC_Xenium/Fibroblast/spatial.h5ad`。复制后内容 SHA-256 与原件一致；来源和哈希保存在 `data/P2CRC_Xenium/sample.yaml`。临时 SVC 是普通单对象空间标签载体，不是新输入模式。未读取 `expr.h5ad`、未投影表达、未改写原件。

确认 Raw 340,837 × 422，Fibroblast 26,152；SVC 17,455 × 343，完整 10 类 `SVC_cluster`。两侧表达 identity/scale 均 unknown，Raw 无 Level2。真实可验收范围为完整主标签的支持、多样性、State、组成及与完整 Raw Anatomy 的关系。

| 关键能力 | 分类 | 新旧口径 / 本轮验收方式 |
| --- | --- | --- |
| 主重建标签 | 恢复旧定义 | 旧 P2 即使用 SVC_cluster；核对 17,455 单位与 10 标签，不以新 Leiden 替代 |
| 多样性与 bootstrap 阈值 | 保持算法 | min 4、200 draws、500 bootstrap；仅在输入范围/网格/随机过程均一致时比较数值 |
| Anatomy 与 parent 网格 | 保持旧语义 | 独立尺度；观测点落格后汇总，测试跨尺度、未知覆盖与并列 |
| Interface | 收紧权威定义 | 只由窗口 Tumor/Normal 共存生成；字面来源标签不直接生成 Interface |
| Raw/SVC 独立单位范围 | 必要适配 | 原生标签空间使用全量 SVC；共同 ID 只用于 membership，历史强配对结果不要求一致 |
| 尺度推荐 | 方法修正 | 平坦/退化曲线不报推荐；显式尺度有支持仍计算，不优化 State/Gain 结果 |
| Raw Leiden HVG | 方法修正 | 正确处理 log 输入，记录 seurat 与预处理；不要求复现旧错误 HVG 数值 |
| ΔNeff vs Raw Leiden | 对照解释修正 | 相同 parent 物理网格内共同有效窗口；完整正负差值，降低因果解释 |
| Raw K-control | 新增控制 | 固定 cohort/feature/PCA/kNN，有限 resolution，精确/最近 K 如实记录 |
| AUCell 空间聚合 | 明确评分群体 | 固定单位级评分重用；Region 数量不增加 provider 调用；缺失不补零 |
| Integrated Spatial Evidence | 新增关系表 | 从保存底表重算标签组成、区域与 Anatomy、changed units 和 program，明确分母 |
| 历史 P2 数值/图形 parity | 未比较 | 未找到可直接对应本轮范围和参数的历史窗口表；不填通过、不假定所有旧图一致 |
| P2 Gain / K-control / Moran / program | 待真实表达验收 | 方法行为可测试；真实输入表达条件未确认，当前不声明真实科学完成 |

软件行为、方法对照、真实标签空间链、真实表达分析和用户对区域定义/科学阅读的判断分别记录。合成成功、更多 State/Gain 或符合预期的现象不等于科学验收。正式单对象 SVC 到达后应更新来源与能力，重跑受影响结果。

## 前一轮计算修订的验收记录

以下记录保留此前的结果与方法证据。Notebook 随后增加了阅读转换和最终摘要展示，因此当时的执行版不能证明后续源码已执行；当前源码的重新执行结果见下方“科学阅读重组验收”。

| 层次 | 实际结果 | 证据 / 限制 |
| --- | --- | --- |
| 软件与回归 | 70 项测试通过 | 包括尺度独立、阈值失败、缺 Anatomy、固定评分、异常隔离、未知输入、报告只读、旧产物不冒充当前结果 |
| 方法执行 | 真实 Scanpy、native Moran、OmicVerse AUCell 在新版合成样本运行 | 当时曾验证非整数线性表达与旧 log1p 路径；本轮固定线性输入契约已取代旧 log1p 接口，旧结果不再证明当前输入声明。不是生物学验收 |
| Raw K-control | 固定 cohort/graph 的真实合成执行通过 | 精确与 nearest 的选择、seed 重现和主标签不变有测试；真实 P2 表达尚不能运行 |
| P2 主标签空间链 | 17,455 Fibroblast、10 类主标签，0 个阶段错误 | 原样输入来源见样本配置，不代表表达身份已确认 |
| P2 40 μm State | 4,653 个占据窗口；1,547 个有效；514 个 State 窗口 | Neff 阈值 2.4272768171；500/500 bootstrap 有效；95% CI 2.3346561680–2.5363448625 |
| P2 单位分母 | State 内 4,617；有效 State 外 7,830；支持不足未知 5,008 | 三者相加为 17,455；未知单位不混入 State 外 |
| 尺度推荐 | SVC 支持曲线建议 32 μm，计算仍用显式 40 μm | Raw/common 曲线当前只代表坐标潜在支持，不能称已计算的表达 baseline 支持 |
| Notebook/batch | 旧版执行版曾运行 10 个代码单元；此前 21 份科学 CSV/JSON 逐字节相同 | 后续展示修改曾使源码与执行版不一致；本行仅保留历史验收记录 |
| 保存结果的可追溯性 | 标签组成可由保存单位底表重算；每份报告 72 个本地链接通过 | 报告重绘前后科学文件哈希不变；关键图完成图像检查，未宣称完整浏览器交互验收 |
| 待确认项 | P2 Raw Leiden/Gain/K-control/Moran/AUCell/membership 的真实表达分析 | Raw 表达确认可推进 Raw baseline、Gain、K-control、membership；Moran/AUCell 按侧判断，共享基因比较需要两侧结果。未填通过 |
| 用户科学审阅 | 待阅读当前 State、Anatomy 与组成结果 | 是否最终采用绝对多样性或差值定位区域，仍不是自动化测试决定 |

重跑入口：

```bash
.venv/bin/revise-analysis batch --config configs/p2_project.yaml
```

- 正式结果：`output/P2CRC_Xenium/reconstruction_impact/report.html`。
- Notebook 源码：`notebooks/01_reconstruction_impact.ipynb`；默认使用同一 P2 样本参数。
- 实际执行版：`output/notebook/P2CRC_Xenium/executed_revised.ipynb`。
- 底表一致性与派生事实验收记录：`output/verification/p2_acceptance.json`。
- 新合成示例：`data/example_reconstruction/`；旧 `data/example/` 及其 H5AD 保留，不覆盖。

本轮表格中的数值用于说明当前定义及支持，不构成生物学改善结论。P2 总状态为 `partial`：标签空间链已运行，尚缺表达前提的组件没有冒充完成。

## 2026-09-19 科学阅读重组验收

网页按“重建后特征与两侧差异 → 状态与差异的空间位置 → 空间上的分子与成员关联”组织。执行阶段保持原有依赖；固定单位评分与已有区域在展示层汇合，关系表不作为新的独立验证。

| 分类 | 实际完成与证据 |
| --- | --- |
| 阅读层级 | 固定左侧问题/子问题/scope 目录；单 scope 不重复目录层级，多 scope 按配置展开；顶部最多五条事实链接到原证据 |
| 关联展示 | 保存表生成 Region × Anatomy、标签组成、program 分数地图/分组摘要、changed-unit 汇总及 shared-gene Moran 图；缺失不补零，各侧/范围分母分开 |
| 软件验证 | 全套 78 项测试通过；最终展示调整后 11 项相关测试再次通过。保留原有 4 条依赖弃用/稀疏矩阵效率警告 |
| 科学结果不变 | P2 正式结果 21 份、Notebook 21 份、合成结果 123 份科学表/JSON 与展示重组前逐字节一致 |
| 当前 Notebook | 使用 nbclient 实际执行最终源码的全部 10 个代码单元，无错误输出；源码/执行版代码单元完全一致，P2 Notebook/batch 的 21 份科学产物逐字节相同 |
| 浏览器验证 | 本机 Playwright/Chromium 实测 1440px 与 390px；目录、直接锚点展开、浏览历史、键盘折叠、图片加载通过，无整页横向溢出与页面脚本错误 |
| 产物完整性 | P2 11 张图、合成 81 张图各嵌入一次；P2 两份报告各 109 个本地/内部引用、合成 615 个引用存在，无重复 ID 或悬空锚点 |
| 仍未完成 | 真实 P2 表达分析、历史 P2 同口径数值 parity、用户对科学解释和阅读效果的认可；本轮没有将它们填为通过 |

证据入口：

- [源码/执行版与科学底表核验](../output/verification/reading_acceptance.json)
- [浏览器检查结果](../output/verification/browser_reading_qa.json)
- [P2 桌面截图](../output/verification/p2-desktop.png)、[窄屏截图](../output/verification/p2-mobile.png)
- [重组后的 P2 报告](../output/P2CRC_Xenium/reconstruction_impact/report.html)
- [当前执行版 Notebook](../output/notebook/P2CRC_Xenium/executed_revised.ipynb)

当前沙箱中的 Jupyter socket 与浏览器进程限制，通过本地执行权限机制解决；不依赖 nbconvert。浏览器技术验收不代表用户已认可阅读效果，也不代表生物学改善。

## 2026-09-19 上游生产核验与交互修订验收

本轮重新核对了上游当前实现，但没有修改上游仓库。`REVISE/reconstruct.py:40-68` 现有两条 Raw 路径：显式 `raw_adata` 会在 pipeline 前复制快照；文件路径会在 finalize 阶段重新读取原 source，并检查读取前后的来源文件身份一致。`delivery.prepare_raw` 只补充 obs aliases 和 uns，不写 `.X`；发布的 SVC 来自 pipeline 结果。因此，本仓库继续把 Raw `.X` 称为上游交付的原始侧矩阵，把 SVC `.X` 称为重建侧矩阵，并保持只核验、不改写输入。

sST 当前 runner `sc_svc_super_resolution_application.py:297-310` 将 parent-spot 校正后的 `SVC_X` 构造成 AnnData 保存，最终保存前没有再次 log；pipeline 内部仍包含 normalize，所以不能把它称为 raw counts。当前 `delivery.sample_document` 仍为该路径生成 `scale: unknown`：上游后续需按固定线性接口迁移已确认对象的配置、移除旧 scale 声明，并确认 Raw identity；分析端不增加新的 scale 选项。本轮证据只描述当前生产代码，不能反推历史 P2 文件已经确认，也不能据此放行 P2 Raw/SVC 表达分析。

本轮最终源码已完成以下验收。上面的旧 78 项测试与 10 个代码单元记录保留为历史；本表对应当前交互修订。

| 分类 | 本轮实际结果 |
| --- | --- |
| 输入与参数 | 固定线性接口；旧 log 声明拒绝、unknown 不升级；项目/sample/override 共用解析；应用参数和发布前检查 |
| 重跑正确性 | 当前记录替换，旧 extent/错误/图撤销；局部参数 patch；baseline 错误不阻断标签 State；阈值改变不重复 Moran/AUCell 或成员比较 |
| 科学修正 | baseline 先 scope 后抽样；Gain 保留 Raw/SVC 独立单位分母；新增 SVC 单位口径 Anatomy 条件表，未知与零分母保持缺失 |
| 逐节图形 | 输入、支持、多样性、Region、Anatomy、integration 在各自单元产生 1/4/2/3/2/3 张图；末尾不重复展示全部图片 |
| 软件验证 | 103 项测试通过，4 条依赖弃用/稀疏矩阵效率警告；最后绘图局部兼容调整的 14 项聚焦测试亦通过 |
| Notebook | nbclient 实际执行全部 17 个代码单元，无错误输出；源码与执行版代码相同；PROJECT_YAML 使用正式 P2 project |
| 表格核验 | P2 Notebook/batch 23 份表/JSON 逐字节一致；原 19 份科学产物与本轮前一致，2 份输入说明按接口更新，新增支持网格与 Anatomy 条件比例 |
| 页面 | P2 与合成报告在 1440px/390px 验证图片、目录、直接锚点、前进后退及键盘折叠；无页面错误或整页横向溢出。P2 16 张图、合成 100 张图均只嵌入一次 |
| 独立待验收项 | 真实 P2 表达分析、历史同口径数值 parity、用户对科学解释和阅读效果的认可 |

证据与交付：

- [本轮执行与逐节图形核验](../output/verification/interactive_acceptance.json)
- [测试日志](../output/verification/interactive_tests.log) · [Notebook 执行日志](../output/verification/interactive_notebook.log)
- [轻量验收附件](../output/verification/interactive_acceptance_bundle.zip)：代码版本/文件指纹、配置、环境摘要、核验结果、代表截图和关系图；不含 H5AD。
- [P2 报告](../output/P2CRC_Xenium/reconstruction_impact/report.html) · [执行 Notebook](../output/notebook/P2CRC_Xenium/executed_revised.ipynb)

本轮使用现有 nbclient 与本地浏览器完成执行；Jupyter socket 与浏览器启动按环境权限机制处理，不依赖 nbconvert。报告刷新只读保存结果。新增 Anatomy 条件比例是描述性关系，不能解释为富集或生物学改善。
