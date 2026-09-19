# Reconstruction Impact

主线是输入重建标签的空间证据。Raw 与 SVC 保持独立，只有成员对应使用明确共同 ID；坐标或行顺序从不用于推断单位身份。

## 三个阅读问题

### 1. 重建后呈现什么，两侧差异在哪里？

先读取输入提供的 `columns.reconstruction` 标签，确认每个 SVC 单位的重建状态及其空间载体。Raw 表达满足声明条件时，再在 Raw scope 独立建立表达 baseline；这一步描述 Raw 的表达分区，不能覆盖 SVC 主标签。两侧的 Moran 与 AUCell 结果保持各自的原生分支，单侧可用时不伪造另一侧结果。

### 2. State 与差异出现在哪里？

State 使用完整有效 SVC 主标签在 parent 窗口中的局部多样性。Gain 只在 Raw 与 SVC 都有有效窗口时比较 `ΔNeff`，并保留正负差值、支持数和阈值状态。Anatomy 由完整 Raw broad 标签独立建立，再按 Raw/SVC 的实际观测点组成解释 parent 窗口；它提供位置背景，不定义 State 或 Gain。

### 3. 位置与分子/成员有什么关系？

Moran 是与 State/Gain 平行的原生空间图分支，不依赖 Region，也不因 Region 再次建图。AUCell 先对固定单位集、基因轴和 gene set 做一次单位级评分，再把已保存分数按 window、Anatomy 或 Region 聚合。Membership 只在 Raw-defined cohort 与完整 SVC 标签的 shared IDs 上比较。Integrated Spatial Evidence 只连接已保存的 State/Gain、Anatomy、program 和 membership 底表；标签组成用于解释空间多样性，不构成独立验证。

### 表达确认决定哪些分支可读

确认 Raw 的表达 identity 和固定线性输入契约后，Raw baseline、Gain、Raw K-control 和 Raw-defined membership 可按各自的标签、空间与共同 ID 前提运行。Moran 与 AUCell 每侧独立：一侧表达条件满足即可保存该侧结果。跨侧比较（shared-gene Moran、Raw/SVC 共同窗口差值等）需要比较涉及的两侧条件都满足；否则保留单侧证据并明确比较不可用。确认的非负线性浮点值允许小数和小于 1 的值，流程不做隐式四舍五入；旧 `log`/`log1p` 声明被拒绝，identity unknown 仍不能放行表达分支。

## 标签、baseline 与解释

`columns.reconstruction` 默认 `SVC_cluster`，直接读取输入。称为“重建状态标签/重建分群”，不自动称 CAF biological subtype。State 使用 scope 内完整有效 SVC 主标签，表达抽样、HVG、Raw baseline 或可选 SVC Leiden 不改变它。Fibroblast 为 P2 首轮重点；其他 parent 使用同样能力。

表达条件满足时，在 Raw scope 独立计算 Leiden baseline。Raw Level2 只读取已提供列，缺失不映射、不推断。SVC Leiden 仅为可选 baseline，永不覆盖主标签。

主差值是 **ΔNeff vs Raw Leiden / Gain candidate**。两侧标签意义、群数、群大小与空间分布均可能不同，因此差值不能独立证明生物学改善。Raw Level2 的描述性差值另存。State 为主要阅读视图，最终选择绝对多样性或差值定位区域仍需真实结果审阅。

## 独立空间尺度

完整 Raw tissue 构建 Anatomy；每个 parent 有自己的 State/Gain 网格。共享物理坐标系和完整 Raw 原点，但尺度独立：Anatomy `anatomy_window_side_microns` 与 parent 默认 `window_side_microns` 初始均 40 μm。候选尺度为 16、24、32、40、56、80 μm。

Anatomy 唯一权威定义是窗口共存：Tumor 单独出现为 Tumor，配置的 Normal 来源（默认 Intestinal Epithelial）单独出现为 Normal，两者共存为 Interface，两者均无为 Other。上游单个单位的字面 `Interface` 仅保留来源，不直接生成窗口 Interface。

跨尺度关系按实际观测点：每个 Raw/SVC 点独立落入 Anatomy 格，再在该侧 parent window 汇总 Anatomy 数量和比例。保存主导类别、并列与未知覆盖；Other 与 Unknown 分开。State/Gain 阅读用 SVC 点组成，同时保存完整 Raw parent 组成及实际 Raw Leiden cohort 的独立组成（存在表达抽样时两者分母不同）。不能用中心点、面积权重或不同网格恰好相同的 window_id 替代此关系。

## 支持与尺度推荐

SVC parent 有效主标签单位为尺度推荐依据。候选曲线保存有效窗口数、保留单位比例；并列审计 Raw baseline 与 common-valid 支持及共同有效窗口占 SVC 有效窗口比例。Raw baseline 未计算时，只能称坐标潜在支持。

采用 log(scale) 与保留单位比例的支持曲线拐点，排除零有效窗口候选。平坦、退化、无明确拐点返回无推荐；等价候选确定性选较小值。不使用 State/Gain 面积、差值或 program 结果选择尺度，推荐不覆盖配置。无推荐并不阻断支持充分的显式尺度。

## 多样性与区域

每窗最少 4 单位、等量抽样 200 次，保存 Kobs、entropy、Neff、evenness。H = −Σp log(p)，Neff = exp(H)，evenness = Neff/Kobs；各侧使用自己的原生单位和确定性随机流。支持不足保留坐标与数量，指标保持缺失。

State 对 SVC Neff 选择稳定阈值。主差值在共同有效窗口保存全部正负 ΔNeff，仅正差值进入 Gain candidate 阈值。若配置 `sample_n_units`，Raw 与 SVC 先使用 `random_state` 和 `random_state + 1` 独立抽样；每个 scope 的 Gain 再以两侧实际入组的共同有效窗口为比较集合，同时记录 Raw、SVC 和共同有效支持分母，不能把任一侧抽样数冒充共同分母。沿用最少 200 个值、500 次 bootstrap、有效 bootstrap 至少 80%、95% 区间宽度不超过观测范围 25% 的规则。失败仍保存连续场、分布、支持、bootstrap 和原因；不把失败转成任意分位数或空的“成功区域”。

## Raw K-control 与成员对应

Raw K-control 为可选粒度敏感性分析，target K 来自该 scope 实际主 SVC 标签。固定 Raw cohort、预处理、HVG、PCA/kNN 和 seed，仅扫预配置有限 resolution 列表；精确 K 优先，多个精确候选优先接近主 Raw resolution。无精确候选可保存 nearest，明确实际 K、target K 与差值，不使 State 失败。它不替换主 Raw baseline，不修改主 SVC 标签，另存 ΔNeff 结果。它只控制群数，不控制标签意义或群大小等差异。

Membership 使用 Raw baseline cohort 与完整 SVC 主标签的共同 ID，保留 changed-unit map。即使 SVC broad label 改变，Raw-defined cohort 仍按 ID 比较。图表示分区对应差异，不等于生物身份改变；只用实际参与比较单位作 State/Interface 内外比例分母，不参与 State 定义。

## 分子与集成结果

Moran 在每个配置 scope 分别构建每侧原生图，保存完整基因表，共同基因对照另存 Raw/SVC/差值、单位数与基因分母。没有全局配对图。

AUCell 每侧固定评分单位集、基因轴、gene set 和参数，每个 program 只调用一次评分。保存单位 score、坐标、scope、覆盖率、provider 和实际 rank cutoff。window、Anatomy、Region 只聚合固定分数；保存有效评分数、总体分母、缺失，不补零生成差值。

Integrated Spatial Evidence 汇总 State/Gain × Anatomy、State 内外与 Anatomy 中主标签组成、program × window/Anatomy/Region、changed-unit × Region/Anatomy，并指回底表。changed-unit 的 State/Anatomy 上下文在 integration 阶段基于已完成 membership 与区域结果生成，membership 本身保持原生共同 ID 证据。Anatomy 图除主导类别外，还在 Raw/SVC 共享物理坐标系中显示 Tumor/Normal/Interface/Other/Unknown 的条件比例；比例使用各自 parent 窗口内实际观测点分母，不借用另一侧或另一 scope 分母。标签组成用于解释多样性，不作为独立验证。无需合并 Raw/SVC 表达矩阵。

## 共用执行与缺失

Notebook 和 batch 使用 `resolve_analysis_parameters`、`effective_parameters` 与 `ImpactWorkflow` 同一分阶段流程。显式参数按 sample < project < override 合并；方法默认值最后由 analysis 补齐。Notebook 显示每层来源和中间对象，不在样本路径错误时退回 example；多样本项目必须显式选择样本。

`ImpactWorkflow` 保留完整阶段状态表。`apply_parameters(...)` 的小字典是 patch 语义；Notebook 会传入重新解析后的完整最终参数快照，因此删除 override 也能恢复配置/default 值。它比较新旧最终参数，清除受影响阶段的登记产物和状态并返回失效阶段列表；`run_stage(name)` 拒绝消费 `pending/error` 前置。阈值 bootstrap 改变从 regions 起失效，不重算 Moran、AUCell 或原生 membership；窗口参数从 support 起失效，分子参数只失效 molecular 及其集成消费者。`render_figures(section)` 只读取当前登记表，为相应阅读小节登记图，不执行科学计算。Notebook 每个阶段完成后就近渲染和展示，末尾保存结果与综合事实，不自动补跑 pending 阶段。

每侧声明 expression identity；正式 `.X` 契约固定为 finite、nonnegative、unlogged linear。unknown 允许标签空间分析，但表达消费者保持 unavailable。已确认线性非负表达在方法内部复制后归一化和 log1p；HVG 使用相应方法并记录。预期缺前提单独记录；意外代码错误在 Notebook 抛出，在 batch 保存 stage error 和已完成产物。网页只读这些产物，不运行科学计算。
