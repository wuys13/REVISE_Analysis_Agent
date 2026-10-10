# REVISE Analysis Agent 文档总览导航

本目录管理下游分析 Agent 的系统架构、输入输出规范、交付物标准及历史审计演进记录。为避免规范与历史备忘混淆，文档划分为以下三个清晰层级：

---

## 1. 📘 当前生效规范 (Active Specifications - 开发者必读)

日常开发、添加新分析或维护交付管道时，以此层级文档为准：

| 文档 / 目录 | 核心主题 | 说明 |
|---|---|---|
| [`architecture.md`](architecture.md) | **系统架构与模块边界** | 模块分层、Raw/SVC 默认不配对原则、Agent 行为准则 |
| [`input-output.md`](input-output.md) | **输入输出与数值契约** | sample.yaml/project.yaml 格式、非负线性数值契约、CLI 规范 |
| [`delivery_framework/`](delivery_framework/README.md) | **标准交付框架规范** | 三阅读槽位（生信事实、重建影响、假说发现）与 L0-L4 倒金字塔交付标准 |
| [`analyses/`](analyses/README.md) | **内置分析能力库** | 空间重建影响、AUCell 通路活性、Moran's I 空间自相关详解 |

---

## 2. 🤝 跨仓生产消费协议 (Cross-Repo Contracts)

定义与上游生产仓 `REVISE` 之间的代码映射与验收边界：

| 文档 | 核心主题 | 说明 |
|---|---|---|
| [`cross-repo-review/contract-and-code.md`](cross-repo-review/contract-and-code.md) | **跨仓符号与代码对齐** | 上游生产与下游消费的代码行号、符号与配置对应关系 |
| [`cross-repo-review/evidence-and-gaps.md`](cross-repo-review/evidence-and-gaps.md) | **验收基准与边界备忘** | 冻结的双仓验收证据基线及当前未决边界 |

---

## 3. 📦 历史审计与演进留档 (Archived Audits - 仅供溯源)

阶段性审阅过程备忘与历史迁移记录，供溯源核查，不作为当前日常开发规范：

| 目录 / 文档 | 属性 | 说明 |
|---|---|---|
| [`cross-repo-review/20261010-architecture-review/`](cross-repo-review/20261010-architecture-review/README.md) | *[Archived Review]* | 2026-10-10 架构审阅与脚本整改过程记录（注：其中提及的硬编码旧脚本已在 commit `2759f93` 中清退，本目录仅供审计凭证保留） |
| [`migration.md`](migration.md) | *[Archived Log]* | 从单体仓向独立 Analysis Agent 拆分时的技术检查与 Playwright 验收记录 |
