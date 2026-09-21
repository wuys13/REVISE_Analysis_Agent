# REVISE Analysis Agent

独立的重建后分析基础设施。输入 Raw 与统一 SVC，输出可继续计算、探索和审阅的结果。REVISE／上游负责重建和输入准备；这里负责分析。

## 双仓库工作区

本仓库与上游 `REVISE` 在共同父目录下并列放置。`configs/*project.yaml` 使用
相对路径消费兄弟仓库中的交付样本；远程工作区沿用同名的
`REVISE/` 与 `REVISE_Analysis_Agent/` 目录。公共工作区背景见
[根目录 README](../../README.md)，本地与远程协同规则见
[协同协议](../../docs/collaboration.md)。

## 从这里开始

| 要做什么 | 入口 |
| --- | --- |
| 让 GPT 协同审阅 REVISE 与本仓库 | [双仓库审阅入口](docs/cross-repo-review/README.md) |
| 理解仓库边界、核心设计 | [架构](docs/architecture.md) |
| 准备输入或调用 API／批量命令 | [输入输出协议](docs/input-output.md) |
| 找到某个科学问题对应的能力 | [分析能力索引](docs/analyses/README.md) |
| 看连续分析逻辑、调整参数和图 | [Impact Notebook](notebooks/01_reconstruction_impact.ipynb) |
| 核对迁移来源、真实验收限制 | [迁移与验证](docs/migration.md) |

## 安装与快速运行

使用 Python 3.10 或 3.11，在本仓库建立独立环境：

```bash
python3.11 -m venv .venv
.venv/bin/python -m pip install -e '.[dev,pathway]'
```

`pathway` 提供 OmicVerse AUCell；核心安装不要求通讯、轨迹或旧 Squidpy 工具依赖。可选能力见能力索引。

使用带 `SVC_cluster` 重建主标签的合成样本验证安装（已有目录和项目配置均不会覆盖）：

```bash
.venv/bin/python scripts/create_example.py
.venv/bin/revise-analysis batch --config configs/example_project.yaml
```

脚本默认新建 `data/example_reconstruction/`；旧的 `data/example/` 保持原样，不再作为修订后流程的输入。若该新目录已随交付提供，可直接运行第二条。需要并存多个样本时可传入安全的单段名称，例如 `--sample-id example_reconstruction_2`；已有 `configs/example_project.yaml` 会保留，需要自行把它指向新样本，或用 `revise-analysis run` 直接运行该样本。合成数据只用于软件行为验证，不能支持真实生物学结论；阈值支持不足等预期限制会在结果中逐项列明。

真实样本从 [sample.template.yaml](configs/sample.template.yaml) 和 [project.template.yaml](configs/project.template.yaml) 开始。P2 标签空间验收配置见 [p2_project.yaml](configs/p2_project.yaml)，分项证据见 [迁移与验证](docs/migration.md)。表达来源未知时仍能分析主标签与空间关系；正式表达接口固定为未 log 的非负线性值。

结果导航：`output/<sample>/index.json → <analysis>/result.json → tables/figures/report.html`。

Impact 报告按三个问题逐层阅读：**重建后特征与两侧差异 → 状态与差异的空间位置 → 空间上的分子与成员关联**。左侧目录定位问题与 scope，正文先给事实和关键图，支持诊断、参数和完整表格可继续展开。Moran 与单位级 program 评分是并行证据；program 的空间聚合与成员变化定位随后连接已有区域结果，不构成独立的生物学验证。

## 核心约定

- Raw `.X` 是上游交付的原始侧矩阵，SVC `.X` 是重建侧矩阵；本仓库核验并消费，不改写输入 H5AD。正式 `.X` 固定为非负、有限、未取 log 的线性契约；两侧 identity 分别声明，unknown 不阻止加载但不能放行表达计算。
- `SVC_cluster` 为重建主标签；State 用完整主标签，Raw Leiden 是表达条件满足时的独立 baseline。
- Anatomy 与 parent/State 网格共用物理坐标和原点，各自配置尺度；关系按观测点组成汇总。
- 默认分别分析 Raw 和 SVC，不强制观察单位配对；仅成员变化分析需要明确同单位关系。
- Notebook 展开科学主线、参数来源和中间结果；样本、项目与临时 override 显式合并，参数变化只失效受影响阶段。批量执行共享相同解析与计算模块。
- Web 只组织已保存的结果。算法改变、结构检查通过、真实科学验收分别说明。
- 不重建、不补 Raw Level2、不兼容旧 `spatial/expr` 载体、不实现 Agent 自动编排。

代码按 MIT 从 REVISE 选择性迁移，保留 [许可](LICENSE)。源仓库不受本次迁移修改。
