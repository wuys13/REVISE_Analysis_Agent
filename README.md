# REVISE Analysis Agent

独立的重建后分析基础设施。输入 Raw 与统一 SVC，输出可继续计算、探索和审阅的结果。REVISE／上游负责重建和输入准备；这里负责分析。

## 从这里开始

| 要做什么 | 入口 |
| --- | --- |
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

使用合成样本验证安装（已存在时脚本不会覆盖）：

```bash
.venv/bin/python scripts/create_example.py
.venv/bin/revise-analysis batch --config configs/example_project.yaml
```

若交付目录已附带 `data/example/`，直接运行第二条。它是合成线性表达，不能支持任何真实生物学结论。Impact 可能因缺少 Level2 或没有足够窗口确定区域阈值而报告 `partial`，具体原因在结果中列明。

真实样本从 [sample.template.yaml](configs/sample.template.yaml) 和 [project.template.yaml](configs/project.template.yaml) 开始。**P1 历史文件的表达尺度尚未确认，不能把模板声明当成输入证明。**

结果导航：`output/<sample>/index.json → <analysis>/result.json → tables/figures/report.html`。

## 核心约定

- Raw/SVC 的 `.X` 是未归一化、未取对数的非负表达；重建值允许连续浮点数。
- 默认分别分析 Raw 和 SVC，不强制观察单位配对；仅成员变化分析需要明确同单位关系。
- Notebook 展开科学主线、参数和中间结果；批量执行共享相同计算模块。
- Web 只组织已保存的结果。算法改变、结构检查通过、真实科学验收分别说明。
- 不重建、不补 Raw Level2、不兼容旧 `spatial/expr` 载体、不实现 Agent 自动编排。

代码按 MIT 从 REVISE 选择性迁移，保留 [许可](LICENSE)。源仓库不受本次迁移修改。
