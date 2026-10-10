# 连续重建影响 Notebook

> [!IMPORTANT]
> **Notebook 双层公约：以哪个为准？**
> 1. **正式连续 Notebook 源码真理 (Source of Truth)**：本目录下正式连续分析 Notebook 的 `.ipynb` 是**唯一可编辑的 Notebook 源码模板**（Outputs 均为空，受 Git 管理）。该 Notebook 的阅读逻辑、展示参数和 Notebook 内调用编排只在这里维护；通用计算逻辑仍归 `methods/`、`analyses/`、`runner.py` 和 `batch.py` 等正式代码入口，外部独立探索不纳入本公约。
> 2. **运行时只读证据 (Execution Artifacts)**：`output/<run>/01_reconstruction_impact.executed.ipynb` 是自动化批处理或验收时由 `nbclient` 生成的**带图运行凭证**；`<run>` 由执行环境的 `OUTPUT_DIR` 决定。它受 `.gitignore` 忽略，切勿直接作为源码编辑。

本目录只保留当前正式分析的连续阅读入口；正式连续 Notebook 的可编辑源码只从这里维护。计算责任仍在
`revise_analysis/analyses/`、`runner.py` 与 `batch.py`；Notebook 调用同一套参数解析、阶段状态和结果保存协议，方便检查科学主线、参数来源和中间结果。

## 当前入口

| Notebook | 用途 | 结果入口 |
| --- | --- | --- |
| [`3_standard_pipeline/01_reconstruction_impact.ipynb`](3_standard_pipeline/01_reconstruction_impact.ipynb) | 读取一个已声明的 Raw/SVC 样本，连续展示重建影响、空间位置、分子程序和成员关系 | `output/<sample_id>/reconstruction_impact/result.json` 及其登记的 tables、figures、report |

正式流程之外的科学问题可以在外部工作区建立独立 Notebook。它们必须注明输入身份、表达前提、单位和基因范围、参数、输出位置及局限；独立探索不进入本目录的源码真理，也不会自动成为正式分析或结果入口。上游重建方法比较继续归属并列的 `REVISE` 仓库，本仓库不复制其 Notebook。

## 执行与参数

源码 Notebook 有意保持无输出。需要生成执行副本时，在仓库根目录使用项目环境中的 `nbclient` 和 `nbformat`：

```bash
SAMPLE_YAML=/path/to/sample.yaml \
PROJECT_YAML=/path/to/project.yaml \
OUTPUT_DIR=/path/to/notebook-output \
.venv/bin/python - <<'PY'
import os
from pathlib import Path

import nbformat
from nbclient import NotebookClient

source = Path("notebooks/3_standard_pipeline/01_reconstruction_impact.ipynb")
executed = Path(os.environ["OUTPUT_DIR"]) / "01_reconstruction_impact.executed.ipynb"
notebook = nbformat.read(source, as_version=4)
NotebookClient(
    notebook,
    timeout=None,
    kernel_name="python3",
    resources={"metadata": {"path": str(Path.cwd())}},
).execute()
executed.parent.mkdir(parents=True, exist_ok=True)
nbformat.write(notebook, executed)
PY
```

上述命令从仓库根目录运行；`SAMPLE_YAML`、可选的 `PROJECT_YAML` 和 `OUTPUT_DIR` 由执行环境决定。参数优先级和失效规则见[输入输出协议](../docs/input-output.md)。Notebook 执行副本只写入 `OUTPUT_DIR`；默认生产结果通过 `revise-analysis batch` 或 `revise-analysis run` 发布到 `output/<sample_id>/<analysis>/`，不由 Notebook 另建第二套结果目录。

## 阅读边界

- Raw 与 SVC 默认独立；只有明确的共同单位、共同基因或成员问题才做对应比较。
- 缺少表达、标签、坐标或可选 provider 时保留 `partial`、`unavailable` 及原因，不补零、不用随机 fixture 替代。
- `State`、`Gain`、Moran、program 和 membership 是不同证据视图；它们不能单独证明生物学改善。
- Web 和报告只读取本次结果 manifest 登记的文件；目录里残留的旧文件不自动成为当前结果。
