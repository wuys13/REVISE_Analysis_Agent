# 连续重建影响 Notebook

`01_reconstruction_impact.ipynb` 是针对单个 Raw/SVC 样本的可读科学工作台。阶段执行顺序仍由 `STAGE_ORDER` 控制，参数解析仍复用 `effective_parameters`；Notebook 的新增内容只整理阅读路径，不改变计算对象或科学定义。

阅读按三个问题推进：

1. **重建后呈现什么，两侧差异在哪里？** 输入小节读取 SVC 提供的重建状态标签，baseline 小节读取 Raw 独立表达分群与已提供的 Raw Level2。两侧各自的 Moran 与 AUCell 结果保留原生分支。
2. **State 与差异出现在哪里？** State 使用 SVC 主标签的局部多样性，Gain 只在共同有效窗口上对照 Raw Leiden；Anatomy 由完整 Raw 的 broad 标签独立建立，再按实际观测点解释 parent 窗口。
3. **位置与分子/成员有什么关系？** Moran 与 State/Gain 平行；AUCell 先固定单位级评分，再做窗口、Anatomy 或 Region 聚合；membership 只在 shared IDs 上比较；集成表连接已保存事实，不是独立验证。

Notebook 会显示每个阶段的中间产物、不可用原因和最终有效参数。`result.outputs` 是当前输出 manifest，阅读辅助函数只读取 manifest 中声明且实际存在的文件，不把目录中的旧文件当作本次结果。源码 Notebook 有意保持未执行状态。

默认样本是 `data/P2CRC_Xenium/sample.yaml`，默认 Notebook 输出是 `output/notebook/P2CRC_Xenium`。可通过 `SAMPLE_YAML` 和 `OUTPUT_DIR` 覆盖；`RECONSTRUCTION_IMPACT_OVERRIDES` 接受 JSON 参数映射。流程不会猜测表达尺度、替换输入标签或修改 H5AD。

如果需要生成执行副本，在仓库根目录使用项目环境中的 `nbclient` 和 `nbformat`，不要求 `nbconvert`：

```bash
SAMPLE_YAML=/path/to/sample.yaml \
OUTPUT_DIR=/path/to/notebook-output \
.venv/bin/python - <<'PY'
from pathlib import Path

import nbformat
from nbclient import NotebookClient

source = Path("notebooks/01_reconstruction_impact.ipynb")
executed = Path("output/notebook/01_reconstruction_impact.executed.ipynb")
notebook = nbformat.read(source, as_version=4)
NotebookClient(
    notebook,
    timeout=None,
    kernel_name="python3",
    resources={"metadata": {"path": str(source.parent.parent)}},
).execute()
executed.parent.mkdir(parents=True, exist_ok=True)
nbformat.write(notebook, executed)
PY
```

小型合成 fixture 适合检查调用顺序。真实样本必须根据上游证据声明表达尺度；非负线性浮点值（包括很小的值）可以保留，流程不因数值大小四舍五入。缺少 Level2、某一侧表达或 AUCell provider 时，只影响依赖它们的阶段，并在结果中保留不可用原因。正式结果与静态报告由批量运行器负责；报告只读取已保存结果。
