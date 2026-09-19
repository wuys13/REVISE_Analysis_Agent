# 连续重建影响 Notebook

`01_reconstruction_impact.ipynb` 是针对单个 Raw/SVC 样本的可读科学工作台。阶段执行顺序由 `STAGE_ORDER` 控制，参数通过 `resolve_analysis_parameters` 与 `effective_parameters` 解析；Notebook 只组织交互、阅读和显示，不改变计算对象或科学定义。

阅读按三个问题推进：

1. **重建后呈现什么，两侧差异在哪里？** 输入小节读取 SVC 提供的重建状态标签，baseline 小节读取 Raw 独立表达分群与已提供的 Raw Level2。两侧各自的 Moran 与 AUCell 结果保留原生分支。
2. **State 与差异出现在哪里？** State 使用 SVC 主标签的局部多样性，Gain 只在共同有效窗口上对照 Raw Leiden；Anatomy 由完整 Raw 的 broad 标签独立建立，再按实际观测点解释 parent 窗口。
3. **位置与分子/成员有什么关系？** Moran 与 State/Gain 平行；AUCell 先固定单位级评分，再做窗口、Anatomy 或 Region 聚合；membership 只在 shared IDs 上比较；集成表连接已保存事实，不是独立验证。

Notebook 会显示每个阶段的 `pending/completed/partial/unavailable/error` 状态、中间产物、不可用原因、关键参数和最终参数来源。`result.outputs` 是当前输出 manifest，阅读辅助函数只读取 manifest 中声明且实际存在的文件，不把目录中的旧文件当作本次结果。源码 Notebook 有意保持未执行状态。

默认样本是 `data/P2CRC_Xenium/sample.yaml`，默认 Notebook 输出是 `output/notebook/P2CRC_Xenium`。可通过 `SAMPLE_YAML`、可选 `PROJECT_YAML` 和 `OUTPUT_DIR` 覆盖；`RECONSTRUCTION_IMPACT_OVERRIDES` 接受 JSON 参数映射。项目只有一个样本时可从项目声明解析；项目含多个样本时必须显式设置 `SAMPLE_YAML`，不会静默选择第一个。

参数优先级为包默认值 < sample YAML < project YAML < Notebook override。修改 `OVERRIDES` 后先重新运行“应用参数并查看失效阶段”单元。该单元重新解析所有来源，把完整 `EFFECTIVE_PARAMETERS` 快照传给 `workflow.apply_parameters(...)`；这样删除 override 也会恢复 sample/project/default 值。已完成但依赖改动参数的阶段回到 `pending`，无关阶段保持当前。阶段运行前 Notebook 会再次解析并比较待应用参数，若与 `workflow.parameters` 不一致则要求先应用。`run_stage` 会拒绝消费失效前置，成功或科学不可用后按节调用 `render_figures(section)`，所以图与解释就近出现；末尾只保存清单、报告和综合事实，不补跑科学阶段。

例如只改变区域 bootstrap 阈值会失效 regions 及其下游 integration/figures，不要求重算 Moran、AUCell 或原生 membership。改变窗口/支持参数则从 support 起按真实依赖失效；改变表达抽样或分子参数只重算对应 baseline/molecular 分支及其集成消费者。

如果需要生成执行副本，在仓库根目录使用项目环境中的 `nbclient` 和 `nbformat`，不要求 `nbconvert`：

```bash
SAMPLE_YAML=/path/to/sample.yaml \
PROJECT_YAML=/path/to/project.yaml \
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

小型合成 fixture 适合检查调用顺序。Raw `.X` 表示上游交付的原始侧矩阵，SVC `.X` 表示重建侧矩阵；本仓库不改写它们。正式 `.X` 契约固定为 finite、nonnegative、unlogged linear，真实样本必须根据上游证据声明 identity。确认的线性浮点值（包括小数和小于 1 的值）可以保留；旧 `log`/`log1p` 声明会被拒绝，`identity: unknown` 也不会放行 P2 表达消费者。缺少 Level2、某一侧已确认表达或 AUCell provider 时，只影响依赖它们的阶段，并保留不可用原因。正式结果与静态报告由批量运行器负责；报告只读取已保存结果。

## 独立问题的 Notebook

`01_reconstruction_impact.ipynb` 是 Impact 的连续主线，不是所有下游实验的唯一容器。Agent 可针对额外测试建立独立 Notebook，优先复用包内方法，并明确问题、输入来源及表达前提、预处理、比较范围、随机参数、资源、输出位置和解释限制。直接方法调用的前提与结果记录由 Notebook 显式承担；不修改输入，不将缺失补为零。稳定后再下沉计算和流程。

上游已有[四方法 assembly 比较 Notebook](../../REVISE/reproduce/case/assembly_comparison.ipynb)，本轮只引用，不复制或迁移。它的共同 ID/基因比较是局部研究约定，不改变本库独立 Raw/SVC 规则。真实输入与协同审阅见[跨库入口](../docs/cross-repo-review/README.md)。
