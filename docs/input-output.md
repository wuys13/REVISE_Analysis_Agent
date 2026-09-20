# 输入输出与调用

每个样本为 `sample.yaml + raw.h5ad + SVC.h5ad`。Raw 与 SVC 的行数、顺序、ID 和基因集合可以不同。

```yaml
schema_version: 1
sample_id: example
files:
  raw: raw.h5ad
  svc: SVC.h5ad
expression:
  raw: {matrix: X, identity: unknown}
  svc: {matrix: X, identity: unknown}
columns:
  broad: Level1
  subtype: Level2
  reconstruction: SVC_cluster
spatial:
  key: spatial
  unit: pixel
  microns_per_coordinate: 0.27380817798463214
```

相对路径相对于声明它的 YAML。读取 cell-type 标签时统一将 `/` 归一为 `_`，保留缺失值且不改写输入文件；重建标签按原值读取。Raw `.X` 是上游交付的原始侧矩阵，SVC `.X` 是重建侧矩阵；“原始侧”不等于未经预处理的整数 counts。物理尺度仅在需要微米的分析中使用。

表达声明必须有来源依据。正式输入契约固定为 finite、nonnegative、unlogged linear `.X`；identity 已知且未声明旧 scale 时直接使用该契约。兼容旧 `untransformed` / `untransformed_nonnegative`，旧 `scale: unknown` 保持表达不可用；`log`、`log1p`、`log1p_nonnegative` 等旧 log 声明全部拒绝。identity 应描述真实矩阵来源。两侧声明互相独立，`.X` 缺失也可加载标签空间对象。确认的非负线性浮点值允许小数和小于 1 的值，不因数值很小而四舍五入、截断或改写；负值和非有限值由相应表达消费者拒绝。`identity: unknown` 允许标签与空间分析，但不能据此放行表达消费者。

公开接口：`load_sample(sample_yaml)`、`resolve_analysis_parameters(sample_yaml, analysis, *, project_yaml=None, overrides=None)`、`run_analysis(sample_yaml, analysis, output_dir, parameters)`、`run_batch(project_yaml)`。参数解析器按 sample < project < overrides 合并显式层；方法默认值仍由具体 analysis 所有。资源路径先相对声明它的 YAML 解析，再参与覆盖。若提供项目，样本必须真实列在该项目 `samples` 中。输出根下按 `<sample_id>/<analysis>/` 保存 `result.json`、表、图和报告，样本 `index.json` 负责导航。

`result.json` 保存实际参数、完成状态、未完成部分及原因、产物相对路径。Impact 的 `stages` 始终包含完整阶段表，尚未执行或因参数变化失效的阶段为 `pending`。失败不会把历史结果标为本次成功。可执行合成示例见 `configs/example_project.yaml`。

## Python 与命令行

```python
from revise_analysis import load_sample, run_analysis, run_batch
sample = load_sample('data/my_sample/sample.yaml')
result = run_analysis('data/my_sample/sample.yaml', 'spatial_autocorrelation',
                      'output', {'moran_n_neighbors': 6})
batch = run_batch('configs/project.yaml')
```

```bash
revise-analysis run --sample data/my_sample/sample.yaml \
  --analysis spatial_autocorrelation --output output
revise-analysis batch --config configs/project.yaml
revise-analysis report output/my_sample/spatial_autocorrelation
```

`run_analysis` 的 `output_dir` 是所有样本的输出根。直接调用中的 `geneset_path` 相对样本 YAML；批量配置中的资源路径相对项目 YAML。`geneset_path` 可配 `gene_set_names` 筛选 GMT；也可用 `gene_sets: {NAME: [GENE1, GENE2]}` 明确提供列表，两种方式不能混用。

项目配置使用显式样本列表与分析参数：

```yaml
schema_version: 1
output_dir: ../output
samples:
  - ../data/my_sample/sample.yaml
analyses:
  spatial_autocorrelation:
    moran_n_neighbors: 6
    sample_n_units: 30000
    random_state: 42
  pathway_activity:
    geneset_path: ../resources/h.all.v2025.1.Hs.symbols.gmt
    gene_set_names: [HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION]
    pathway_auc_threshold: 0.01
```

## 状态与重跑

- `succeeded`：请求部分完成。
- `partial`：有结果，但某些部分因为缺少标签、资源或方法支持而不可用。
- `skipped`：当前条件下没有可完成的分析结果。
- `failed`：执行错误；查看 `stage_errors` 或 `error.txt`。Impact 保留已完成阶段产物及 traceback，不把错误标成缺前提。Notebook 遇到意外错误直接抛出。

批量按样本和分析顺序执行，单项失败不阻断其他项。CLI 全部成功返回 0，部分完成／失败返回 1，非法调用返回 2。退出码不能替代对结果限制的阅读。

指定分析先在临时目录计算，再发布当前结果；历史目录保留在样本的 `.previous/`，以保护手写备注并避免旧表冒充新结果。该目录不进入结果导航，可在不再需要时由用户清理。不自动删除其他分析、输入或非本程序拥有的目录。不支持同时向同一样本目录写入。

每个分析自行检查科学前提。输入 loader 不猜测归一化、不建立全局匹配、不调用重建。负值和非有限表达由消费它们的正式计算拒绝。H5AD 的行／基因 ID 必须唯一。

Impact 的 `sample_n_units` 只约束表达 baseline/分子评分群体，不裁剪完整主 SVC 标签的 State。抽样分别作用于 Raw/SVC；Raw 用 `random_state`，SVC 用 `random_state + 1`，不要求相同抽样 ID。成员变化例外：以 Raw 选定范围中的共同 ID 读取完整 SVC 对应单位。scope 和基因集名称转换为产物文件名后必须仍唯一；碰撞会显式拒绝，避免覆盖。

## Impact 阅读分支与表达确认

确认 Raw 表达 identity 并满足固定线性 `.X` 契约后，Raw baseline、Gain、Raw K-control 和 Raw-defined membership 可以按各自前提运行；它们仍需要 SVC 的重建标签、空间坐标或共同 ID。Moran 与 AUCell 按侧独立消费表达：一侧满足条件即可保存该侧的原生结果，不能把另一侧缺失写成零值。跨侧比较（例如 shared-gene Moran 或 Raw/SVC 共同窗口的差值）只有在比较涉及的两侧条件都满足时才可用；单侧结果仍应保留并标明比较不可用。

## Impact 参数来源

样本 YAML 可在 `analysis_parameters.reconstruction_impact` 声明正式参数；项目 `analyses.reconstruction_impact` 覆盖样本同名值，显式 overrides 再覆盖项目。Notebook 调用统一解析器和 `effective_parameters`，分别展示 sample、project、override 与最终来源。多样本项目必须在 Notebook 显式选择 `SAMPLE_YAML`。样本资源路径相对于样本 YAML；项目资源路径相对于项目 YAML。`columns.reconstruction` 默认 `SVC_cluster`，不是 `Level2` 或 de novo Leiden 的别名。

Anatomy 的 `anatomy_window_side_microns` 与 parent 默认 `window_side_microns` 独立，初始均为 40 μm。按 parent 设置的尺度及可选控制见 Impact 方法页。推荐尺度只提供支持诊断，不覆盖显式值。

Anatomy 的窗口标签来源是 SVC 的 `columns.broad` 与 `spatial`；Raw origin 继续作为共享物理坐标原点，Raw/SVC 点映射保留 `Other`（已交付 SVC 格内未观察到 Tumor/Normal）与 `Unknown`（无对应 SVC Anatomy 格）的区别。

`expression.<side>.identity` 是上游来源声明字符串（如 measured_expression、reconstructed_expression）；`unknown` 表示尚未确认，不是根据数值自动检测的类别。scale 不用于声明另一种可接受的正式输入；已知 identity 对应的计算契约固定记录为 `untransformed_nonnegative`。
