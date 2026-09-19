# 输入输出与调用

每个样本为 `sample.yaml + raw.h5ad + SVC.h5ad`。Raw 与 SVC 的行数、顺序、ID 和基因集合可以不同。

```yaml
schema_version: 1
sample_id: example
files:
  raw: raw.h5ad
  svc: SVC.h5ad
expression:
  raw: {matrix: X, identity: unknown, scale: unknown}
  svc: {matrix: X, identity: unknown, scale: unknown}
columns:
  broad: Level1
  subtype: Level2
  reconstruction: SVC_cluster
label_aliases:
  Mono_Macro: Mono/Macro
spatial:
  key: spatial
  unit: pixel
  microns_per_coordinate: 0.27380817798463214
```

相对路径相对于声明它的 YAML。标签别名只用于分析视图，不改写输入文件。物理尺度仅在需要微米的分析中使用。表达声明必须有来源依据；不知道历史文件是否已归一化时使用 `unknown`。支持 `untransformed_nonnegative`（归一化并 log1p）或 `log1p`（不再次变换）；identity 应描述真实矩阵来源。两侧声明互相独立，`.X` 缺失也可加载标签空间对象。历史共享 scale 可读取但不证明 identity。非负线性浮点值（包括小于 1 的值）是合法输入，不因数值很小而四舍五入、截断或改写；消费方按声明的 scale 处理，负值和非有限值才由相应表达消费者拒绝。

公开接口：`load_sample(sample_yaml)`、`run_analysis(sample_yaml, analysis, output_dir, parameters)`、`run_batch(project_yaml)`。输出根下按 `<sample_id>/<analysis>/` 保存 `result.json`、表、图和报告，样本 `index.json` 负责导航。

`result.json` 保存实际参数、完成状态、未完成部分及原因、产物相对路径。失败不会把历史结果标为本次成功。可执行合成示例见 `configs/example_project.yaml`。

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

确认 Raw 表达的 identity/scale 后，Raw baseline、Gain、Raw K-control 和 Raw-defined membership 可以按各自前提运行；它们仍需要 SVC 的重建标签、空间坐标或共同 ID。Moran 与 AUCell 按侧独立消费表达：一侧满足条件即可保存该侧的原生结果，不能把另一侧缺失写成零值。跨侧比较（例如 shared-gene Moran 或 Raw/SVC 共同窗口的差值）只有在比较涉及的两侧条件都满足时才可用；单侧结果仍应保留并标明比较不可用。

## Impact 参数来源

样本 YAML 可在 `analysis_parameters.reconstruction_impact` 声明正式参数；显式调用参数/项目参数覆盖同名值。Notebook 使用相同 `effective_parameters` 与阶段对象，并展示 overrides。样本资源路径相对于样本 YAML；项目资源路径先按项目 YAML 解析。`columns.reconstruction` 默认 `SVC_cluster`，不是 `Level2` 或 de novo Leiden 的别名。

Anatomy 的 `anatomy_window_side_microns` 与 parent 默认 `window_side_microns` 独立，初始均为 40 μm。按 parent 设置的尺度及可选控制见 Impact 方法页。推荐尺度只提供支持诊断，不覆盖显式值。

`expression.<side>.identity` 是上游来源声明字符串（如 measured_counts、reconstructed_expression）；`unknown` 表示尚未确认，不是根据数值自动检测的类别。`log1p_nonnegative` 是 `log1p` 的明确同义值，方法记录采用前者。
