# 输入输出与调用

每个样本为 `sample.yaml + raw.h5ad + SVC.h5ad`。Raw 与 SVC 的行数、顺序、ID 和基因集合可以不同。

```yaml
schema_version: 1
sample_id: example
files:
  raw: raw.h5ad
  svc: SVC.h5ad
expression:
  scale: untransformed_nonnegative
columns:
  broad: Level1
  subtype: Level2
label_aliases:
  Mono_Macro: Mono/Macro
spatial:
  key: spatial
  unit: pixel
  microns_per_coordinate: 0.27380817798463214
```

相对路径相对于声明它的 YAML。标签别名只用于分析视图，不改写输入文件。物理尺度仅在需要微米的分析中使用。表达声明必须有来源依据；不知道历史文件是否已归一化时不能把示例声明当作证明。

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
- `failed`：执行错误；查看 `error.txt` 与结果原因。

批量按样本和分析顺序执行，单项失败不阻断其他项。CLI 全部成功返回 0，部分完成／失败返回 1，非法调用返回 2。退出码不能替代对结果限制的阅读。

指定分析先在临时目录计算，再发布当前结果；历史目录保留在样本的 `.previous/`，以保护手写备注并避免旧表冒充新结果。该目录不进入结果导航，可在不再需要时由用户清理。不自动删除其他分析、输入或非本程序拥有的目录。不支持同时向同一样本目录写入。

每个分析自行检查科学前提。输入 loader 不猜测归一化、不建立全局匹配、不调用重建。负值和非有限表达由消费它们的正式计算拒绝。H5AD 的行／基因 ID 必须唯一。

抽样参数 `sample_n_units` 分别作用于 Raw/SVC；Raw 用 `random_state`，SVC 用 `random_state + 1`，不要求相同抽样 ID。成员变化例外：以 Raw 选定范围中的共同 ID 读取完整 SVC 对应单位。scope 和基因集名称转换为产物文件名后必须仍唯一；碰撞会显式拒绝，避免覆盖。
