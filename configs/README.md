# 配置用途与输入状态

相对路径由声明它的 YAML 解析。选择配置前确认实际 sample YAML 和两份 H5AD 存在；配置文件存在不代表输入已经交付。

| 配置 | 用途与使用条件 |
|---|---|
| `sample.template.yaml`、`project.template.yaml` | 新样本的正式接口模板；用实际交付身份填写，不从文件名推断表达来源 |
| `example_project.yaml` | 合成软件验证；先用 `scripts/create_example.py` 创建输入，不能用于科学结论 |
| `p2_project.yaml`、`p1_hd_project.yaml`、`p2_visium_project.yaml` | 2026-09-20 mini 验收配置；仅在对应历史交付存在的环境复用 |
| `*_full_project.yaml` | 全量交付模板；所列 `acceptance-full` 输入需由上游生成，不能当作已部署样本 |
| `p2_xenium_impact_20260922_*.yaml` | 固定历史运行/诊断参数，配合原 run 身份复核；不是当前远程结果的位置映射 |

正式 Notebook 默认使用 P2 mini 配置；其他样本显式设置 `PROJECT_YAML`，或单独设置 `SAMPLE_YAML`。未提供历史输入的远程 checkout 需要新项目配置，不能把另一次运行的矩阵替换到原配置路径。

历史 `data/P2CRC_Xenium/sample.yaml` 的绝对 provenance 路径保留当时含义；当前位置由[数据说明](../data/README.md)及上游证据解析。这个旧 carrier 的表达身份为 unknown，不放行表达计算。

当前综合阅读包在共同工作区的[发布目录](../../output/README.md)，属于派生视图。正式计算仍使用[输入输出协议](../docs/input-output.md)的 sample/index/result 接口。
