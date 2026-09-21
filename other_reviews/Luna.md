目前已完成的是“独立分析仓库的第一版基础设施”，核心代码、文档和示例都已落地。

## 一、已经完成

### 1. 仓库与架构

已建立：

[REVISE_Analysis_Agent](../README.md)

包括：

- `revise_analysis/` 独立 Python 包
- `methods/` 可复用科学计算
- `analyses/` 正式分析流程
- `plotting/` 图表
- `reporting/` 静态报告
- `notebooks/` 连续科学 Notebook
- `configs/` 配置模板
- `data/`、`resources/`、`output/`
- 独立 `.venv`
- `README.md`、`AGENTS.md` 及分层说明文档

新包不依赖 REVISE runtime，也没有迁入 reconstruction backend、旧 batch、旧 `spatial/expr` 兼容层或自动 Agent 编排。

### 2. 输入输出协议

已实现：

- `sample.yaml + raw.h5ad + SVC.h5ad`
- Raw/SVC 独立读取
- 不自动配对、不自动归一化、不自动取对数
- 标签别名，例如 `Mono_Macro → Mono/Macro`
- 空间坐标和微米尺度配置
- `load_sample`
- `run_analysis`
- `run_batch`
- CLI 命令

结果结构已经固定为：

```text
output/<sample>/
├── index.json
└── <analysis>/
    ├── result.json
    ├── tables/
    ├── figures/
    └── report.html
```

失败任务会记录为当前失败，不会把旧结果误标为本次成功；旧结果会移入 `.previous/`。

### 3. 三个正式分析入口

已经完成：

- `reconstruction_impact`
- `spatial_autocorrelation`
- `pathway_activity`

Impact 主线已经包括：

- Raw/SVC 独立概览
- 独立 Leiden 分群
- 可选的 Raw 定义成员变化比较
- Raw Level1 anatomy
- 共同物理窗口
- 两侧独立局部多样性
- SVC State
- Raw→SVC Gain
- Moran
- 通路活性
- 中间表格和诊断图

成员变化只使用明确的共同 ID；不会把全流程变成强制配对。

### 4. 科学口径

已经写入代码和方法文档：

- Moran 两侧分别建立原生坐标 kNN 图，默认 6 邻居
- State 是 SVC 自身的高多样性区域
- Gain 是共同有效窗口中的 `SVC − Raw` 多样性差异
- 缺失窗口不补零
- Raw Level2 缺失只影响 Level2 baseline
- Raw Level2 baseline 与 SVC Leiden 的对照不混入 Gain
- AUCell 使用真实 OmicVerse provider，没有 fallback 评分
- 可选依赖缺失时记录 unavailable 或 partial

相关说明见：

[Impact 方法说明](../docs/analyses/reconstruction-impact.md)

[分析能力索引](../docs/analyses/README.md)

### 5. Notebook 和 Web 审阅

已完成连续 Notebook：

[01_reconstruction_impact.ipynb](../notebooks/01_reconstruction_impact.ipynb)

它保留了：

- `anatomy`
- `partition`
- `windows`
- `state`
- `gain`
- `moran`
- `AUCell`
- 综合 observations

Notebook 已在目标仓库的独立环境中执行过。静态报告只读取已保存产物，不重新读取 H5AD 计算。

### 6. 迁移和测试

已完成：

- 通用分析函数选择性迁移
- 函数签名与源代码核对
- 可选依赖延迟导入
- 输入不改写检查
- 不同 ID、不同数量、不同顺序检查
- 成员变化按 ID 对齐检查
- Level2 缺失检查
- 常量基因和无通路交集检查
- State/Gain 成功和不可计算分支检查
- 输出链接检查

目标仓库独立环境中：

- 39 项自动化测试通过
- 合成数据三个正式流程实际运行
- Notebook 实际执行
- 119 个报告内部链接和产物路径检查通过

示例结果见：

[合成 Impact 报告](../output/example/reconstruction_impact/report.html)

## 二、还没有完成

### 1. P1 真实样本尚未运行

这是当前最重要的未完成项。

P1CRC VisiumHD 的历史 Raw/SVC 文件虽然已经完成结构检查，但目前仍无法从已有元信息证明：

- `.X` 是未归一化的
- `.X` 是未取对数的
- 表达值可以直接作为本仓库约定的输入

因此目前没有：

- 将 P1 H5AD 复制到新仓库 `data/`
- 对 P1 运行三个正式分析
- 声称完成 P1 全量或抽样科学验收
- 根据 P1 结果做生物学结论

对应限制写在：

[迁移与验收](../docs/migration.md)

### 2. 尚未完成真实数据的科学验收

合成数据已经验证了：

- 调用链
- 参数传递
- 产物结构
- 失败和局部不可用处理

但它不能验证：

- P1 的真实空间结构
- P1 的 State/Gain 区域是否稳定
- 真实通路覆盖度
- 真实 Moran 分布
- 任何生物学改善结论

这些需要合规的 P1 输入后才能进行。

### 3. 浏览器交互式审阅还没有完整验收

已经完成静态 HTML 检查、链接检查和代表图视检，但浏览器自动化没有完整跑通。因此目前可以确认：

- 报告文件生成
- 内部链接有效
- 图表和表格存在

还没有把“浏览器中完整点击阅读体验”作为正式验收结论。

### 4. Xenium 尚未接入

Xenium 需要等待上游统一 SVC 和输入协议，目前没有迁入或适配。

### 5. 还没有做版本提交或发布

目标仓库已经有文件和独立环境，但尚未进行：

- Git commit
- tag
- 远程推送
- 版本发布

这不影响当前代码使用，但如果要形成可回溯交付版本，下一步需要单独整理提交。

## 三、下一步最短路径

现在只需要先确认 P1 的表达尺度：

> 能否确认历史 P1 Raw/SVC 的 `.X` 均为未归一化、未取对数的非负表达？

确认后再：

1. 复制 P1 文件到目标仓库；
2. 填写 P1 `sample.yaml`；
3. 用明确抽样规模运行三个流程；
4. 检查中间表、图、报告和 partial 项；
5. 再决定是否进行全量运行。
