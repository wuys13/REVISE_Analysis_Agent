# 工作入口

- 公共工作区背景：[根目录 README](../../README.md)。
- 本地与远程协同协议：[协同协议](../../docs/collaboration.md)。
- 本仓库与 `REVISE` 必须保持在共同父目录下；配置中的 `../../REVISE/` 和交付清单中的
  `REVISE/`、`REVISE_Analysis_Agent/` 是当前本地与远程布局的相对入口。
- 仓库职责与依赖方向：`docs/architecture.md`。
- 文件、调用与结果协议：`docs/input-output.md`。
- 选择分析及方法定义：`docs/analyses/README.md`。
- 来源、验证和未完成工作：`docs/migration.md`。
- 新能力先复用 `methods/`；稳定流程进 `analyses/`；Notebook 保留可观察的连续主线。
- Raw/SVC 默认独立。仅成员变化分析局部配对；不得让它成为其他计算的前置条件。
- Web 只读已保存产物。缺少数据、资源或依赖时如实报告，不制造零值或成功。
- 不导入 `revise`，不改写输入 H5AD，不引入重建 backend 或旧载体兼容层。
- 保护用户和并发修改；变更后只运行与行为相关的验证。

## 接手与交付

- 按[根级接手规则](../../AGENTS.md)读取本仓状态入口 `docs/cross-repo-review/README.md`，核对 `sl`、输入提交和 `REVISE` 的三文件交付身份。
- 上游是并列的 `../REVISE`；结果报告由本仓 `data`、sample index 和 analysis result 登记，公共工作区只索引其入口。
- 结束时记录实际解释器、输入/输出摘要、结果校验、观测时间及 `partial/skipped/failed/unavailable` 状态；报告可供网页读取，但不改变输入对象。
