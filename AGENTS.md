# 工作入口

- 仓库职责与依赖方向：`docs/architecture.md`。
- 文件、调用与结果协议：`docs/input-output.md`。
- 选择分析及方法定义：`docs/analyses/README.md`。
- 来源、验证和未完成工作：`docs/migration.md`。
- 新能力先复用 `methods/`；稳定流程进 `analyses/`；Notebook 保留可观察的连续主线。
- Raw/SVC 默认独立。仅成员变化分析局部配对；不得让它成为其他计算的前置条件。
- Web 只读已保存产物。缺少数据、资源或依赖时如实报告，不制造零值或成功。
- 不导入 `revise`，不改写输入 H5AD，不引入重建 backend 或旧载体兼容层。
- 保护用户和并发修改；变更后只运行与行为相关的验证。
