# 本地输入

每个样本使用 `<sample>/sample.yaml`、`raw.h5ad`、`SVC.h5ad`，H5AD 不进入 Git。

P2CRC_Xenium 为真实标签空间验收：Raw 原样复制，Fibroblast 历史空间对象原样复制并命名为 SVC.h5ad。配置记录临时载体身份、源路径、SHA-256，以及未确认的表达身份/尺度。普通 loader 不分支处理它。正式单对象 SVC 到达后更新来源与能力，再运行受影响结果。不得把旧对象的结果当作新对象的验收。

P1 真实表达前提尚未确认；详见 [迁移验收](../docs/migration.md)。合成 example 只用于软件行为检查。
