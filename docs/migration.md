# 迁移与验收

## 来源与边界

通用分析来源：REVISE `main@c83dc97d25b6d513b59cc301255e5bdc7e9c7cd9`。实施时复核 `revise-2.0@a26d36fe3b11d652f2cd9eeafb92e990518eb3ad` 的通用分析内容相同。

AUCell 增强与 Impact 来源：`reconstruction-impact@e82dd13ce013f3f120475e79892906c97cb114d3`。2026-09-18 再次核对分支指向。选择性复制，源仓库保持不动，保留 MIT 许可。

不迁入 `revise.svc` 服务、重建 backend、Raw Level2 mapping、旧 batch 或 Notebook builder。新包没有 REVISE runtime 依赖。

## 已落地

| 类别 | 实现 | 验收依据 |
| --- | --- | --- |
| 架构与协议 | 独立 Python 包、样本与项目 YAML、公开 API、CLI、分层文档 | 输入原生轴和标签别名检查，输入 H5AD 不改写 |
| 计算工具 | 通用分析函数、Leiden、局部成员比较、窗口多样性、稳定阈值、原生 Moran、OmicVerse AUCell | 原实现回归与新非配对语义检查；可选提供者延迟导入 |
| 正式流程 | Impact、Moran、通路活性；顺序 sample × analysis 批量运行 | 单任务失败隔离；失败重跑发布当前失败记录；已有非生成目录受保护 |
| 结果审阅 | 结果索引、逐分析表图与静态报告 | 报告从已保存结果重绘；产物链接检查与代表图视检 |
| 科学工作台 | 连续 Impact Notebook，保留 `anatomy`、`partition`、`windows`、`state` 等中间变量 | 独立 Jupyter 内核实际执行合成样本，逐阶段检查 availability |
| 独立环境 | 目标目录 `.venv`，Python 3.11 与 `[dev,pathway]` | 不安装 REVISE；Scanpy 分群及 OmicVerse AUCell 实际执行 |

目标仓库独立环境自动化测试 **39 项通过**；三个正式流程与 12 个 Notebook 代码单元已在目标路径实际执行，119 个本地报告链接与产物路径检查通过。测试包含 ID 不相交/重排、成员按 ID 对齐、无 Level2、常量基因、无通路交集/缺少 provider、共同窗口缺失不补零、State/Gain 成功和不可计算分支。通用迁入函数的签名与来源 AST 核对一致，计算函数体除可选 Squidpy 延迟导入外与来源一致；11 项通用回归通过。部分流程测试替换昂贵提供者或受控阈值；它们用于检查组合语义，不作为真实算法或生物学验证。

## 比较口径改变

这次不是仅改包名。Raw/SVC 分群、Moran、通路评分独立运行，允许单位数与 ID 集不同。Moran 每侧原生 kNN 建图，默认 6 邻居；局部多样性共享物理窗口但独立单位与等量抽样。成员变化独立开关，使用 Raw 定义范围内的共同 ID；SVC 的独立抽样不缩小这一局部对应范围。与旧强制配对结果不要求数值相同。

State 是 SVC 父群自身的高多样性区域。Gain 只比较共同有效窗口的 SVC–Raw Leiden 多样性，正差异参与阈值选择。Raw Level2 是可选的独立局部 baseline，与 SVC Leiden 的窗口差异只作描述，不混入 Gain。具体参数与解释边界见[方法索引](analyses/README.md)。

## 合成验收与真实样本的区别

合成样本由 `scripts/create_example.py` 生成：Raw 240、SVC 210 个单位，600 个基因，部分共同 ID、两侧无 Level2。正式示例配置各侧独立抽样 180 个单位，显式 DEMO_PROGRAM；Moran 与通路流程成功，Impact 因缺少 Level2、窗口不足以稳定确定区域阈值等原因保持 `partial`。实际观察和参数记录在 `output/example/`。Notebook 执行版位于 `output/notebook/example/`，源码 Notebook 保持无输出。

这些结果验证调用、计算和产物组织，不能支持真实生物学结论，也不能称为 P1 或全量验收。浏览器自动化未成功启动；报告完成静态链接检查与图像视检，未宣称完整浏览器交互验收。

## 尚未完成：P1 真实输入验收

P1CRC VisiumHD 历史 Raw 与 SVC 合计约 8.8 GB：Raw 507,684 × 18,085，SVC 424,433 × 12,926。SVC 全部 ID 可对应 Raw，顺序不同；按 ID 对齐坐标一致；两侧无 Level2。

现有元信息和可追溯生成记录**不能证明这些历史文件的 `.X` 是未归一化、未取对数的线性表达**。当前生产代码的行为不能反向证明历史产物。未将它们复制为已核验输入，未对它们执行分析；不猜测、不静默逆变换。

下一步只需要上游确认该事实或交付合规 Raw/SVC，然后：

1. 复制至 `data/P1CRC_VisiumHD/`（H5AD 已忽略），从样本模板明确填写标签、坐标和尺度。
2. 使用项目模板的明确抽样规模运行三个流程，检查中间图与局部不可计算项。
3. 将实际结果与限制补到本页；抽样运行保持标注，只有实际运行过才声称全量验证。

Xenium 等待上游统一 SVC 后接入。
