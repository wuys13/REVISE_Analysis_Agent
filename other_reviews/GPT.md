我检查了当前 `sl` 分支的 **`66b7fb4`**，对照了新仓库的输入代码、分析流程、完整主 Notebook、绘图与网页代码、测试，以及原 `reconstruction-impact` 中对应的方法说明。**你的感觉有依据：当前落实得比较好的是“独立分析仓库的运行骨架”，但你最看重的“保留已有科学逻辑、通过连续 Notebook 和中间可视化驱动开发”，只落实了一部分。**

我的总体判断是：

> **不需要推倒重做目录架构，但需要纠正几个边界，并把开发重心从“几个流程能运行、能生成文件”转回“真实样本上的完整分析链能被看懂、检查和继续迭代”。**

先说明审查范围：下面是源码与 Notebook 的静态审查；我没有在本轮重新运行 P1 数据或打开其实际生成的网页。仓库中“39 项测试通过”等是实施文档记录的结果，不是我本轮复跑的结果。该文档也明确承认真实 P1 验收、完整浏览器交互验收尚未完成。

---

## 一、哪些方向落实了，哪些没有达到要求

| 核心要求                                  | 当前落实情况                                               | 我的判断                        |
| ------------------------------------- | ---------------------------------------------------- | --------------------------- |
| 独立包，不依赖 REVISE reconstruction backend | 已形成独立包；Raw Level2 直接读取上游标签，不再调用 OT mapping           | **方向正确，保留**                 |
| 统一 Raw + SVC，不强制全局配对                  | 输入独立读取；成员变化作为局部可选比较                                  | **基本正确，保留**                 |
| 底层计算、流程、绘图、报告分层                       | 已有 `methods/`、`analyses/`、`plotting/`、`reporting/`   | **这套分层可以继续，不必改回我们示例中的目录名字** |
| 少数连续主 Notebook，不按算法拆碎                 | 确实只有一份主要 Reconstruction Impact Notebook，分阶段执行并保留中间变量 | **形式上落实，但科学串联和诊断能力不足**      |
| 继承旧 Reconstruction Impact 的空间分析逻辑     | anatomy、局部多样性、State/Gain 等已有；跨这些结果的定位、分层和解释链没有完整接上   | **核心缺口**                    |
| Notebook 调整后能进入正式 batch               | 两者存在参数、随机种子、输入范围和实现重复                                | **需要纠正**                    |
| 网页按科学问题组织，支持连续审阅                      | 有问题标题，但主要仍是按文件名分类的图表浏览器                              | **尚未达到要求**                  |
| 真实 H5AD 复制与真实案例验收                     | 文档明确记载 P1 尚未复制、尚未执行；Xenium 暂缓                        | **P1 未完成；Xenium 暂缓符合约定**    |

上述判断分别来自当前输入与主流程、包目录、Notebook 和报告实现，而不是仅看 README 的描述。

下面具体说问题。**我不会把“没有 Agent planner”“没有复杂 manifest”“没有旧 sc-SVC compatibility”算作缺点——这些本来就是我们决定暂时不做的。**

---

## 二、最明确的边界偏差：Input loader 新增了不必要的全局表达门槛

### 当前代码实际上要求什么

`revise_analysis/io.py::read_sample_config()` 中有一项硬检查：

```python
if config.get("expression", {}).get("scale") != "untransformed_nonnegative":
    raise ValueError(
        "Declare expression.scale: untransformed_nonnegative "
        "only for verified linear, unnormalized inputs"
    )
```

也就是说，**必须声明 Raw/SVC 都是未归一化的线性表达，整个样本才能加载**。这不是某个表达分析内部的要求，而是所有分析共同的入口要求。

这与我们最后确定的原则不一致：

> 输入层只负责读取两个独立对象；某项分析需要什么表达尺度、标签或空间信息，由该分析处理。

例如，只检查 full Raw 的 broad annotation 空间分布、生成 anatomy context，或检查已有标签的局部组成，**不应该先被整个 `.X` 的归一化历史阻断**。

### 这个新增门槛已经产生实际影响

迁移文档明确写道：P1 历史 Raw/SVC 的 `.X` 是否为未归一化、未取对数的表达尚未确认，因此**未复制成已核验输入，也未执行分析**。

这里需要分清两件事：

* **不猜测历史 `.X` 的语义，是对的。**
* **因此让所有分析和数据迁入都停住，不是我们要求的设计。**

正确调整不是强行给文件贴上“未归一化”标签，而是：

> **原样复制数据并记录当前已知信息；输入 loader 不强制一种表达尺度。需要 normalization/log/HVG 的具体计算，明确使用哪一个矩阵、是否还需要变换；不受表达尺度影响的检查和分析先正常进行。**

这只需要轻量修改，不需要再建立复杂 preflight 或新的大协议。

---

## 三、Notebook 最大的问题不是“数量”，而是“章节连续，但证据链没有完整连起来”

这一点是你最近反复强调的重点，我认为也是当前最值得纠正的地方。

### 1. 它确实没有被拆成一堆小 Notebook

当前 Notebook 有连续章节，也保留了：

```text
partition_results
anatomy_context
anatomy_windows
window_tables
state_region_tables
gain_tables
moran_tables
pathway_scores
```

这些中间对象。这部分是按照我们最后的方向做的，不能说完全理解错了。

**但“对象都出现了”不等于“对象之间的科学关系被展示出来”。**

### 2. 目前更像按顺序完成几个计算

现在的实际主线大致是：

```text
分群
→ anatomy
→ window diversity
→ State/Gain
→ Moran
→ AUCell
→ availability/status 汇总
```

Notebook 最后一节虽然叫 **Integrated observations and boundaries**，实际主要是在汇总：

```text
stage / status / reason / detail
```

以及“有多少个 summary”“保留了多少 Raw units”“Level2 是否存在”等执行信息。**这里不是跨结果的空间整合分析。**

而我们要保留的主线是：

```text
重建前后出现了什么差异？
        ↓
这些差异在空间上发生在哪里？
        ↓
它们位于 Tumor / Interface / Normal 的什么位置？
        ↓
局部多样性和 State Region 能不能解释这些位置？
        ↓
这些区域里是什么细胞状态、什么分子程序？
```

当前后半段还没有真正接起来。

### 3. 具体缺在哪里

| 当前已有结果                                | 还缺少的关键串联                                         | 为什么重要                 |
| ------------------------------------- | ------------------------------------------------ | --------------------- |
| `raw_anatomy_context`、anatomy windows | 把局部多样性和 State/Gain 结果按 anatomy 分层、比较             | 否则 anatomy 只是画过的一张背景图 |
| State/Gain window tables              | Region × anatomy 的空间叠加及分层范围、分母                   | 才能回答“高多样性区域主要在哪里”     |
| Raw Level2 baseline 与 SVC 局部多样性表      | 双 baseline 的并列阅读、空间分层，以及 Kobs/Neff/evenness 的对应图 | 不只是多保存一个 baseline CSV |
| Moran tables                          | 默认三个 parent 的结果组织，以及明确的比较口径                      | 当前主要展示全对象的两侧直方图       |
| AUCell scores                         | 与空间坐标、parent、anatomy/Region 的连接                  | 否则只知道分数分布，不知道程序在哪里    |
| 最后的 observations                      | 基于以上连接形成的事实总结与待判断问题                              | 才能支持你看结果后决定下一步改什么     |

这些缺口可以从当前主流程和 Notebook 看出来：anatomy 单独写出；State/Gain 单独写出；Moran 和 pathway 主要按 `raw/svc` 计算，通路图是分数直方图；最后并没有把它们组合成区域分层结果。

这里也要注意：**我不是要求马上新增完整 Treg/TLS/CCI 分析。**先恢复旧 Reconstruction Impact 已经具备的 anatomy 分层、局部多样性、Region 和 EMT/program localization 串联，就已经能显著改善主 Notebook。Subtype-enriched region 的新方法，仍可在对应主线位置继续探索，不必为了目录完整而仓促实现。

---

## 四、可视化诊断没有充分继承：有些能力已经迁了，却没有接入主线

### 最典型的例子：`select_window_scale`

新仓库的 `methods/regions.py` **已经存在** `select_window_scale()`，能够返回候选窗口的支持情况和尺度选择结果。

但是正式 Impact 流程和 Notebook 当前主要使用一个固定窗口：

* Notebook 默认 **32 μm**；
* batch API 默认 **100 μm**。

方法文档也明确写出了这个区别。

这说明问题不是“没写这个函数”，而是：

> **迁入了尺度选择能力，却没有把“候选尺度 → 支持曲线 → 当前采用尺度 → 空间表现”接回你要看的连续开发主线。**

原分支中，窗口支持曲线、尺度选择、双 baseline、anatomy 分层、threshold reliability 和空间 mask，本来就是分析组织的一部分，而不是可有可无的装饰。

### State threshold 的检查也存在类似不足

当前 Notebook 会保存 threshold 和 bootstrap 表，也会在阈值有效时画 State Region；但没有完整呈现“分布/阈值诊断 → 连续场 → mask → anatomy 关系”这一套检查链。尤其阈值不可用时，parent State 图会被跳过；前面的 local diversity 图主要展示 `All`。

这对你的开发方式很不方便。因为你最需要看图的时候，往往恰好是：

> 为什么这个 parent 没有稳定 threshold？是窗口支持不足、分布形态不合适，还是当前分析范围不合理？

**只给 `no_stable_threshold`，不把帮助判断原因的关键图放在附近，不能满足“结果驱动调参”的要求。**

因此下一轮不该先补更多算法，而应把已有的尺度选择、支持情况、连续场和阈值诊断接到对应 Notebook 小节中。

---

## 五、Notebook 与 batch 已经发生实质性漂移，不只是展开方式不同

我们原来的要求是：

> Notebook 分阶段展开，batch 一次执行；但稳定阶段应共享科学计算和有效参数。

当前不是完全这样。

| 项目                                  | Notebook                         | batch Impact                           |
| ----------------------------------- | -------------------------------- | -------------------------------------- |
| 窗口默认值                               | 32 μm                            | 100 μm                                 |
| SVC 分群 seed                         | 默认 43                            | 分群调用使用 `random_state`，默认 42            |
| 局部 membership 的 SVC 重分群 seed        | 默认 44                            | 默认 42                                  |
| State/Gain threshold bootstrap seed | 使用 SVC window seed，默认 43         | 默认 42                                  |
| Moran 的数据范围                         | 直接使用完整 `sample.raw / sample.svc` | 使用已经过 `sample_n_units` 抽样的 `side_data` |
| 绘图实现                                | 大量图在 Notebook 内重新写               | 另用 `plotting/` 中的函数                    |

这些差异能直接在当前 Notebook 和正式流程中对应到。

**不同 seed 本身不是错误，不同分析范围也不一定不合理。问题在于：你在 Notebook 看好、调好的结果，没有一个可靠方式直接成为正式 batch 的同一个结果。**

这会破坏整个开发循环：

```text
Notebook 看结果、调参数
        ↓
批量运行
        ↓
结果却变了
        ↓
无法判断是参数、抽样、算法还是实现不同
```

### 应该怎么改，但不增加复杂框架

只需要做三件事：

1. **一份有效参数。**Notebook 读入后允许局部修改，但正式运行要使用这份修改后的参数，而不是另一组默认值。
2. **稳定阶段共用函数。**例如形成双 baseline 对比表、State/Gain 表的稳定逻辑，不要在 Notebook 和 batch 中维护两份。
3. **Notebook 保留展开。**仍然明确看到 anatomy、windows、state、program 等对象和短调用；不是把整本 Notebook 压成一个 `run()`。

这不是要求两者执行形式相同，而是要求：

> **同一阶段、同一输入范围、同一参数，计算含义和结果应相同。**

---

## 六、网页目前仍偏“按文件分类浏览”，没有真正成为科学阅读主线

当前 `reporting/report.py` 确实有四个科学问题标题，而且不读取 H5AD，这两点是好的。

但实际分组逻辑是 `_group(relative)`：

```text
文件名包含 membership → membership
包含 moran/pathway/program → genes
包含 anatomy/window/state/gain → spatial
其余 → overall
```

每个大节再先集中放图，其他表格统一收进 “Saved evidence”。

所以它目前更像：

> **带科学问题标题的文件浏览器。**

不是我们想要的：

```text
这个问题是什么
↓
用了什么空间定义
↓
关键定义是否合理
↓
对应 Raw/SVC 结果
↓
按 parent/anatomy 进一步展开
↓
完整表格入口
```

### 一个很具体的例子

现在只要文件名含有 `anatomy`、`state`、`gain`，就会被放进同一个 spatial 大节。

但你真正需要的是：

```text
空间定位差异
├── 组织背景：Tumor / Interface / Normal
├── 局部状态：Raw Leiden / Raw Level2 / SVC
├── State：连续场、可靠性、mask
├── State 在不同 anatomy 中的分布
└── 程序在这些区域中的定位
```

这两种组织方式的差别，不是 CSS 好不好看，而是**读者是否需要自己在几十个文件间重新拼出科学逻辑**。

### 下一轮应如何修正

不需要恢复旧 2000 多行 renderer，也不需要通用组件系统。只需为 Reconstruction Impact 写一个小而明确的阅读结构：

* 某个小问题展示哪些已保存的表和图；
* 三个 parent 如何并列或展开；
* 支持诊断放在哪里；
* 哪些是正式结果，哪些只是辅助判断；
* 没有结果时在原位置说明原因。

**轻量化应该减少实现复杂度，而不是减少科学组织。**

---

## 七、复制迁移与科学改造混在一起，保真验收因此被弱化了

这一点需要特别注意。

迁移文档明确写了：

> “这次不是仅改包名。”

并列出了独立分群、每侧 native kNN Moran、独立抽样和窗口比较等改变，同时说明与旧强制配对结果不要求数值相同。

其中一部分改变确实是你要求的，例如：

* 不再全局要求 Raw/SVC 对应；
* 不在分析侧做 Raw Level2 mapping；
* 不使用旧双载体输入。

**但这些必要改动，并不自动授权把其他已有科学处理一起替换。**

现在至少需要分别核对：

| 改动类型                             | 可以怎样处理                  |
| -------------------------------- | ----------------------- |
| namespace、文件路径、去除 REVISE backend | 属于迁移，原则上应保持计算行为         |
| 双载体改成完整 SVC 输入                   | 属于必要适配，应明确新旧表达对象的关系     |
| 不完整对应时增加分布级分析                    | 属于新支持范围，单独验证            |
| 特征选择、QC、抽样范围、空间图、窗口尺度改变          | 属于科学分析口径变化，不能统一算作“解耦”   |
| 旧方法本身存在问题                        | 单独记录、单独修复，不应静默改，也不应盲目保留 |

例如当前 `prepare_side()` 是**先对整侧抽样，再在抽样结果中切 parent**；随后这些对象又用于局部窗口多样性。这种范围选择会影响 parent 数量和窗口支持，不能与旧方案的 scope-specific 分析当作同一口径。

同样，**允许在不同空间图上分别计算 Moran，不代表这些数值自然就成为同口径的“恢复提升”比较**。当前实现可以保留为 native spatial description，但需要区分它与同一空间尺度或明确对应条件下的比较；不应因此把全局 pairing 加回来。

我也需要补正前面交接中的一个容易导致偏移的地方：

> **“统一输入、不强制配对”不等于“所有分析统一改成两侧独立计算”；“删除旧 orchestration”也不等于“删去旧分析的条件、诊断和串联关系”。**

下一轮应恢复一张很小的“原能力 → 当前实现 → 保留/必要适配/科学改动/暂缓”核对表，而不是再写一份庞大的新架构文档。

---

## 八、有几个具体实现问题，应单独修正

### 1. 高变基因选择的输入顺序不符合所选方法要求

当前 `methods/partition.py::compute_partition()` 是：

```python
normalize_total(...)
log1p(...)
highly_variable_genes(
    ...,
    flavor="seurat_v3",
    check_values=False,
)
```

Scanpy 官方文档明确指出，`seurat_v3` 这一 flavor 预期的是 count 数据，而不是已经 log 的数据。关闭 `check_values` 不会改变算法的输入要求。([scanpy.readthedocs.io][1])

**但这不是此次迁移才新引入的错误：我回查原分支，旧 `select_raw_hvg_feature_names()` 也存在同样顺序。**因此应把它列为一个继承下来的方法问题，单独确定修复方式和结果影响，不能全部归因于新仓库实现。

这也说明为什么“新旧一致”与“方法正确”要分开验收。

### 2. Notebook 在指定样本路径不存在时，可能静默切回合成示例

开头存在这样的逻辑：

```python
if not SAMPLE_YAML.exists() and Path("data/example/sample.yaml").exists():
    SAMPLE_YAML = Path("data/example/sample.yaml")
```

这没有区分“默认路径解析失败”和“用户明确指定的真实样本不存在”。在相应工作目录下，用户请求的路径错误可能被替换成 example。默认输出目录也仍指向 `output/notebook/example`，没有自然绑定实际 sample ID。

这两项应直接修正：

> **用户明确指定的样本不存在就报错；默认输出目录从实际加载的 sample ID 派生。**

这是非常必要的小检查，不是我们反对的复杂 preflight。

### 3. 大量 `except Exception → unavailable` 不利于开发调试

Notebook 的分群、membership、anatomy、window 等阶段多处捕获所有异常，再继续记录为 `unavailable`。

对可预期的缺标签、缺资源，这样处理合理；但代码错误也被归成同一种“不可用”，会削弱你想要的调试能力。

应区别：

> 缺少可选输入，可以跳过；实现错误，应让你看到 traceback，而不是看起来只是“这个样本不支持”。

### 4. Impact 与独立 Moran/pathway 任务存在重复计算

Impact 内部已经调用 Moran 和 pathway 计算；示例项目又把 `spatial_autocorrelation`、`pathway_activity` 作为独立分析再运行一遍。

这不是必须马上建立 DAG/cache 的理由。第一版只要明确一种默认运行方式：**主系列已经算过的结果，报告和后续章节直接消费；不要默认以两个入口重复算同一份证据。**

---

## 九、为什么测试通过了，你仍然觉得“不像我们要的东西”

因为当前验收重点与产品要求之间有落差。

仓库文档记录了合成样本运行、39 项测试和链接检查，同时明确承认：

* P1 真实数据尚未复制和执行；
* 合成数据不能当作真实生物学验收；
* 浏览器自动化没有成功启动；
* 部分流程测试使用了替代提供者或受控阈值。

我检查的测试也符合这个定位。例如 State/Gain 输出测试会把分群和 threshold 函数替换成受控结果，用来验证文件和 mask 能否写出；报告测试主要检查 HTML 是否生成、是否包含预期文本。**这些测试有价值，但验证的是调用和产物组织，不是你所要求的完整科学阅读与调试体验。**

因此，目前应该区分三种完成度：

| 验收层                | 当前证据能支持什么                        |
| ------------------ | -------------------------------- |
| 软件骨架验收             | 已有相当进展：独立运行、输入输出、部分错误隔离和测试       |
| 复制迁移保真验收           | 尚不足：有科学口径改变，真实案例对照未完成            |
| Notebook/Web 工作台验收 | 尚不足：诊断链、跨结果串联、批量一致性和真实阅读效果没有完整验证 |

**不是应该停止写测试，而是不能再把“文件产生了、测试通过了”当作主要完成标准。**

---

## 十、我建议下一轮给 Codex 的整改顺序

我不会让它重新设计仓库。建议保留当前分层，集中完成下面几项：

| 顺序 | 修改重点                                       | 主要位置                                         | 验收方式                                  |
| -- | ------------------------------------------ | -------------------------------------------- | ------------------------------------- |
| 1  | 取消统一线性表达的 loader 硬门槛；表达需求回到具体计算；修正显式样本路径回退 | `io.py`、Notebook 开头                          | 真实样本可读取；未知表达语义不被伪造；不相关分析不被阻塞          |
| 2  | 接入真实 P1 工作数据，先做明确范围的运行                     | 数据目录、sample/project 配置                       | 真实输入已复制；哪些部分运行、哪些受限均有实际依据             |
| 3  | 统一 Notebook/batch 有效参数和稳定阶段实现              | Notebook、`analyses/reconstruction_impact.py` | 相同输入范围与参数下，关键中间表、Region 和分数一致         |
| 4  | 恢复完整空间主线及必要诊断                              | `methods/regions.py`、绘图、主 Notebook           | 能连续检查尺度支持、局部多样性、State、anatomy 分层和程序定位 |
| 5  | 把网页从文件分类改成明确的科学阅读结构                        | `reporting/report.py`                        | 不靠读者猜文件名，就能沿同一问题下钻到对应表图               |
| 6  | 分开记录迁移改动与科学方法修复，包括 HVG 问题                  | `methods/partition.py`、迁移说明、测试               | 保留项做等价核对；改变项明确理由、影响和验证，不混为一谈          |

其中第 4 项的最小目标不是“把所有未来分析补完”，而是把这一条做通：

```text
组织背景与分析范围
        ↓
窗口支持与尺度选择
        ↓
Raw Leiden / Raw Level2 / SVC 的局部状态
        ↓
State 连续场、阈值可靠性和 Region
        ↓
Region 在 Tumor / Interface / Normal 中的位置与分布
        ↓
对应 parent 的程序空间定位
        ↓
当前能支持什么、还需要检查什么
```

这条链可以仍然放在**同一个主 Notebook**中。代码继续模块化；参数探索就在对应小节中；网页按阅读顺序重新组织同一批结果。完全不需要再拆很多 Notebook，也不需要加 Agent planner。

---

### 最后的判断

**当前不是“架构全部错了”，而是“技术骨架落实得比科学工作台更完整”。**

最需要纠正的三个理解是：

> **简洁架构，不等于删减已经形成的科学分析逻辑。**
> **一个连续 Notebook，不等于把几个独立计算依次排在一个文件里。**
> **轻量网页，不等于把图表按文件名装进几个大章节。**


可以。基于当前仓库 `sl` 分支的实际实现，我认为下一轮不应该再做“架构扩建”，而应该进入一个很明确的阶段：

> **把 Reconstruction Impact 这条主线真正做完整，并以它为样板，把“模块 → 连续 Notebook → batch → Web”闭环打通；然后再以同样模式扩展 Recoverability 和 biological follow-up。**

下面这份可以直接作为下一轮 Codex 的整改与开发设计。

# REVISE_Analysis_Agent 下一阶段完整整改与开发设计

## 0. 下一阶段的核心目标

当前仓库已经搭出了基本骨架：

```text
Raw + SVC
   ↓
methods/
   ↓
analyses/
   ↓
batch
   ↓
tables / figures / result.json
   ↓
Notebook / Web
```

但当前主要完成的是：

> **“分析可以被调用并产生文件”**

下一阶段要提升到：

> **“一条完整科学分析主线可以在真实样本上被连续观察、调参、验证、批量执行，并在网页中按照同样的科学逻辑被阅读。”**

因此下一阶段不要继续优先增加：

* 更多目录；
* 更多 schema；
* Agent planner；
* trigger graph；
* 复杂 provenance；
* 更多 synthetic demo。

优先把已有核心能力真正**串起来**。

---

# 一、下一阶段需要形成的三个层次

最终开发模式固定为：

```text
                    reusable methods
                           │
                           ▼
                  scientific workflow
                           │
              ┌────────────┴────────────┐
              │                         │
              ▼                         ▼
     Canonical Notebook              Batch
 连续分析、调参、可视化诊断       稳定批量执行
              │                         │
              └────────────┬────────────┘
                           ▼
                    saved results
                           │
                           ▼
                       Web report
                  连续科学阅读与审阅
```

关键原则：

### Methods

回答：

> 这个具体东西怎么算？

### Notebook

回答：

> 为什么现在要算这个？
> 它依赖前面什么？
> 结果长什么样？
> 参数合理吗？
> 下一步为什么进入另一个分析？

### Batch

回答：

> 已经稳定的这条分析链如何规模化运行？

### Web

回答：

> 已经稳定的分析结果如何按照科学问题连续阅读？

---

# 二、当前最优先需要修正的基础问题

这些应先完成，否则后面的主线会不断漂移。

---

## 2.1 去掉 Input loader 对统一表达尺度的硬门槛

当前：

```python
expression.scale == "untransformed_nonnegative"
```

是整个 sample loader 的硬要求。

这过强。

应该改成：

```text
loader:
    只加载 Raw / SVC
    不修改输入
    不替用户推断表达语义
```

具体 method 如果要求某种表达：

```text
partition
AUCell
DEG
...
```

由 method 自己处理或拒绝。

### sample.yaml

可以保留：

```yaml
expression:
  scale: ...
```

作为 metadata。

但：

> 不应该阻止与表达尺度无关的分析运行。

---

# 三、Notebook 与 batch 必须统一参数来源

当前存在：

```text
Notebook window = 32 µm
batch default = 100 µm
project template = 40 µm
```

以及 seed、sampling 等不同。

这会直接破坏：

```text
Notebook 调好
↓
batch 跑出来却不是同一个结果
```

---

## 3.1 以后正式参数只有一个来源

推荐：

```text
project.yaml
```

作为 batch 默认参数。

Notebook 开始时：

```python
params = load_analysis_parameters(
    project_yaml,
    sample_id,
    "reconstruction_impact",
)
```

然后：

```python
params
```

直接显示。

---

## 3.2 Notebook 允许临时 override

例如：

```python
PARAMETER_OVERRIDES = {
    "window_side_microns": 32,
}
```

得到：

```python
working_params
```

但 Notebook 应明确显示：

```text
formal config: 40 µm
current notebook: 32 µm
```

最好提供非常轻量的：

```python
show_parameter_diff(base_params, working_params)
```

这样你一眼知道自己正在试什么。

---

## 3.3 Notebook 验证成功后

再把：

```text
32 µm
```

写回 `project.yaml`。

不需要复杂自动同步机制。

---

# 四、sampling 逻辑需要修正

当前：

```text
先对整侧 sample_n_units
↓
再从这个抽样结果里取 Fibroblast / T / Mono
```

容易导致：

```text
某些 parent 被不必要缩小
```

尤其不同侧 composition 不同时，影响局部空间分析。

---

## 4.1 应改为 scope-aware sampling

更合理：

```text
Raw full object
├─ All → independently sample up to N
├─ Fibroblast → independently sample up to N
├─ Mono/Macro → independently sample up to N
└─ T → independently sample up to N

SVC 同理
```

也就是说：

> sample limit 属于分析 scope，不属于整个 object 的全局预过滤。

---

## 4.2 full tissue context 不要抽样

例如：

```text
anatomy
Tumor / Normal / Interface
```

应始终使用：

```text
full Raw Level1 context
```

因为它是组织背景。

不要因为 partition sampling 改变 anatomy。

---

# 五、错误处理也需要调整

开发 Notebook 最重要的是：

> 出 bug 时让开发者真正看到 bug。

因此：

### 可以捕获

```text
Level2 缺失
pathway resource 缺失
某个 scope 不存在
支持不足
```

这些是预期 availability condition。

### 不应吞掉

```text
AttributeError
unexpected indexing bug
column typo
shape bug
programming error
```

Notebook 中这类错误应该直接 traceback。

不要统一变成：

```text
unavailable
```

否则结果驱动开发会非常困难。

---

# 六、先解决 Partition 的方法学问题

当前：

```text
normalize
→ log1p
→ highly_variable_genes(flavor="seurat_v3")
```

这需要重新审视。

而且这是旧实现继承的问题，不应当在迁移中悄悄改变。

---

## 6.1 下一步应该做一个小型方法核对

至少比较：

### Option A

```text
linear expression
→ seurat_v3 HVG
→ normalize/log
→ PCA/graph
```

### Option B

使用适合 log-expression 的 feature selection。

然后在真实 P1 上看：

```text
cluster stability
cluster number
spatial coherence
Raw/SVC representation
```

再确定正式方法。

---

## 6.2 不要同时混入其他修改

这应该是一个独立 methodological correction。

记录：

```text
旧方法是什么
新方法是什么
为什么改
结果变化多少
```

而不是藏在迁移里。

---

# 七、第一条必须完整做完的主线：Reconstruction Impact

这是当前最高优先级。

它应该真正回答两个总问题：

> **Reconstruction changed what?**

以及：

> **Where did those changes occur?**

所以主线应该形成：

```text
Representation
      ↓
Spatial context
      ↓
Local state structure
      ↓
Regions
      ↓
Gene / program localization
      ↓
Integrated spatial interpretation
```

---

# 八、主线 1：Representation change

回答：

> Raw 与 SVC 能解析出的 cellular/state representation 有什么变化？

---

## 8.1 Broad composition overview

增加简单：

```text
Level1
Raw count / fraction
SVC count / fraction
```

输出：

```text
representation/broad_composition.csv
```

这不是为了说谁更好。

只是先明确：

> 两个对象总体 population structure 是否发生明显变化。

---

## 8.2 Independent partition

继续保留当前：

```text
Raw independently cluster
SVC independently cluster
```

每个：

```text
All
Fibroblast
Mono/Macro
T
```

输出：

```text
n units
n features
n clusters
cluster-size distribution
```

---

## 8.3 Same-resolution comparison

固定相同：

```text
Leiden resolution
```

观察：

```text
Raw K
vs
SVC K
```

这是：

> representation complexity diagnostic。

---

## 8.4 Membership change 仍然是 optional paired analysis

只有：

```text
真实共享 unit ID
```

时才运行。

不要让它成为 Reconstruction Impact 的入口。

---

## 8.5 matched-K 继续作为 control

它解决：

> 如果强制相同 cluster 数，membership difference 是否仍然存在？

保留为 secondary diagnostic。

---

## 8.6 后续建议增加一个 controlled-feature representation check

不是第一轮阻塞项，但值得补。

现在 independent partition：

```text
Raw 用 Raw 自己的 features
SVC 用 SVC 自己的 features
```

它回答：

> 每个 representation 自己最大限度能解析什么。

后续再增加：

```text
same common feature space
```

作为 control：

> representation 差异是否只是 feature selection 导致？

不要取代 native partition。

两者回答不同问题。

---

# 九、主线 2：建立 Spatial Context

这是当前真正没有串完整的关键部分。

在进入：

```text
State Region
EMT
Treg
Moran
```

之前，必须先建立一个统一空间背景。

---

## 9.1 Full tissue anatomy

继续使用完整 Raw Level1。

构建：

```text
Tumor
Normal
Interface
Other
```

注意：

```text
Interface
```

应该是 window context：

```text
Tumor source + Normal source
```

共同存在的窗口。

---

# 十、Window scale 不应该再只是一个 hidden parameter

它应该成为 Notebook 中明确的一节。

当前已经有：

```python
select_window_scale()
```

但没有真正接进主线。

---

## 10.1 增加候选 window sizes

例如：

```yaml
candidate_window_sides_microns:
  - 16
  - 24
  - 32
  - 40
  - 56
  - 80
```

---

## 10.2 Notebook 展示

每个候选尺度至少看：

```text
valid window fraction
retained unit fraction
n valid windows
```

按：

```text
Raw
SVC
Fibroblast
Mono/Macro
T
```

必要时并排。

---

## 10.3 `select_window_scale()` 只给 recommendation

非常重要：

> 不要自动把推荐结果隐藏起来直接进入后面。

Notebook 应该显示：

```text
recommended scale: 40 µm
current configured scale: 32 µm
```

由你看图决定是否修改。

这正符合“可视化驱动开发”。

---

# 十一、Spatial Context Notebook 中至少要有这些图

### A. Raw tissue map

```text
Level1
```

### B. Anatomy window map

```text
Tumor
Normal
Interface
Other
```

### C. Window support curves

### D. Selected grid overlay

至少检查：

> window size 和 tissue geometry 是否合理。

---

# 十二、主线 3：Local state structure

这是 Reconstruction Impact 的核心之一。

需要明确形成三个 local-state representations：

```text
Raw Leiden
Raw Level2
SVC Leiden
```

---

## 12.1 为什么是三套

### Raw Leiden

回答：

> Raw 自己从表达中解析出的 local state complexity。

### Raw Level2

回答：

> 上游 annotation 已经给定的 biological subtype baseline。

### SVC Leiden

回答：

> Reconstruction 后能解析出的 local state complexity。

三者不能混成一个指标。

---

# 十三、Local diversity 继续保留 4 个指标

每个 window：

```text
Kobs
Entropy
Neff
Evenness
```

并保留：

```text
n_units
valid_window
```

---

# 十四、恢复真正的“双 baseline”阅读

当前已经计算：

```text
Raw Level2
vs
SVC Leiden
```

但没有充分组织。

需要生成一张统一 common-window table：

```text
window_id
anatomy
raw_leiden_neff
raw_level2_neff
svc_neff
delta_vs_raw_leiden
delta_vs_raw_level2
...
```

同样包括：

```text
Kobs
entropy
evenness
```

---

## 14.1 这张表很重要

它以后既服务：

```text
Notebook
Web
Agent
```

也使你不用在几个 CSV 之间自己 join。

---

# 十五、Local state 要按 Anatomy 分层

这是当前最明显缺失的连接之一。

增加：

```text
diversity_by_anatomy
```

例如：

| scope      | anatomy   | representation | Neff median | IQR | n windows |
| ---------- | --------- | -------------- | ----------: | --: | --------: |
| Fibroblast | Tumor     | Raw Leiden     |         ... | ... |       ... |
| Fibroblast | Interface | SVC            |         ... | ... |       ... |

这样才能回答：

> SVC 增加的复杂状态主要在哪里？

而不是只知道：

```text
overall Neff higher
```

---

# 十六、主线 4：State Region

这是当前必须完整恢复的一条链。

不要直接：

```text
Neff
→ threshold
→ mask
```

Notebook 应明确展示：

```text
SVC local Neff continuous field
         ↓
Neff distribution
         ↓
threshold breakpoint
         ↓
bootstrap threshold distribution
         ↓
CI / stability
         ↓
Region mask
         ↓
Region extent
         ↓
Anatomy localization
```

---

# 十七、State threshold 必须重新成为可视化诊断

增加：

```text
threshold point
bootstrap distribution
CI
valid bootstrap fraction
observed range
```

如果：

```text
no_stable_threshold
```

仍然显示：

```text
continuous field
distribution
bootstrap result
```

而不是：

```text
没有 region → 什么都不画
```

因为恰恰这个时候最需要判断：

> 为什么没有稳定 Region？

---

# 十八、State Region 必须进一步按 anatomy 分解

增加：

```text
State Region extent by anatomy
```

例如：

| scope | anatomy   | valid windows | state windows | fraction |
| ----- | --------- | ------------: | ------------: | -------: |
| T     | Interface |            52 |            31 |     0.60 |
| T     | Tumor     |           140 |            20 |     0.14 |

这才真正回答：

> High-diversity State 在哪里？

---

# 十九、Gain Region 保留，但定位要明确

继续坚持旧设计：

```text
State Region
≠
Gain Region
```

### State

回答：

> SVC 中哪里存在高度复杂的 reconstructed state。

### Gain

回答：

> 哪些共同空间窗口的 SVC local diversity 高于 Raw。

Gain 应作为：

```text
audit / supporting view
```

而不是 Reconstruction Impact 的唯一 headline。

---

# 二十、增加 subtype-enriched spatial field

这是下一阶段非常重要的扩展。

例如：

```text
Treg-enriched
Mreg-enriched
CAF5-enriched
```

但不要一开始设计复杂 composite score。

---

## 20.1 推荐方式

先计算：

```text
每个 window:
    subtype count
    subtype fraction
    global/parent background fraction
    enrichment
```

例如：

```text
T parent:
Treg fraction in window
/
Treg fraction in all T
```

---

## 20.2 先形成独立 enrichment field

不要直接定义：

```text
High-diversity-Treg-Interface super score
```

而是分别拥有：

```text
State mask
Treg enrichment field
Anatomy class
```

然后 Notebook 做：

```text
State
∩
Treg enriched
∩
Interface
```

这比构造一个无法解释的综合分数更合理。

---

## 20.3 subtype enrichment threshold

第一阶段可以先在 Notebook 探索。

等真实样本上确认：

```text
什么 support
什么 enrichment threshold
```

合理以后再固化进 batch。

这正是 Notebook 应该发挥作用的地方。

---

# 二十一、主线 5：Gene spatial structure / Moran

当前 Moran 太孤立，而且只跑：

```text
Raw whole side
SVC whole side
```

需要扩展成：

```text
All
Fibroblast
Mono/Macro
T
```

分别计算。

---

# 二十二、Moran 应区分两种阅读

## A. Native-side Moran

Raw：

```text
Raw coordinates
Raw native graph
```

SVC：

```text
SVC coordinates
SVC native graph
```

回答：

> 各自 representation 中 gene spatial structure 是什么。

这是 descriptive。

---

## B. Shared-gene comparison

对于两边都存在的 gene：

```text
gene
raw_moran
svc_moran
delta
```

形成 scatter / distribution。

但必须明确：

> 这是同 gene 的 cross-representation descriptor，不代表 paired spatial reconstruction error。

---

## 22.1 输出至少有

```text
moran_raw_<scope>.csv
moran_svc_<scope>.csv
moran_shared_<scope>.csv
```

---

## 22.2 Notebook 图

至少：

```text
Raw vs SVC Moran scatter
Moran distribution
top positive delta genes
top negative delta genes
```

然后可以选择代表 gene 画空间图。

---

# 二十三、主线 6：Pathway / program

当前 AUCell 主要是：

```text
score distribution
```

还不够。

真正重要的是：

> **program 在哪里发生。**

---

## 23.1 对每个 parent 分开计算

例如：

```text
Fibroblast EMT
Mono/Macro inflammatory / suppressive program
T cell state programs
```

第一阶段可以仍然从 EMT 开始。

---

## 23.2 保存 unit-level spatial field

至少：

```text
side
scope
unit_id
x
y
program
score
```

而不是只有一列 score。

---

## 23.3 再聚合到共同 physical windows

分别：

```text
Raw unit scores
→ windows

SVC unit scores
→ same windows
```

可以得到：

```text
raw_window_score
svc_window_score
delta
```

完全不需要 ST-unit pairing。

---

# 二十四、Pathway 要进一步与 Anatomy / Region 结合

增加：

```text
program_by_anatomy
```

例如：

```text
EMT:
Tumor
Interface
Normal
```

以及：

```text
program_in_state_region
```

例如：

> EMT 是否集中于 Fibroblast high-diversity State Region？

这才把：

```text
program
```

接回：

```text
spatial localization
```

主线。

---

# 二十五、主线 7：Integrated Spatial Localization

这一节是当前最缺的。

它不是新算法，而是：

> **把前面已经计算的空间事实组合起来。**

应该至少形成：

```text
Anatomy
    ×
State Region
    ×
Subtype enrichment
    ×
Program field
```

必要时：

```text
× membership change
```

---

# 二十六、Integrated section 要产生几张标准事实表

例如：

## A. State × Anatomy

```text
scope
anatomy
state_window_fraction
```

## B. Subtype × Anatomy

```text
scope
subtype
anatomy
enrichment
```

## C. State × Subtype

```text
scope
subtype
state_vs_nonstate enrichment
```

## D. Program × State

```text
scope
program
score_in_state
score_outside_state
```

## E. Program × Anatomy

```text
scope
program
anatomy
score
```

这些还不是 biological conclusion。

它们是：

> 可以进一步解释的事实证据。

---

# 二十七、Notebook 最终应形成下面这一条真正连续的故事

```text
1. Raw 与 SVC 是什么？
        ↓
2. Reconstruction 改变了能解析出的 cellular representation 吗？
        ↓
3. 要谈“在哪里”，先建立怎样的 tissue / window spatial context？
        ↓
4. Raw 和 SVC 的 local state complexity 在哪里不同？
        ↓
5. Reconstruction 后哪里形成 high-diversity State？
        ↓
6. 这些 State 位于 Tumor、Interface 还是 Normal？
        ↓
7. 哪些 subtype 在这些区域富集？
        ↓
8. 哪些 gene/program 的空间组织发生变化？
        ↓
9. program 是否与 State / subtype / anatomy 共定位？
        ↓
10. 当前事实支持我们进一步追问什么？
```

这个才是：

```text
01_reconstruction_impact.ipynb
```

最终要呈现的主线。

---

# 二十八、Notebook 的建议详细结构

## Section 0 — Configuration

显示：

```text
sample
formal params
notebook overrides
parameter diff
```

---

## Section 1 — Input overview

显示：

```text
Raw n
SVC n
genes
Level1 composition
coordinates
```

不做复杂 audit。

---

## Section 2 — Representation

显示：

```text
broad composition
Raw partitions
SVC partitions
cluster numbers
cluster size
optional membership
optional matched-K
```

---

## Section 3 — Spatial framework

显示：

```text
candidate window support
recommended scale
selected scale
Raw tissue map
anatomy map
selected grid
```

这一节之后才正式进入：

```text
local spatial analyses
```

---

## Section 4 — Local state

针对：

```text
Fibroblast
Mono/Macro
T
```

每个看：

```text
Raw Leiden
Raw Level2
SVC Leiden
```

再展示：

```text
Kobs
Neff
entropy
evenness
```

以及：

```text
by anatomy
```

---

## Section 5 — State / Gain Regions

每个 parent：

```text
continuous SVC Neff
threshold diagnostic
State mask
State extent
State × anatomy
```

Gain：

```text
common-window delta
Gain mask
Gain extent
```

放次级。

---

## Section 6 — Subtype enrichment

例如：

```text
Treg
Mreg
CAF subtype
```

显示：

```text
continuous enrichment field
```

再和：

```text
State
Anatomy
```

叠加。

---

## Section 7 — Gene spatial structure

Moran：

```text
Raw
SVC
shared genes
delta
representative genes
```

---

## Section 8 — Programs

例如 EMT：

```text
coverage
Raw/SVC distribution
Raw/SVC spatial field
window-level field
program × anatomy
program × State
```

---

## Section 9 — Integrated spatial picture

将核心 Region 叠在同一个逻辑图中：

```text
Anatomy
State
Subtype
Program
```

不是必须一张巨大图。

可以是：

```text
2×2 aligned panels
```

---

## Section 10 — Current factual observations

只总结：

```text
已经看到什么
哪些条件成立
哪些还不确定
```

不要自动写 biological causality。

并列出：

```text
下一步值得进一步分析的问题
```

这部分以后自然可被 Agent 消费。

---

# 二十九、Notebook 不应该把整个 workflow 包成一个 run()

Batch 可以：

```python
run_reconstruction_impact(...)
```

Notebook 不要。

Notebook 应展开：

```python
representation = ...
spatial_context = ...
local_state = ...
state_regions = ...
programs = ...
integration = ...
```

这几个变量本身就是主逻辑。

---

# 三十、但是 Notebook 与 batch 必须共享这些阶段函数

这是当前应该做的一次关键重构。

当前很多：

```python
_write_xxx()
```

同时：

```text
计算
+
写文件
```

不利于 Notebook 共用。

---

## 30.1 建议逐步改成

```python
compute_xxx(...)
    ↓
DataFrame / dict

save_xxx(...)
```

而不是：

```python
_write_xxx(...)
```

一口气全做。

---

## 30.2 Reconstruction Impact 最少形成这些稳定 stage functions

不需要 class。

简单函数即可：

```python
compute_representation(...)
build_spatial_context(...)
compute_local_state(...)
compute_state_regions(...)
compute_gene_spatial(...)
compute_program_activity(...)
compute_spatial_integration(...)
```

---

## 30.3 Batch

调用：

```text
compute
↓
save
```

---

## 30.4 Notebook

调用相同：

```text
compute
↓
display
↓
plot
```

这样真正解决：

> Notebook/batch 两套逻辑漂移。

---

# 三十一、现有 methods/regions.py 可以继续用，但需要补几个小能力

目前已有：

```text
assign_square_windows
compute_window_diversity
select_window_scale
select_region_threshold
assign_anatomy_candidates
flag_region_windows
```

很好，不要推翻。

建议增加：

```text
summarize_diversity_by_anatomy
summarize_region_by_anatomy
aggregate_feature_to_windows
compute_subtype_window_enrichment
```

这些都是明确会反复复用的能力。

不要再加复杂 Region class。

DataFrame 就够。

---

# 三十二、plotting/impact.py 需要重点补齐

目前主要有：

```text
partition sizes
gain histogram
anatomy
membership
window field
```

下一轮至少补：

```text
plot_window_scale_diagnostic
plot_local_diversity_comparison
plot_local_diversity_by_anatomy
plot_threshold_diagnostic
plot_region_anatomy_overlay
plot_region_extent_by_anatomy
plot_subtype_enrichment
plot_moran_raw_svc
plot_program_spatial
plot_program_by_anatomy
plot_program_region_overlap
```

这些比增加新 metric 更重要。

因为当前瓶颈是：

> **算出来了但看不出逻辑。**

---

# 三十三、输出目录也应该顺着科学主线整理

当前：

```text
tables/
```

越来越平。

建议只增加一层，不做复杂 artifact system：

```text
reconstruction_impact/
├── result.json
│
├── tables/
│   ├── representation/
│   ├── spatial_context/
│   ├── local_state/
│   ├── regions/
│   ├── programs/
│   └── integration/
│
├── figures/
│   ├── representation/
│   ├── spatial_context/
│   ├── local_state/
│   ├── regions/
│   ├── programs/
│   └── integration/
│
└── report.html
```

这样：

```text
Notebook
Web
Agent
```

都更容易理解。

---

# 三十四、result.json 不需要复杂 schema，但 outputs 要更有逻辑

可以仍然保持简单。

例如：

```json
{
  "analysis": "reconstruction_impact",
  "status": "succeeded",
  "sample_id": "P1CRC",
  "parameters": {...},
  "outputs": {
    "representation": "tables/representation/",
    "spatial_context": "tables/spatial_context/",
    "local_state": "tables/local_state/",
    "regions": "tables/regions/",
    "programs": "tables/programs/",
    "integration": "tables/integration/",
    "report": "report.html"
  }
}
```

不需要逐个 artifact ontology。

---

# 三十五、Web 必须从“文件浏览器”改成真正的科学阅读页面

Web 应完全对应科学问题。

---

## Web Overview

建议五个入口：

```text
1. Representation
2. Local state
3. Spatial localization
4. Gene / Program
5. Integrated spatial evidence
```

而不是：

```text
根据文件名 substring 自动分类
```

---

# 三十六、每个 Web section 用固定阅读模式

例如：

## Spatial localization

```text
Question
↓
Spatial definition / scale
↓
Diagnostic
↓
Main map
↓
Scope comparison
↓
Anatomy stratification
↓
Detailed tables
```

---

# 三十七、Web 不需要自动解释科学结论

可以自动给：

```text
Raw n
SVC n
region fraction
median score
```

这种 factual summary。

但不要自动写：

```text
reconstruction restored biology
```

这种 interpretation。

---

# 三十八、Standalone Moran / Pathway analyses 不应默认和 Impact 重复跑

当前 project template：

```text
reconstruction_impact
spatial_autocorrelation
pathway_activity
```

都启用。

但 Impact 内部已经计算：

```text
Moran
AUCell
```

会重复。

---

## 建议

默认 project：

```text
reconstruction_impact
```

即可。

Standalone：

```text
spatial_autocorrelation
pathway_activity
```

保留为：

> 用户单独想跑这项分析时的入口。

不要默认三次计算同一个东西。

---

# 三十九、第二条后续主线：Recoverability

当 Reconstruction Impact 做完整后，下一条大主线应该不是 CCI。

而是：

> **Recoverability across datasets/platforms**

因为这是后续大规模 10+ technology 分析真正需要的标准体系。

建议形成：

```text
02_recoverability.ipynb
```

但它仍然是一条大主线，不按 metric 拆 Notebook。

---

# 四十、Recoverability 建议分成三个层级

## 1. Cellular / State Recoverability

回答：

> reconstruction 后细胞状态是否更可分、更一致、更有生物学可解释性？

可包括：

```text
partition complexity
ARI / AMI（有 reference 时）
TMP / MER
ASW
subtype recoverability
marker specificity
```

---

## 2. Program Recoverability

回答：

> biological programs 是否更加可恢复？

例如：

```text
pathway activity
EMT
immune programs
CAF programs
program spatial structure
```

---

## 3. Spatial Recoverability

回答：

> spatial organization 是否更加可恢复？

例如：

```text
Moran
local spatial coherence
Region continuity
anatomy localization
state localization
```

---

# 四十一、Recoverability 与 Reconstruction Impact 的关系

不要复制代码。

应该：

```text
Reconstruction Impact
    ↓
产生单样本分析 evidence

Recoverability
    ↓
复用其中的标准指标
    ↓
跨 sample / platform 汇总
```

也就是说：

> Impact 更像 within-sample “发生了什么”。

> Recoverability 更像 cross-dataset “这种变化是否系统性存在”。

---

# 四十二、第三条未来主线：Spatial Biological Follow-up

等前两条稳定以后，再做：

```text
03_spatial_biology.ipynb
```

这条才负责：

```text
Treg
Mreg
CAF states
CCI
TLS
Niche
Trajectory
```

---

# 四十三、这条主线不是所有数据默认全跑

它应该来自前面的发现。

例如：

```text
T State Region
+
Treg enriched
+
Mreg enriched
+
Interface localization
```

然后才进入：

```text
Treg ↔ Mreg CCI
```

再看：

```text
TLS / niche
```

---

# 四十四、未来 Agent 应建立在这个层次上

Agent 最终应该：

```text
read Reconstruction Impact
↓
read Recoverability
↓
发现值得追的问题
↓
调用 Spatial Biological Follow-up methods
```

但当前只留接口。

不用现在实现 planner。

---

# 四十五、跨样本 Web 是 Recoverability 阶段再做

当数据多起来以后增加：

```text
output/index.html
```

可以按：

```text
Sample
Platform
Tissue
Reconstruction Impact
Recoverability
```

浏览。

现在先把单 sample 页面做好。

---

# 四十六、真实 P1 必须成为下一轮第一验收样本

不要继续以 synthetic fixture 为主要科学验收。

需要：

```text
复制真实 Raw
复制真实 SVC
写 sample.yaml
```

加入：

```text
data/P1CRC_VisiumHD/
```

H5AD 继续 gitignore。

---

# 四十七、真实 P1 验收必须回答下面这些问题

### Representation

旧结果还能否观察？

### Anatomy

Tumor / Interface / Normal 是否合理？

### Window scale

候选尺度中哪个合理？为什么？

### Local state

双 baseline 是否恢复？

### State Region

旧 high-diversity Region 是否基本能复现？

### Program

EMT spatial localization 是否能恢复？

### Integration

State / Interface / EMT 之间能否重新串起来？

---

# 四十八、验收不能只有 CSV parity

应有三层。

## Layer 1：method correctness

Unit tests。

---

## Layer 2：numeric parity / deliberate difference

对于未改变定义的结果：

```text
应基本 parity
```

对于明确修改的方法：

```text
记录 old/new difference
```

---

## Layer 3：visual/scientific logic parity

检查：

```text
以前能看见的空间现象
现在是不是仍能看见
```

这是当前最重要的一层。

---

# 四十九、建议专门建立一张 Migration Difference Table

不要再混在长文档里。

例如：

| Component              | Old              | New                   | Type                         | Action             |
| ---------------------- | ---------------- | --------------------- | ---------------------------- | ------------------ |
| Raw Level2 mapping     | Analysis-side OT | upstream annotation   | required architecture change | keep new           |
| paired iST carrier     | spatial+expr     | unified SVC           | required architecture change | keep new           |
| window scale           | support-guided   | fixed default         | unintended loss              | restore diagnostic |
| anatomy stratification | available        | partial               | unintended loss              | restore            |
| EMT spatial field      | available        | distribution mainly   | unintended loss              | restore            |
| sampling               | scope-specific   | side-first            | scientific change            | review/fix         |
| HVG                    | old problematic  | inherited problematic | old method issue             | independently fix  |

这个表以后非常有用。

---

# 五十、建议 Codex 下一轮按下面几个工作包推进

## Work Package 1 — Foundation correction

修改：

```text
io.py
config
parameter loading
sampling
exception handling
```

验收：

```text
P1 可以真实读取
Notebook/batch 参数一致
```

---

## Work Package 2 — Spatial context

补：

```text
window scale diagnostics
anatomy
common grid
```

并在 Notebook 看图。

---

## Work Package 3 — Local state

补：

```text
Raw Leiden
Raw Level2
SVC
dual baseline
anatomy stratification
```

---

## Work Package 4 — Regions

补：

```text
State
threshold diagnostics
State × anatomy
Gain audit
subtype enrichment
```

---

## Work Package 5 — Gene / Program localization

补：

```text
Moran by scope
shared-gene comparison
program spatial field
window aggregation
program × anatomy
program × State
```

---

## Work Package 6 — Integration

补：

```text
Region × Anatomy
Subtype × Region
Program × Region
```

和主 Notebook 第 9–10 节。

---

## Work Package 7 — Batch parity

让：

```text
Notebook stage functions
=
batch stage functions
```

正式 batch 重跑 P1。

---

## Work Package 8 — Web rewrite

不要文件名分类。

按照主科学问题显式组织。

---

# 五十一、这一阶段不应该做的东西

继续不要：

```text
Agent planner
trigger graph
complex finding schema
full provenance engine
universal pairing
legacy sc-SVC compatibility
huge report framework
many tiny notebooks
```

也不要因为想到以后 CCI/TLS 就提前把所有 biological analysis 塞进 Reconstruction Impact。

---

# 五十二、最终希望仓库达到的状态

完成这一阶段后：

```text
revise_analysis/methods/
    └─ 足够模块化、可复用

revise_analysis/analyses/reconstruction_impact.py
    └─ 稳定组合这些 methods 供 batch 使用

notebooks/01_reconstruction_impact.ipynb
    └─ 把整条科学主线展开给人看
       可局部调参数
       可看关键中间图
       可判断实现是否合理

output/<sample>/reconstruction_impact/
    └─ 稳定 tables + figures + result.json

report.html
    └─ 按科学问题连续阅读
```

之后再自然扩展：

```text
02_recoverability.ipynb
03_spatial_biology.ipynb
```

而不是继续扩展软件框架本身。

---

# 五十三、给 Codex 的最高优先级判断规则

遇到一个新需求时先问：

### 这是 reusable computation 吗？

是：

```text
→ methods/
```

### 这是稳定的 Reconstruction Impact 组合逻辑吗？

是：

```text
→ analyses/reconstruction_impact.py
```

### 这是正在探索的逻辑关系或参数吗？

是：

```text
→ 01_reconstruction_impact.ipynb 对应位置
```

### 这是已经稳定的结果阅读组织吗？

是：

```text
→ report.html renderer
```

### 这是未来 Agent 的智能决策吗？

是：

```text
→ TODO
```

当前不要实现。

---

# 五十四、这一轮真正的开发原则

> **不要继续把“功能存在”当成完成。一个 Reconstruction Impact 能力只有在下面四件事都成立时才算真正完成：**

```text
① module 能算
② Notebook 能沿主逻辑看到并检查
③ batch 能用同一逻辑稳定跑
④ Web 能沿科学问题读懂结果
```

如果只有：

```text
CSV 生成成功
```

还不能算完成。

---

# 五十五、下一轮的最高优先级路径

如果只能按顺序做，我建议严格：

```text
真实 P1 输入
    ↓
参数 / sampling / partition 基础修正
    ↓
window scale + anatomy
    ↓
双 baseline local state
    ↓
State Region + anatomy
    ↓
EMT/program spatial localization
    ↓
Integrated spatial localization
    ↓
Notebook/batch parity
    ↓
Web scientific reading
```

**先把这一条完整跑通，再扩展 Recoverability。**

这是当前最能把仓库从“已经搭了不少代码”推进到“真正实现我们最初设计目标”的路径。

其中我认为最重要的一个调整是：**下一轮不要以“还缺哪些分析函数”为中心开发，而要以 `01_reconstruction_impact.ipynb` 这条科学主线缺哪一环为中心开发。** 每补一环，底层有可复用部分就沉到 `methods/`，稳定后再接入 batch 和网页。这样开发方向会明显更接近你最初想要的结果可视化驱动方式。
