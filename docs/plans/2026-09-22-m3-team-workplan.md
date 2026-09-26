# M3 团队分工与推进计划

- 状态：待启动；M2未正式关闭，仅允许代码与计划草稿准备，不执行正式M3实验
- 负责人：A（项目管理与证据统筹）
- 适用阶段：M3——ESM-2 表示与候选功能区域发现
- 制定日期：2026-09-22
- 跟踪任务：[Issue #6](https://github.com/Lithium9767/MLP/issues/6)

## 1. 阶段目标

在不使用 validation 集参与发现或调参的前提下，从 discovery 集生成可追溯的 ESM-2 全长及滑动窗口表示，完成候选区域发现、混杂检查和参数敏感性分析，并为 M4 的 HMM/MSA、结构、实验突变和残基扰动验证准备统一坐标。

本阶段不训练天然 GvpA 功能分类器，也不把 PF00741 命中、聚类结果或 attention 单独解释为功能证明。

## 2. 开始条件

2026-09-26核实main为 `5891c2f7a6eb79e48dd7b75175a3af8176e410a2`。PR #13/#14已归档M2技术结果，但图表错误、运行来源、共享原始输出和真实独立复核尚未解决，见[Issue #15](https://github.com/Lithium9767/MLP/issues/15)及[收尾PR #16](https://github.com/Lithium9767/MLP/pull/16)。本计划不能作为M2已关闭的证据。

1. M2独立复核报告进入main、最终SHA测试通过、README/台账/验收报告准确、PR #9经真实非作者审查合并、原始输出位置和哈希可追溯，并由A本人在Issue #15留下最终验收评论。
2. 不整体合并 PR #4。仅准备ESM-2入口、最小依赖、ESM相关测试和迁移审计的替代草稿；M2关闭后再完成非作者审查并合并。HMM与结构使用main已归档实现，不用旧代码覆盖。
3. 所有成员从合并后的最新 `main` 新建任务分支；复用已有Issue #5/#6，相同任务不重复创建。
4. 使用已冻结的数据版本 `gvpa-recognition-c7f6f005d717`。
5. 使用已冻结的划分版本 `homology-8b9005e2d9-s42`。

## 3. 共同约束

- discovery 用于发现候选区域和选择参数；候选区域冻结前不得查看 validation 结果。
- primary 队列用于主分析；partial 和类型冲突记录只进入单独的敏感性分析。
- PF00741 用于家族确认和同源坐标，不作为功能正负标签。
- UMAP 只用于可视化，聚类应在原始表示或保留足够信息的 PCA 空间中进行。
- attention 只作为 sequence-dependency evidence，不能单独作为关键位点归因。
- 原始数据、模型权重、PDB 文件、大型 embedding 和缓存不提交普通 Git。
- 每次正式实验登记数据版本、split 版本、配置、随机种子、软件版本、完整 Git commit、输出哈希和限制。
- 所有图表和最终统计量由脚本生成，不手工修改最终数字。

## 4. 工作依赖

```text
A：冻结范围、验收规则、Issue/PR 与实验登记
                         |
                         v
B：冻结数据输入与 split 交接
             +-----------+-----------+
             |                       |
             v                       v
C：PF00741/HMM/MSA/结构坐标      D：ESM-2 表示与聚类
             |                       |
             +-----------+-----------+
                         |
                         v
                  候选区域冻结
                         |
                         v
             E：扰动、随机基线和可视化
                         |
                         v
                 A：M3 汇总与验收
```

## 5. 成员分工

### A：项目管理与证据统筹

职责：

1. 组织ESM-2最小迁移替代PR的非作者审查；保留PR #4及其历史，不整体合并。
2. 为 B、C、D、E 建立任务 Issue，明确输入、输出、依赖和验收条件。
3. 冻结 M3 的 discovery/validation 使用边界和候选区域冻结规则。
4. 建立实验编号，检查 `experiments/registry.csv` 登记是否完整。
5. 组织阶段检查、PR 审查和结果复现。
6. 汇总 M3 报告，区分发现结果、验证结果和未验证假设。
7. 向教师确认无监督候选功能区域路线是否满足课程目标。

交付物：

- M3 Issue 与验收清单；
- 实验编号和证据矩阵；
- `reports/M3/m3_summary.md`；
- 阶段汇报材料；
- PR 与结论审查记录。

验收标准：

- 每项结论能追溯到实验编号、配置、数据版本和 Git commit；
- 每个 PR 至少有一名非作者审查；
- validation 未在候选区域冻结前参与发现或调参；
- 报告不把计算关联表述为湿实验功能证明。

### B：数据工程与输入冻结

职责：

1. 独立复核 M2 数据数量、同源簇和 discovery/validation 隔离。
2. 为 C 和 D 生成使用同一 `internal_id` 的稳定输入 manifest。
3. 检查 accession、序列、sequence hash、analysis cohort 和 split 一致性。
4. 分别准备 primary、partial、type conflict 和组合敏感性分析清单。
5. 记录本地大型文件的位置、版本和 SHA-256，不上传受限数据。

建议分支：`feature/<issue-id>-m3-data-handoff`

交付物：

- M3 输入 manifest；
- 数据完整性检查脚本；
- 数据与 split receipt；
- discovery/validation 和各 cohort 数量汇总；
- 输入文件哈希清单。

验收标准：

- 共 2076 条 sequence-QC-eligible 序列进入冻结划分；
- discovery 为 1453 条、335 个同源簇；
- validation 为 623 条、143 个同源簇；
- primary discovery 为 1202 条，primary validation 为 519 条；
- 同源簇跨集合泄漏为 0；
- 不重新随机生成或覆盖已冻结 split。

审查：C 审查字段和坐标接口，A 审查数据版本和实验边界。

### C：PF00741、MSA 与结构坐标

职责：

1. 审查远端 C 分支中的 A3M、MSA、Pfam、区域映射和 workflow 模块。
2. 从最新 `main` 选择性迁移可复用模块，不整体合并旧 C 分支。
3. 交接M2已有PF00741扫描和坐标，补齐原始文件访问；仅在复核发现问题时按新实验编号重跑，不无故重复已归档工作。
4. 建立天然序列位置、窗口位置、PF00741 HMM match state、文献位置和 PDB 7R1C 编号之间的显式映射。
5. 计算 HMM 列 occupancy、consensus、identity、conservation 及 insertion/deletion 情况。
6. 完成 7R1C 结构和二级结构映射。

建议分支：`analysis/<issue-id>-pf00741-coordinate-system`

交付物：

- `pf00741_run_receipt.json`；
- `hmm_coordinate_map.csv`；
- `hmm_conservation.csv`；
- `structure_7r1c_map.csv`；
- 坐标映射测试和运行说明。

验收标准：

- Pfam/HMM、软件版本和输入哈希完整；
- 插入、缺失和 partial 不会造成静默错位；
- 任一天然序列窗口可以显式映射到 HMM 列；
- PF00741 命中不被解释为功能标签；
- PDB 原始文件不进入普通 Git。

审查：B 审查输入 ID 和序列，E 审查残基证据接口。

### D：ESM-2 表示与候选区域发现

职责：

1. 对合并后的 ESM-2 入口执行小样本冒烟测试。
2. 为 discovery primary 的 1202 条序列生成全长 mean pooling 表示。
3. 以 30 aa 窗口、步长 5 生成局部残基窗口表示。
4. 检查维度、token 对齐、NaN、缓存和重复运行一致性。
5. 实现 PCA 和 HDBSCAN；UMAP 仅用于二维展示。
6. 检查长度、物种、partial、类型冲突和 C 端延伸等混杂因素。
7. 完成窗口、PCA 和 HDBSCAN 参数敏感性分析。
8. 生成候选窗口表并映射到原始位置和 HMM 坐标。
9. 在 A 批准后冻结候选区域，再对 validation 进行独立验证。

建议分支：`model/<issue-id>-esm2-region-discovery`

交付物：

- ESM-2 运行 receipt 和 embedding manifest；
- PCA 解释方差表；
- HDBSCAN 聚类指标；
- 参数敏感性和混杂检查表；
- 候选窗口表与候选区域冻结文件；
- PCA/UMAP 图及生成脚本；
- validation 独立验证摘要。

验收标准：

- 发现阶段只读取 discovery；
- 模型名称、权重版本、层、池化、设备和软件版本完整；
- 大型 embedding 不提交 Git；
- 候选区域能映射到原始序列位置和 HMM match state；
- 至少与随机窗口或简单序列统计基线比较；
- validation 未参与超参数选择。

审查：E 审查窗口和残基接口，A 审查是否回答 RQ1。

### E：残基扰动、随机基线与可视化

职责：

1. 在等待 D 冻结候选区域期间，用合成序列实现通用扰动接口。
2. 至少实现 residue masking 或 alanine substitution 中的一种方法。
3. 建立相同长度、相似位置分布的随机区域基线。
4. 定义原始表示与扰动表示之间的变化分数。
5. 接收 D 的候选区域和 C 的 HMM/结构坐标，生成统一证据表。
6. 制作序列轨迹、HMM 保守性、扰动分数、突变证据和结构位置图。

建议分支：`attribution/<issue-id>-residue-perturbation`

交付物：

- 残基扰动脚本；
- 随机基线脚本；
- token/序列/HMM/PDB 坐标测试；
- 候选区域和残基证据表；
- 可视化脚本与图表说明。

验收标准：

- 扰动前后 token 和残基坐标没有偏移；
- 扰动分数和聚合规则明确；
- 候选区域与随机区域使用相同统计方法；
- attention 只作为补充列；
- 图表由脚本生成并绑定实验编号；
- 结果可以映射到 HMM 列和 PDB 编号。

审查：D 审查模型计算，C 审查坐标和结构解释。

## 6. 一周推进节奏

以下“第1天”从M2正式关闭、M3启动条件满足后开始计时；此前不启动embedding或聚类实验。C沿用M2结果，正式重跑仅用于解决复核差异。

| 时间 | A | B | C | D | E |
| --- | --- | --- | --- | --- | --- |
| 第 1 天 | 审查ESM-2替代PR，复用已有Issue | 复核数据版本 | 交接main已有坐标实现 | 环境和小样本冒烟测试 | 建立扰动接口 |
| 第 2 天 | 冻结 M3 验收标准 | 生成输入 manifest | 核对M2坐标与receipt | discovery 全长 embedding | 设计随机基线 |
| 第 3 天 | 检查实验登记 | 完成 C/D 数据交接 | HMM/MSA 坐标 | 滑窗 embedding | 完成合成测试 |
| 第 4 天 | 组织中期检查 | 支持数据问题复核 | 保守性和 7R1C 映射 | PCA/HDBSCAN | 建立可视化模板 |
| 第 5 天 | 审查候选区域 | 复核结果 ID | 审查区域坐标 | 参数敏感性和混杂检查 | 接入候选区域 |
| 第 6 天 | 批准候选区域冻结 | 支持独立复跑 | 输出生物信息学证据 | validation 验证 | 扰动与随机比较 |
| 第 7 天 | 汇总报告和组织验收 | 提交说明与 PR | 提交结果与 PR | 提交结果与 PR | 提交结果与 PR |

## 7. M3 最终验收包

M3 结束时至少应具备：

1. 冻结的数据和 split 版本；
2. ESM-2 全长与窗口表示 receipt；
3. PCA/HDBSCAN 候选区域结果；
4. 参数敏感性和混杂因素检查；
5. PF00741 HMM 统一坐标；
6. discovery 中冻结的候选区域；
7. validation 上的独立复现结果；
8. 至少一种残基扰动方法和随机基线；
9. 完整实验台账；
10. M3 阶段报告、图表和限制说明。

## 8. 完成判定

只有在候选区域已于 discovery 阶段冻结、validation 结果由固定脚本独立生成、所有结果绑定版本与提交、关键结论通过非作者复核后，M3 才可标记完成。候选区域和候选残基均使用“潜在”“候选”表述，直至获得独立实验验证。
