# GvpA 中期汇报 v2：讲稿、要求核对与问答

汇报人：史晨皓。9 页主讲 + 3 页隐藏问答备用页；建议主讲约 6 分钟。

## 按教师要求核对

依据 `intro-mlproj.pdf` 第 8、10 页及用户提供截图，中期检查关注模型定制设计、Baseline 实现、初步实验验证，以及 Baseline 代码库和实验记录。课程文件同时注明不同项目可有不同内部路径。

| 要求 | PPT 对应 | 已有证据与限制 |
|---|---|---|
| 模型定制设计 | 第 4 页 | 上下文残基表征、滑动窗口、PCA/HDBSCAN、双坐标；已实现 |
| Baseline 实现 | 第 5 页 | ESM-2 与氨基酸组成基线均已运行；没有功能准确率 |
| 初步实验验证 | 第 6–7 页 | 12 组参数、50 次随机诊断与长度混杂；独立 validation 尚未运行 |
| 代码库 + 实验记录 | 第 9 页、附录 A | 配置、代码、台账、真实运行 SHA 与哈希；共享访问与来源仍需闭环 |

**判断：新版覆盖中期汇报要求的各个内容项，但 PPT 不能代替未完成的独立验证、团队审查或 M3 最终验收。不能承诺教师已经认可无监督路线；需在问答中如实说明路线调整依据。**

## 约 6 分钟逐页讲稿

### 第 1 页｜开场（约 20 秒）

各位老师好，我是史晨皓。我们研究 GvpA 的潜在功能区域识别。这次按中期要求，汇报模型设计、Baseline 实现和初步实验，以及代码与实验记录的交付情况。

### 第 2 页｜任务定义（约 40 秒）

最初我们考虑功能分类，但所有可扫描序列都命中同一家族，缺少可靠的功能差异标签。因此改为区域识别：第一，表征能否发现可定位模式；第二，它是否只是长度或位置效应；第三，能否得到独立生物学证据。输出是有证据等级的候选区域，不把聚类直接称为功能证明。

### 第 3 页｜数据与 EDA（约 40 秒）

M2 已验收。2078 条去重候选按质量分层，其中 1721 条进入主分析；partial 和类型冲突保留作敏感性队列，另有两条被质量排除。2076 条可扫描序列按 478 个同源簇划分，没有跨集合簇泄漏。M3 只在 1202 条主分析 discovery 序列上拟合；519 条主分析 validation 继续隔离。

### 第 4 页｜模型定制设计（约 45 秒）

方法先把整条序列送入 ESM-2 35M，取得第十二层、每个残基 480 维的上下文表示。再对三十个残基的窗口做平均池化，而不是重新把短片段当成独立序列编码。主配置使用 PCA 二十维和 HDBSCAN。窗口保留原序列及 HMM 坐标，排除特殊 token；UMAP 只作展示，全部拟合在 discovery 内完成。

### 第 5 页｜Baseline 对照（约 40 秒）

我们实现了二十维氨基酸组成基线，采用相同窗口及聚类流程。ESM 主配置四簇，噪声比例百分之三十点六八；组成基线七十四簇，噪声百分之六十八点三三。二者 ARI 为负，说明划分差异较大。不过簇数和噪声比例没有功能真值含义，不能凭这些数字宣称 ESM 更准确。

### 第 6 页｜探索性候选（约 45 秒）

初步分析中，簇三包含 2216 个窗口和 265 条序列，多数覆盖 HMM 十四到三十九区间，其中百分之六十三点九四覆盖状态二十五。按长度及相对位置分箱匹配的五十次随机诊断，九十五分位为百分之五十九点九一。这提供了一点探索性支持，但窗口重叠，而且区域是在发现后选择的，不能当作功能显著性结论。

### 第 7 页｜稳健性与混杂（约 55 秒）

稳健性检查显示，修改最小簇大小影响较小，但改变 PCA 维度或窗口尺度后结果变化明显。更关键的是，短序列中有 245 条进入簇三，而一百三十一到一百六十残基的 845 条序列只有一条进入。只用长度区分参与情况，AUC 达到零点九六八八。这是混杂诊断，不是功能预测成绩，所以候选尚未冻结，validation 尚未运行。

### 第 8 页｜成员分工（约 30 秒）

成员分工是：我负责统筹、冻结规则和证据整合；韩涛鸿负责数据来源和划分；燕相楠负责 HMM 与结构坐标；李若凡负责表征、聚类及后续验证；杨傲翔负责匹配基线、扰动归因和图表。后续各模块使用同一个冻结候选与坐标版本，避免各做各的、结果无法对应。

### 第 9 页｜中期交付与下一步（约 40 秒）

对照课程要求，模型定制设计和 Baseline 已实现，也完成了参数、随机和混杂等初步实验。代码、配置、台账及带真实 SHA 和哈希的收据已有交付。当前满足这些内容的展示，但尚未完成独立验证和 M3 验收。下一步先控制长度混杂，再冻结候选，执行锁定验证，并补充生物学与扰动证据。

## 2–4 分钟问答

**为什么不用监督分类？** 天然数据没有可靠的功能差异标签。全部可扫描序列的 PF00741 命中只提供家族身份；我们把目标调整为提出可检验区域假设。

**这次是否完成 M3？** 已有模型与 Baseline、发现集初步实验和可追溯运行记录。候选未冻结、独立 validation 未运行、正式扰动未完成，因此 M3 未验收。

**ESM-2 比组成基线更好吗？** 当前只能说聚类结果不同，ESM 噪声比例更低。没有功能真值，不能用簇数、噪声或 ARI 证明功能识别更好。

**AUC 0.9688 是什么？** 仅用序列长度区分是否参与簇 3 的诊断 AUC，反映强混杂。它不是功能预测准确率。

**63.94% 超过随机 59.91% 是否显著？** 这是 50 次分箱匹配随机诊断的探索性比较。窗口重叠、区域后选择和残余混杂尚在，不能给出功能显著性或已验证结论。

**HMM 区间能直接和论文突变位置比较吗？** 不可以。HMM 状态、原序列位置、结构残基编号必须通过映射表对应；还要核对突变亲本。见附录 C。

**为什么 validation 未运行？** 先在 discovery 修订和冻结候选与评价规则，再一次性固定模型迁移。利用 validation 调候选会破坏独立验证。

**代码与实验记录在哪？** 仓库 `features/esm2_embed.py`、`scripts/m3_discovery.py`、`scripts/m3_candidate_robustness.py`，配置 `configs/m3_discovery.json`，台账 `experiments/registry.csv`，收据 `reports/M3/runs/<运行号>/run_receipt.json`。附录 A 列出运行号。

## 附录页的使用

第 10–12 页默认隐藏：A 为可复现运行证据，B 为成员交接与证据缺口，C 为坐标和验证解释边界。只在问答需要时跳转，不计入 6 分钟主讲。

## 实际运行 SHA 与演示命令

- EMBED / DISCOVERY-001：`dc066b76c637cdb1b90b7c8490e28ad81c4cff52`。
- COHORT-001：`eb2292e41e638d8efd33c8bc84bbfb71e9368d22`。
- ROBUST-001：`8a0bf23715783de9fc1c5fe9fd3e6a5bb83644ee`。

原始命令含确切输入路径与环境，见各收据 `provenance.command`。正式演示优先打开已有真实结果与收据；临时换环境不保证可即时下载模型或重跑。

## 来源与版式检查说明

所有页面的来源也写入 PPT 备注。下列路径均为仓库实际证据；成员姓名取自用户确认。旧报告文字与原始 CSV 有差异时，本版优先引用原始 CSV：组成基线簇数为 **74**。

- 第 1 页：intro-mlproj.pdf，第8、10页；用户提供的中期要求截图
- 第 2 页：MLP/results/bioinformatics/pf00741_scan_summary.json；MLP/reports/M3/candidate_freeze_decision.md
- 第 3 页：MLP/results/data_audit/gvpa_v1_audit_summary.json；MLP/results/data_audit/gvpa_v1_split_summary.json
- 第 4 页：MLP/features/esm2_embed.py；MLP/scripts/m3_discovery.py:window_rows / cluster_fit；MLP/configs/m3_discovery.json
- 第 5 页：MLP/reports/M3/runs/M3-D-DISCOVERY-001/simple_baseline.csv；MLP/reports/M3/runs/M3-D-DISCOVERY-001/parameter_sensitivity.csv
- 第 6 页：MLP/reports/M3/runs/M3-D-DISCOVERY-001/candidate_regions.csv；MLP/reports/M3/runs/M3-D-DISCOVERY-001/window_umap.png
- 第 7 页：MLP/reports/M3/runs/M3-D-DISCOVERY-001/parameter_sensitivity.csv；MLP/reports/M3/runs/M3-D-ROBUST-001/summary.json；MLP/reports/M3/runs/M3-D-ROBUST-001/length_bins.csv
- 第 8 页：本次用户确认的姓名与角色；MLP/docs/plans/2026-10-07-m3-midterm-sprint.md
- 第 9 页：intro-mlproj.pdf，第8、10页；MLP/experiments/registry.csv；MLP/reports/M3/candidate_freeze_decision.md
- 第 10 页：MLP/experiments/registry.csv；MLP/reports/M3/runs/M3-D-DISCOVERY-001/run_receipt.json；MLP/reports/M3/runs/M3-D-ROBUST-001/run_receipt.json
- 第 11 页：MLP/reports/M3/c_handoff_technical_check.md；MLP/reports/M3/e_handoff_technical_check.md；MLP/docs/plans/2026-10-07-m3-midterm-sprint.md
- 第 12 页：MLP/reports/M3/candidate_freeze_decision.md；MLP/reports/M3/c_handoff_technical_check.md

本环境缺少 artifact-tool 运行包，且 Office COM 无法创建会话；已通过 OOXML/ZIP、关系目标、文本几何和近似预览检查，未完成原生 PowerPoint 渲染。请在汇报电脑上放映确认字体和换行。
