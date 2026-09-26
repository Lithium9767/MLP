# M2 小组协作与收尾计划

- 状态：技术结果已归档，收尾与真实独立复核待完成；本计划PR仍需非作者审查
- 负责人：A（项目管理与证据统筹）
- 适用阶段：M2——数据基础、同源划分与 PF00741/HMM 坐标冻结
- 制定日期：2026-09-23
- 跟踪任务：[Issue #8](https://github.com/Lithium9767/MLP/issues/8)

## 1. 阶段目标

本计划用于关闭 GvpA 项目的 M2。2026-09-26 核实 main 为 `5891c2f7a6eb79e48dd7b75175a3af8176e410a2`：PR #13/#14 已归档数据复现、PF00741真实扫描、HMM/结构坐标、敏感性和图表。尚缺真实独立复核、准确运行来源、共享原始输出访问与最终验收。全组应先完成 M2 验收，再启动正式 M3 分析。

### 当前半天至一天收尾安排

以下安排取代第6节原始五天日程；原任务清单保留作为交付验收标准，不要求重复已归档工作。

| 顺序 | 负责人 | 当前交付与门槛 |
| --- | --- | --- |
| 0–1小时 | A | 核实角色账号与贡献范围，在[Issue #15](https://github.com/Lithium9767/MLP/issues/15)指派未编写相关改动的B/C；不按PR账号猜测身份 |
| 1–3小时 | C、E | C提供原始输出固定共享位置、哈希与dirty运行来源；E修复cohort图四个0值（应20/305/30/2），从干净提交重跑并登记新编号；旧E运行SHA未知不得冒填归档SHA |
| 2–4小时 | D、非作者审查人 | 复核[收尾PR #16](https://github.com/Lithium9767/MLP/pull/16)，执行完整测试与validator；审查本计划PR #9 |
| 取得证据后 | 独立B/C | 实际下载哈希核对、源数据/图表抽查、泄漏检查，填写reports/M2/independent_review.md并提交PR；缺失项写未验证 |
| 最后 | A | 报告进入main、#9经非作者Review合并、最终SHA测试通过、所有阻塞解决后，亲自在Issue #15评论验收 |

半天至一天是工程准备目标，不是自动关闭期限。助手预检不等于B/C独立复核，合并后复核不冒充历史Review。PR #11已关闭未合并，其C实现与结构摘要由#13整合，不重复合并。真实报告、共享访问、图表修复或A验收未完成时保持待验收。

M2 最终需要交付：

1. 可复现的数据审计；
2. 可复现的 MMseqs2 同源聚类；
3. 冻结的 discovery/validation 划分；
4. 真实 PF00741 全量扫描；
5. 天然序列到 PF00741 HMM match state 的坐标映射；
6. PDB 7R1C 结构坐标映射；
7. 输入、软件、参数和结果哈希；
8. 非作者复核记录；
9. 完整实验台账；
10. M2 最终报告和限制说明。

## 2. 当前冻结基线

| 项目 | 当前结果 |
| --- | ---: |
| 候选代表序列 | 2078 |
| 原始记录 | 17108 |
| 可进入同源聚类 | 2076 |
| Primary 队列 | 1721 |
| MMseqs2 同源簇 | 478 |
| Discovery | 1453 条、335 簇 |
| Validation | 623 条、143 簇 |
| 同源簇泄漏 | 0 |
| 数据版本 | `gvpa-recognition-c7f6f005d717` |
| Split 版本 | `homology-8b9005e2d9-s42` |

以上数字是复现时需要核对的基线。任何差异必须先在 Issue 中说明原因，不得直接覆盖冻结结果。

## 3. 共同规则

- 每项工作先建立 Issue，从最新 `main` 创建独立分支，通过 PR 合并。
- 每个 PR 至少由一名非作者审查。
- 不整体合并旧 B/C 分支，只从最新 `main` 选择性迁移经过审查的代码。
- 原始 `gv.zip`、逐条序列表、HMM 原始输出、PDB 文件和大型中间文件不进入公共 Git。
- 正式运行必须记录数据版本、split 版本、配置、随机种子、软件版本、Git commit、输出哈希和限制。
- discovery/validation 划分保持冻结；M2 期间不重新随机划分。
- PF00741 只用于家族身份、同源坐标和质量核对，不作为天然 GvpA 的功能标签。
- HMM 可以映射全部序列，但 validation 的保守性不能提前用于 M3 候选区域发现。
- 所有最终数字和图表由脚本生成，不手工修改。

## 4. 工作依赖

```text
A：建立 Issue、冻结验收规则、组织 PR 审查
                         |
                         v
B：独立复现数据审计与同源划分
                         |
                         v
C：PF00741 全量扫描、HMM 与 7R1C 坐标
             +-----------+-----------+
             |                       |
             v                       v
D：自动验收与敏感性分析       E：抽查、图表与解释复核
             |                       |
             +-----------+-----------+
                         |
                         v
A：更新台账、冻结报告、关闭 M2
```

B 和 C 可以准备环境和代码，但 C 的正式运行必须使用 B 确认后的输入 manifest。D 和 E 可以提前准备脚本与模板，最终结果必须等待 B/C 的冻结输出。

## 5. 成员任务

### A：项目统筹与最终验收

职责：

1. 复用 Issue #8/#10/#12 及合并后独立复核 Issue #15，不重复建立相同任务。
2. 为每项任务填写负责人、输入、输出、依赖、审查人和验收标准。
3. 保留 PR #4 及其审查记录；不整体合并。M2继续使用main已有HMM/结构实现，ESM-2另做M3最小迁移草稿。
4. 冻结本计划中的数据版本、split 版本和验收数字。
5. 维护任务状态，处理成员之间的交接和阻塞。
6. 检查每次正式运行是否登记到 `experiments/registry.csv`。
7. 汇总结果、方法、限制和复现命令，形成 M2 最终报告。
8. 确认公共仓库中不存在原始数据、论文、PDB、模型权重或凭据。

交付物：

```text
reports/M2/m2_acceptance.md
results/data_audit/m2_release_summary.json
experiments/registry.csv
M2 汇报材料
```

验收标准：

- 每个关键数字都能追溯到脚本、输入和结果文件；
- 每个正式结果绑定完整 Git commit；
- 每个 PR 完成非作者审查；
- M2 报告明确区分已完成、未完成和限制；
- A、B、C 共同确认数据定义后才关闭 M2。

### B：数据审计与同源划分独立复现

输入：

- 老师提供的 `gv.zip`；
- `preprocessing/prepare_gvpa_dataset.py`；
- `preprocessing/homology_split.py`；
- 当前数据审计摘要和 MMseqs2 参数。

职责：

1. 计算并核对 `gv.zip` 和 recognition JSON 的 SHA-256。
2. 从原始 ZIP 重新生成 metadata、member 表和聚类 FASTA。
3. 核对记录数、唯一 ID、序列哈希、partial、类型冲突、异常残基和 primary 队列。
4. 使用 MMseqs2 15-6f452 和冻结参数重新聚类：
   - `--min-seq-id 0.8`
   - `-c 0.8`
   - `--cov-mode 1`
5. 重新生成 discovery/validation 划分。
6. 检查未知 ID、缺失 ID、重复归属和同簇跨集合。
7. 比较新旧结果和输出哈希。
8. 生成独立复现 receipt；如有差异，在 Issue 中说明原因。
9. 协助 A 修正实验台账中仍显示 `split_pending` 的过时状态。

建议分支：`feature/<issue-id>-m2-data-reproduction`

GitHub 交付物：

```text
results/data_audit/m2_data_reproduction.json
results/data_audit/m2_split_reproduction.json
reports/M2/data_reproduction_review.md
```

大型交付物放共享存储：

```text
metadata.csv
members.csv
sequences_for_clustering.fasta
MMseqs2 原始输出
split_manifest.csv
```

验收标准：

- 候选序列为 2078；
- QC 合格序列为 2076；
- primary 为 1721；
- 同源簇为 478；
- discovery 为 1453 条、335 簇；
- validation 为 623 条、143 簇；
- cluster leakage 为 0；
- 新旧结果一致，或所有差异有证据和解释。

审查：C 审查字段和序列输入，A 审查数据版本和结果边界。

### C：PF00741、HMM 与结构坐标

输入：

- B 确认后的 2076 条序列和 manifest；
- PF00741 HMM；
- main 已整合的 `bioinformatics/m2_pf00741.py`、`m2_coordinates.py` 及下载模块；
- PDB 7R1C；
- 冻结的数据和 split 版本。

职责：

1. 审查 main 中 `download_hmms.py`、`download_7r1c.py`、`m2_coordinates.py` 和 `m2_pf00741.py`；不覆盖为PR #4旧实现。
2. 记录 PF00741 来源、数据库版本、下载日期和 SHA-256。
3. 安装并记录 HMMER/PyHMMER 及依赖版本。
4. 对 2076 条序列执行真实 PF00741 扫描。
5. 记录命中、未命中、失败、多命中、覆盖和阈值边界状态。
6. 建立天然序列残基到 PF00741 HMM match state 的映射。
7. 显式处理 insertion、deletion、partial 和异常记录。
8. 建立 PDB 7R1C 残基编号到 HMM 状态的映射。
9. 分别汇总 discovery 和 validation 的覆盖情况。
10. 保留完整原始输出，向 GitHub 只提交小型摘要、代码和哈希。

建议分支：`analysis/<issue-id>-pf00741-hmm-mapping`

GitHub 交付物：

```text
results/bioinformatics/pf00741_scan_summary.json
results/bioinformatics/hmm_coordinate_summary.json
results/bioinformatics/structure_7r1c_summary.json
reports/M2/pf00741_hmm_audit.md
```

共享存储交付物：

```text
PF00741.hmm
全量扫描原始输出
逐序列 HMM 坐标表
7R1C.pdb
逐残基结构映射表
```

验收标准：

- 输入数量、internal ID 和 sequence hash 与 B 的 manifest 一致；
- 软件、数据库、阈值、命令和输入哈希完整；
- 2076 条序列均有成功、失败或未命中的明确状态；
- 坐标约定明确，插入残基不会被错误分配 match state；
- PF00741 没有被解释为功能正负标签；
- PDB 原始文件没有进入普通 Git。

审查：B 审查输入和 ID，E 抽查残基与结构坐标。

### D：工程复现、自动验收与参数敏感性

输入：

- B 的数据和 split 结果；
- C 的 PF00741、HMM 和结构结果；
- 所有正式运行 receipt。

职责：

1. 建立 M2 自动验收脚本。
2. 自动检查数据版本、哈希、internal ID 唯一性和 manifest 完整性。
3. 检查 discovery/validation 完整覆盖和同源簇泄漏。
4. 检查 PF00741 结果覆盖、失败状态和 receipt 必填字段。
5. 检查 HMM/结构坐标范围和一对一关系是否合理。
6. 在干净环境运行全部测试。
7. 检查脚本是否包含个人绝对路径或隐式环境依赖。
8. 对 MMseqs2 阈值进行敏感性分析，但不替换主划分。

建议敏感性设置：

| 实验 | Identity | Coverage |
| --- | ---: | ---: |
| 宽松 | 0.7 | 0.8 |
| 主方案 | 0.8 | 0.8 |
| 严格 | 0.9 | 0.8 |

比较簇数量、最大簇规模、单序列簇数量和 split 稳定性。敏感性结果不得覆盖 `homology-8b9005e2d9-s42`。

建议分支：`fix/<issue-id>-m2-release-validation`

交付物：

```text
scripts/validate_m2_release.py
tests/test_m2_release.py
results/data_audit/mmseqs_sensitivity.json
reports/M2/m2_reproducibility.md
```

验收标准：

- 自动验收脚本退出码为 0；
- 全部测试通过；
- 缺失关键 receipt 字段时验收必须失败；
- 代码不依赖个人绝对路径；
- 敏感性分析与主划分明确分开；
- validation 不进入后续候选区域发现与调参。

审查：B 审查数据校验逻辑，E 审查统计与输出说明。

### E：结果抽查、图表与解释复核

输入：

- B 的数据摘要；
- C 的 PF00741、HMM 和结构结果；
- D 的验收和敏感性结果。

职责：

1. 分层抽查正常 primary、partial、类型冲突、组合冲突、插入/缺失和长度异常记录。
2. 核对天然序列位置与 HMM match state 映射。
3. 抽查 PDB 7R1C 残基编号和二级结构映射。
4. 生成数据筛选流程图、cohort 数量图、序列长度分布和同源簇规模分布。
5. 生成 discovery/validation、PF00741 覆盖和 HMM occupancy 图。
6. 为每张图记录实验编号、数据版本、split 版本、生成脚本和限制。
7. 整理 M2 汇报时可能被问到的问题和证据位置。

建议分支：`analysis/<issue-id>-m2-visual-summary`

交付物：

```text
reports/M2/figures/
reports/M2/figure_manifest.csv
reports/M2/manual_spot_check.md
```

验收标准：

- 图表完全由脚本生成；
- 图中数字与 JSON/CSV 一致；
- 抽查记录包含 internal ID、检查项和结论；
- 不手工修改最终数字；
- 不把 HMM 覆盖或 occupancy 图解释为功能证明。

审查：D 审查统计，C 审查生物信息学解释，A 审查汇报边界。

## 6. 五天执行安排

### 第 1 天：启动与输入冻结

| 成员 | 当日任务 |
| --- | --- |
| A | 建立总 Issue 和子 Issue；指定负责人、审查人；组织 PR #4 审查；发布验收清单 |
| B | 获取 `gv.zip`；核对哈希；准备复跑环境；开始数据构建 |
| C | 准备 PF00741、HMMER/PyHMMER 和 7R1C；审查输入接口 |
| D | 建立验收脚本框架；冻结 receipt 必填字段 |
| E | 建立图表和人工抽查模板；确认 C 的输出字段 |

当日关口：所有成员能访问同一版本输入；任务 Issue、分支和审查人全部确定；PR #4 至少完成一轮非作者审查。

### 第 2 天：数据复现与真实扫描

| 成员 | 当日任务 |
| --- | --- |
| A | 检查 B/C 是否登记实验；处理输入差异和阻塞问题 |
| B | 完成数据审计复跑；运行 MMseqs2；重新生成 split |
| C | 核验 PF00741 HMM；开始全量扫描；保存日志和原始输出 |
| D | 完成数量、哈希和 split 检查；开始阈值敏感性分析 |
| E | 根据 B 的小型摘要生成数据质量图；准备分层抽查清单 |

当日关口：B 能解释全部数量；C 的扫描输入与 B 的 manifest 一致；失败运行也登记原因。

### 第 3 天：坐标和质量检查

| 成员 | 当日任务 |
| --- | --- |
| A | 组织中期检查；划分必须修复项和限制项 |
| B | 比较新旧输出；提交数据复现 PR；支持异常 ID 排查 |
| C | 完成 PF00741 扫描、HMM 坐标和 7R1C 映射 |
| D | 对 B/C 输出运行自动验收；完成敏感性结果 |
| E | 抽查 HMM/PDB 坐标；生成覆盖和 occupancy 图 |

当日关口：2076 条输入均有明确扫描状态；坐标规则完整；自动验收能主动发现错误输入。

### 第 4 天：修复与非作者复核

| 成员 | 当日任务 |
| --- | --- |
| A | 审查实验台账；更新 M2 报告；检查 Issue/PR 关联 |
| B | 修复或解释数据差异；审查 C 的输入和 ID |
| C | 修复扫描或坐标问题；审查 B 的数据复现 |
| D | 在干净环境重新运行测试；审查 E 的统计图 |
| E | 完成图表和抽查报告；审查 C 的坐标输出 |

当日关口：关键 PR 获得非作者审查；结果、图表和报告一致；不存在未解释的关键差异。

### 第 5 天：结果冻结和 M2 验收

全组共同完成：

1. 从干净环境执行最终复现；
2. 核对 Git commit、数据版本和输出哈希；
3. 冻结 M2 结果；
4. 合并通过审查的 PR；
5. 更新实验台账状态；
6. 完成 M2 报告和限制说明；
7. 由 A 组织验收并关闭 M2。

## 7. 成员交接契约

### B 交给 C

```text
数据版本
metadata manifest
internal ID
sequence hash
sequence-QC 清单
discovery/validation 清单
```

C 核对输入哈希后才开始正式扫描。

### B/C 交给 D

```text
运行 receipt
小型输出摘要
输入哈希
参数与软件版本
失败记录
```

D 不接受只有截图或口头说明的结果。

### B/C/D 交给 E

```text
小型统计表
坐标字段说明
实验编号
数据和 split 版本
限制说明
```

E 不从原始大型文件中手工抄写数字。

### 全体交给 A

```text
已审查 PR
复现命令
测试结果
实验台账
图表 manifest
限制说明
```

A 据此编写最终验收结论。

## 8. GitHub 任务组织

建议建立以下任务：

```text
M2-1：独立复现数据审计和同源划分
M2-2：完成 PF00741 全量扫描和 HMM 坐标
M2-3：建立 M2 自动验收和阈值敏感性
M2-4：完成 M2 图表与人工抽查
M2-5：冻结 M2 报告和实验台账
```

任务状态统一为：

```text
Backlog
Ready
In Progress
Blocked
In Review
Done
```

每日在对应 Issue 中按以下格式更新：

```markdown
今日完成：
当前结果：
阻塞问题：
下一步：
需要谁协助：
对应 commit/结果路径：
```

## 9. 共享存储结构

以下只是建议的共享目录结构，尚无已验证共享位置；不能将此目录示例当成访问证据。实际URI、访问方式和下载哈希登记在收尾PR的 `reports/M2/shared_artifacts.json`：

```text
M2_release/
├── source/
│   └── gv.zip
├── processed/
│   └── gvpa_v1/
├── mmseqs/
├── pfam/
├── structure/
├── receipts/
└── README.md
```

每个目录记录负责人、生成日期、数据版本、Git commit、SHA-256 和生成命令。不得在共享目录或仓库中保存个人访问令牌、密码或GitHub凭据。

## 10. M2 完成定义

以下条件全部满足后，M2 才能关闭：

- [ ] 合并后真实独立复核报告进入 main（不能用助手预检替代）；
- [ ] E运行提交有真实证据；结果归档666de2b不能冒充运行提交；
- [ ] PR #9获非作者真实Review并合并，或团队认可等效替代；
- [ ] 最终待验收完整SHA上全部自动测试通过；
- [ ] 原始输出共享访问与实测SHA256核验完成；
- [ ] A本人在Issue #15留下最终验收评论；

- [ ] 数据审计由非原作者独立复现；
- [ ] MMseqs2 同源划分可复现；
- [ ] 同源簇跨集合泄漏为 0；
- [ ] PF00741 全量扫描真实完成；
- [ ] HMM 与 7R1C 坐标映射完成；
- [ ] 所有正式运行具有 receipt；
- [ ] 实验台账状态与实际结果一致；
- [ ] 关键输出具有 SHA-256；
- [ ] 自动测试与验收脚本通过；
- [ ] 图表与结果文件数字一致；
- [ ] 每个 PR 经过非作者审查；
- [ ] 原始数据和大型文件未进入公共 Git；
- [ ] M2 报告记录方法、结果和限制；
- [ ] A、B、C 共同确认数据定义；
- [ ] A 正式关闭 M2 里程碑。

完成上述验收后，B 的冻结 manifest 可以交给 D 进行 ESM-2 表示，C 的 HMM 坐标可以交给 D/E 进行候选区域和残基证据映射，项目再正式进入 M3。
