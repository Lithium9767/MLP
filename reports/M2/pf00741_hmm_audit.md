# M2 C：PF00741、HMM 与 7R1C 审计

## 2026-09-29 交接状态更新

PR #19 已合并。助手已实际下载公开 Release，核对整包及 12 个文件哈希和全部坐标；见 [技术核验报告](../../reports/M2/a_acceptance_precheck.md)。负责人已取消 M2 强制 B/C 独立复核，改由 A 本人最终验收（Issue #15）。以下历史“待非作者核验/独立复核”不再是当前阻塞，也不表示有人补签 Review。旧 dirty 环境和旧 E001 运行 SHA 未知仍保留。

状态（更新于 2026-09-27）：**C 已从干净提交亲自完成 2076 条真实序列全量扫描，四份原始输出和结构参照已打包并作为 GitHub Release 附件共享，摘要与运行收据提交到本分支；非作者独立核验及 A 的 M2 验收仍待完成。** 关联 [Issue #10](https://github.com/Lithium9767/MLP/issues/10)、[后续审计 Issue #15](https://github.com/Lithium9767/MLP/issues/15) 和 [M2 工作计划](https://github.com/Lithium9767/MLP/blob/main/docs/plans/2026-09-23-m2-team-workplan.md)。这次运行不冒充 PR #13 的历史 dirty 运行，也不替 A 或独立复核人验收。

## 输入、版本与结构结果

已从 [EBI InterPro](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00741?annotation=hmm) 下载真实 `PF00741.24` HMM：SHA-256 `e44f933c618426c2d0465657a64f1890f50feaeaea626ecb0b8dedb734f38c2f`，39 个 match state，模型自带 GA 阈值 25/25 bits。已从 [RCSB](https://files.rcsb.org/download/7R1C.pdb) 下载真实 PDB 7R1C：SHA-256 `c83776e7091d1342a4c6dbc40fbcc7be2e2beff87b12f1c52c9d6de5be0be870`。两者下载日期为 2026-09-23，下载收据及原始文件保存在被 Git 忽略的 `data/raw/`。

对 RCSB 的作者链 `N` 建立了 1-based 提交序列位置、实验建模残基编号和 PF00741 match state 的逐位映射。88 个提交位置中 65 个有实验建模的 CA，23 个未建模；39 个位置与 39 个 HMM match state 对齐。局部未对齐的建模残基保留记录，不能把 HMM 不覆盖等同于功能缺失。结构参照没有独立提供天然序列功能标签。真实运行绑定干净代码提交 `a45930c04b4d66be7d34e0a8139b4de4d953ad97`，登记实验 `M2-C-7R1C-001`；逐残基表 SHA-256 为 `c48bd6f1c9f489000872ab3aa50f56d30106c0b56cf3762dbda4ef8a4aa56b7c`。具体复现命令及逐表口径见 [运行说明](../../bioinformatics/M2_PF00741.md) 和 [可入库结构摘要](../../results/bioinformatics/structure_7r1c_summary.json)。

## PR #4 审查与修正

对 [PR #4](https://github.com/Lithium9767/MLP/pull/4) 的 C 相关模块已提交[非作者复核意见](https://github.com/Lithium9767/MLP/pull/4#pullrequestreview-5287185460)，要求修改后再作为 M2 正式入口：

1. 原 `hmm_coordinates.scan_sequences` 用默认搜索阈值且只保留最高分的一个域，无法交付 GA、多命中和边界状态。
2. 原 `write_mapping_outputs` 对全部输入计算保守性，没有 discovery/validation 隔离，会让 validation 提前参与 M3 候选发现。
3. 原结构入口只输出局部比对命中的 39 个建模残基，未显示其余 26 个已建模残基与 23 个未建模提交位点，也没有核对比对残基身份。
4. 原下载入口仅检查 HMMER 文件开头，未将 `PF00741.24` 和 GA 阈值核入收据。

C 从 PR #4 选择性保留比对列与 PDB 的解析思想，建立独立的 [严格运行入口](../../bioinformatics/m2_pf00741.py)：先核对冻结输入哈希、完整序列和 split，再分别做 GA 与宽松搜索，保留多域、每条状态、插入/缺失、分组覆盖和来源收据。结构映射使用完整 88 位提交序列，未建模位置标记清楚。测试和真实 7R1C 运行已验证该入口；PR #4 仍需作者修改或 A 指定经审查的等价入口。

## B 输入与 C 的干净全量扫描

[B 的复现报告](https://github.com/Lithium9767/MLP/blob/feature/m2-1-data-reproduction/reports/M2/data_reproduction_review.md)给出 2078 条候选、2076 条 QC 合格、1721 条 primary、478 个簇和与冻结结果一致的 split manifest SHA-256。其报告注明本地使用 MMseqs2 `18.8cc5c`，而最初冻结记录为 `15-6f452`；原始 cluster TSV 哈希不同，split manifest 哈希相同。C 不把版本差异当成完全同版本复现。

2026-09-27 收到 B 的 `gvpa_v1_reproduction.zip`，ZIP SHA-256 为 `ac6366962bf5e14872127c1383fe27dc62097cd5fa6e0861b598697ce85e9632`。从 ZIP 中实际读取并核对三份文件；哈希均与冻结记录完全一致，未使用测试数据或重新下载的替代序列：

冻结输入哈希：

- `data/processed/gvpa_v1_reproduction/metadata.csv`：`62893e31a2dbc245f29c0e9af0a282dcb5099d0e95b96db1a9cead2ac933004e`；
- `sequences_for_clustering.fasta`：`4ff33c04be1d28f42ef5afa7d67261229d19c8e301784ea59750030196e8b122`；
- `split/split_manifest.csv`：`0a6bd3b44d284a407e1cc9d5be758f0f9439c28395deb6d62d2c787f95341e9c`。

扫描入口在 HMM 计算前检查数量、ID、序列、划分及哈希。C 在代码提交 `1b7b249f986713a3b4c8343b898486be822fc16f`、工作区干净时执行新编号 `pf00741_scan_002`；新 [扫描收据](../../results/bioinformatics/pf00741_scan_002_summary.json)记录 `dirty=false`、Python 3.14.3、PyHMMER 0.12.3、PF00741.24、GA 25/25 bits、完整命令与输入/输出哈希。真实结果：2076 条 `accepted`、0 失败、0 多命中；discovery 1453、validation 623。[坐标收据](../../results/bioinformatics/hmm_coordinate_002_summary.json)记录 80,219 行坐标、2076 条映射序列、0 插入，并分别汇总两组覆盖。此命中只说明 PF00741 家族坐标，不是功能正负标签。

四份新生成的原始输出 `failures.json`、`hmm_coordinate_map.csv`、`per_sequence_status.csv`、`raw_hits.jsonl` 均已实际读取并计算 SHA-256；四项哈希**逐一等于** PR #13 的历史摘要记录。这证明新干净运行复现了相同字节的输出；不证明旧运行的完整 dirty 工作区来源已查明。两份结构参照和 HMM 原始文件也已核验。

原始 B 输入 ZIP、四份扫描输出、两份新摘要、HMM/PDB、下载收据和结构逐残基表已整理为 [GitHub Release 交接包](https://github.com/Lithium9767/MLP/releases/download/c-m2-scan-002-handoff/MLP_C_M2_scan_002_handoff.zip)，ZIP SHA-256 为 `7c65d1d4ca33918e14aed1d01e0b45b002a90e6e47a55562064ffb382698691e`。C 已从公开链接重新下载并核对实际字节；[交接清单](c_scan_002_handoff.json)列出每个文件的长度和 SHA-256。**非作者下载核验仍待完成**；原始序列和大输出不进入 PR #19 的 Git 提交。Release 附件是公开的，数据卡要求的来源与再分发条件仍需团队确认。旧的结构-only ZIP 保留为历史包，不再作为完整 C 交接。

PR #13 的旧扫描收据记录 `c2a88bd5854859a61c778168db045ccdf9d30253` 和 `dirty=true`；该提交尚不包含两个 C 实现文件。旧运行时工作区的完整变更仍未取得，见 [工程预检](cohort_fix/verification.md)。本次干净运行提供独立、可复查的替代证据，不改写旧收据。

## 验收和限制

本分支恢复了原 C 测试，覆盖冻结交接检查、2076 条合成序列端到端扫描、GA/边界/多命中状态、插入位、不一致残基与 PDB `SEQRES` 完整性；并核对真实 7R1C 的 88/65/39 数量关系及逐表哈希。2026-09-27 运行 `python -m pytest tests -q -p no:cacheprovider`：**45 passed、2 subtests passed**。合成测试无生物学发现。新扫描使用真实 B 输入且输出已实际核对；M2 最终验收仍需要非作者从共享位置下载文件复核、检查图表及 A 亲自确认。validation 只可做冻结后覆盖/质量汇总，不得进入后续候选区域或保守性参数选择。
