# M2 C：PF00741、HMM 与 7R1C 审计

状态（更新于 2026-09-27）：**结构映射和 2076 条真实序列扫描的摘要已进入 main；四份扫描原始输出及其固定共享位置仍未取得，扫描 receipt 记录 dirty=true，M2 独立复核未完成。** 关联 [Issue #10](https://github.com/Lithium9767/MLP/issues/10)、[后续审计 Issue #15](https://github.com/Lithium9767/MLP/issues/15) 和 [M2 工作计划](https://github.com/Lithium9767/MLP/blob/main/docs/plans/2026-09-23-m2-team-workplan.md)。本报告区分 C 已亲自运行的结构映射和由 PR #13 归档的全量扫描，不替 A 或独立复核人验收。

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

## B 交接与后续全量扫描

[B 的复现报告](https://github.com/Lithium9767/MLP/blob/feature/m2-1-data-reproduction/reports/M2/data_reproduction_review.md)给出 2078 条候选、2076 条 QC 合格、1721 条 primary、478 个簇和与冻结结果一致的 split manifest SHA-256。其报告注明本地使用 MMseqs2 `18.8cc5c`，而最初冻结记录为 `15-6f452`；原始 cluster TSV 哈希不同，split manifest 哈希相同。C 不把版本差异当成完全同版本复现。

**公共 Git 中只有小型摘要，没有 2076 条序列的三个输入文件。** 2026-09-23 的 C 工作区没有收到实际文件，因此 C 分支未运行全量扫描。随后 [PR #13](https://github.com/Lithium9767/MLP/pull/13) 将真实扫描的 [摘要](../../results/bioinformatics/pf00741_scan_summary.json) 和 [坐标摘要](../../results/bioinformatics/hmm_coordinate_summary.json) 归档到 main：2076 条均为 `accepted`，discovery 1453、validation 623，失败 0、多命中 0；摘要中的三个输入 SHA-256 与冻结记录一致。该归档不等于本轮已独立取得并核验原始文件。

冻结输入哈希：

- B 本地 `data/processed/gvpa_v1_reproduction/metadata.csv`，SHA-256 应为 `62893e31a2dbc245f29c0e9af0a282dcb5099d0e95b96db1a9cead2ac933004e`；
- `sequences_for_clustering.fasta`，SHA-256 应为 `4ff33c04be1d28f42ef5afa7d67261229d19c8e301784ea59750030196e8b122`；
- `split_manifest.csv`，SHA-256 应为 `0a6bd3b44d284a407e1cc9d5be758f0f9439c28395deb6d62d2c787f95341e9c`；
- B 的实际运行 receipt、对应目录位置或带哈希的共享路径，供 C/B 交叉核对。

扫描入口在 HMM 计算前检查数量、ID、序列、划分及哈希；但正式扫描摘要记录代码提交 `c2a88bd5854859a61c778168db045ccdf9d30253` 和 `dirty=true`。该提交尚不包含两个 C 实现文件。现存实现文件字节哈希与 receipt 相符，**原运行时工作区的完整变更及四份原始输出仍不可核验**。详见 [工程预检](cohort_fix/verification.md)；不能把摘要哈希或本地替代文件当成原始输出核验。

当前 C 工作区实际可取得且 SHA-256 匹配的文件：`PF00741.hmm`、`7R1C.pdb`、`7r1c_residue_hmm_map.csv` 以及结构运行的下载收据。四份扫描输出 `failures.json`、`hmm_coordinate_map.csv`、`per_sequence_status.csv`、`raw_hits.jsonl` 和三个 B 输入仍不在本工作区；全部远端分支也未提交这些被 Git 忽略的文件。共享存储的固定 URI 和访问说明仍为空，见 [交付清单](shared_artifacts.json)。

上述五份本地可用文件已整理为工作区外的 `MLP_C_M2_reference_artifacts_20260927.zip`，内含逐文件哈希清单；ZIP SHA-256 为 `abf265a989b25199d8781791b46081bcaa6b8d9a9c539b6850cded1510ac3387`。这是结构参照交接包，**不包含四份扫描原始输出，也不是固定共享 URI**；需由团队上传到实际共享位置并由复核人自行下载核验。

## 验收和限制

本分支恢复了原 C 测试，覆盖冻结交接检查、2076 条合成序列端到端扫描、GA/边界/多命中状态、插入位、不一致残基与 PDB `SEQRES` 完整性；并核对真实 7R1C 的 88/65/39 数量关系及逐表哈希。2026-09-27 运行 `python -m pytest tests -q -p no:cacheprovider`：**45 passed、2 subtests passed**。合成测试无生物学发现。main 已归档的真实扫描摘要不能替代四份原始输出的独立核验；复核人须从固定共享位置取件、实测哈希，或在取得冻结输入后从干净提交重新运行并保留新编号。validation 只可做冻结后覆盖/质量汇总，不得进入后续候选区域或保守性参数选择。
