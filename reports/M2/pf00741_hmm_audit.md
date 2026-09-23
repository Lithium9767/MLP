# M2 C：PF00741、HMM 与 7R1C 审计

状态：**真实 7R1C 参照映射完成；2076 条天然序列全量扫描等待 B 的实际输入文件。** 日期 2026-09-23。关联 [Issue #10](https://github.com/Lithium9767/MLP/issues/10) 和 [M2 工作计划](https://github.com/Lithium9767/MLP/blob/docs/8-m2-team-workplan/docs/plans/2026-09-23-m2-team-workplan.md)。本报告只陈述 C 角色的核查与结果，不替 A 关闭 M2。

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

## B 交接的已知事实与缺失文件

[B 的复现报告](https://github.com/Lithium9767/MLP/blob/feature/m2-1-data-reproduction/reports/M2/data_reproduction_review.md)给出 2078 条候选、2076 条 QC 合格、1721 条 primary、478 个簇和与冻结结果一致的 split manifest SHA-256。其报告注明本地使用 MMseqs2 `18.8cc5c`，而最初冻结记录为 `15-6f452`；原始 cluster TSV 哈希不同，split manifest 哈希相同。C 不把版本差异当成完全同版本复现。

**公共 Git 中只有这些小型摘要，没有 2076 条的 `metadata.csv`、`sequences_for_clustering.fasta`、`split_manifest.csv`。** 当前 C 工作区及常见用户文件目录也未找到真实文件，无法核对 B 的每个 internal ID/sequence hash，不能运行或声称完成全量扫描。所需交接：

- B 本地 `data/processed/gvpa_v1_reproduction/metadata.csv`，SHA-256 应为 `62893e31a2dbc245f29c0e9af0a282dcb5099d0e95b96db1a9cead2ac933004e`；
- `sequences_for_clustering.fasta`，SHA-256 应为 `4ff33c04be1d28f42ef5afa7d67261229d19c8e301784ea59750030196e8b122`；
- `split_manifest.csv`，SHA-256 应为 `0a6bd3b44d284a407e1cc9d5be758f0f9439c28395deb6d62d2c787f95341e9c`；
- B 的实际运行 receipt、对应目录位置或带哈希的共享路径，供 C/B 交叉核对。

拿到这些文件后，C 的真实扫描入口会在任何 HMM 计算前拒绝数量、ID、序列、划分或哈希不一致的输入。应完成的两个小型摘要是 `results/bioinformatics/pf00741_scan_summary.json` 和 `hmm_coordinate_summary.json`；本次没有生成占位文件或合成数字。

## 验收和限制

`python -m pytest tests -q -p no:cacheprovider` 为 **28 passed、2 subtests passed**。新增测试覆盖冻结交接检查、2076 条合成序列端到端扫描、GA/边界/多命中状态、插入位、不一致残基与 PDB `SEQRES` 完整性；并对真实 7R1C 输出核对 88/65/39 的数量关系及逐表哈希。合成测试无生物学发现。真实 PF00741 扫描尚未完成，因此 discovery/validation 覆盖、天然序列的 HMM 映射、D/E 的抽查和最终 M2 结果冻结仍未验收。validation 只可做冻结后覆盖/质量汇总，不得进入后续候选区域或保守性参数选择。
