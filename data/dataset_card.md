# 数据卡（v0.3，候选数据已审计）

## 数据对象

主候选来自 `gv.zip` recognition JSON：2078条按完全一致序列合并的天然GvpA代表记录。数据版本为 `gvpa-recognition-c7f6f005d717`。重新核验得到2078个唯一序列ID和2078个唯一序列哈希。`cluster_statistics` 合计17108条原始记录，其中JSON另列出15030条冗余成员，二者满足 `17108 = 2078 + 15030`。

通过基础序列质量检查、可进入同源聚类的记录为2076条；无partial注释和成员类型冲突的primary队列为1721条。两个含非标准残基 `X` 的序列不进入同源聚类。

## 研究目标

无监督的GvpA潜在功能区域识别与多源验证。当前没有天然序列功能分类标签，不把数据库家族注释或未知值强制转换为正负类。

## 标签调查

- PF00741表示Gas vesicle protein family，只用于身份、同源坐标和区域验证；
- PF01132是EF-P OB domain，不是GvpA家族；
- NCBI/UniProt功能文本缺少足够的样本间差异；
- 文献突变仅在候选区域冻结后作为外部证据。

开发分支报告的2078/2078 PF00741、0/2078 PF01132和数据库覆盖数字必须独立复跑后才能冻结。

## 最低元数据字段

```text
sequence_id
accession
sequence
sequence_hash
sequence_length
organism
description
source_file
is_partial_representative
is_partial_any_member
has_type_conflict
sequence_qc_eligible
primary_analysis_eligible
homology_cluster
analysis_cohort
exclusion_reason
```

## 纳入与排除

- 已复核recognition JSON中不存在重复序列或重复序列ID；
- 2条含非标准残基 `X` 的序列进入序列质量排除队列；
- 代表注释含partial的有17条，任一原始成员注释含partial的代表簇有50个；
- 含非GvpA成员类型的代表簇有335个，均包含GvpJ成员；
- 来源缺失为0条，代表物种缺失为8条；
- partial只匹配注释中的独立单词，不把来源文件名 `GvpA_not_partial_merged.fasta` 误判为partial；
- 每个排除样本保留明确原因；
- 主分析、partial敏感性和类型冲突敏感性分别报告。

## 发现与验证划分

按MMseqs2/CD-HIT等同源簇划分discovery和validation，同簇不得跨集合。无监督研究也必须隔离候选发现和验证；候选区域、参数和阈值在discovery冻结后才可查看validation结果。

## 数据偏差

需要检查长度、C端延伸、物种、partial状态、GvpA/GvpJ冲突和来源文件是否驱动聚类。UMAP图形不作为功能证据。

## 数据来源与许可

最终版本必须记录原始来源、查询/下载日期、数据库版本、许可和SHA-256。许可确认前不在公共仓库重新分发原始序列；Git只提交数据卡、小型清单、划分ID和可复现脚本。

## 版本

`dataset_version=gvpa-recognition-c7f6f005d717`。`split_version` 仍为 `pending_mmseqs2`；得到真实同源簇前不能开始正式discovery/validation比较。输入和输出哈希见 `results/data_audit/gvpa_v1_audit_summary.json`。
