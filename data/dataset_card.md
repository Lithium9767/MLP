# 数据卡（v0.2，待冻结）

## 数据对象

主候选来自 `gv.zip` recognition JSON：2078条按完全一致序列合并的天然GvpA代表记录。该数字是候选集规模，不是质量控制后的主分析规模。

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
is_partial
valid_residues
gvp_type_conflict
homology_cluster
analysis_cohort
exclusion_reason
```

## 纳入与排除

- 完全重复合并状态需复核并记录代表关系；
- 非标准残基进入排除或单独敏感性队列；
- partial和Gvp类型冲突先标记，不默认删除；
- 每个排除样本保留明确原因；
- 主分析、partial敏感性和类型冲突敏感性分别报告。

## 发现与验证划分

按MMseqs2/CD-HIT等同源簇划分discovery和validation，同簇不得跨集合。无监督研究也必须隔离候选发现和验证；候选区域、参数和阈值在discovery冻结后才可查看validation结果。

## 数据偏差

需要检查长度、C端延伸、物种、partial状态、GvpA/GvpJ冲突和来源文件是否驱动聚类。UMAP图形不作为功能证据。

## 数据来源与许可

最终版本必须记录原始来源、查询/下载日期、数据库版本、许可和SHA-256。许可确认前不在公共仓库重新分发原始序列；Git只提交数据卡、小型清单、划分ID和可复现脚本。

## 版本

冻结时填写 `dataset_version`、`split_version`、输入哈希和变更记录。旧版本不得覆盖。
