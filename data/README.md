# 数据获取与字段规范

当前输入为老师提供的 `gv.zip`，读取其中 `gv/recognition/data/gvp/gvpa/GvpA_sequences.json`。脚本直接读取ZIP成员，不修改或重新打包原文件。输入、成员和输出均记录SHA-256。逐条序列和成员关系位于 `data/processed/`，不进入Git；仓库只保留数据卡和不含序列的审计摘要。

运行：

```powershell
python preprocessing/prepare_gvpa_dataset.py `
  --input-zip C:\path\to\gv.zip `
  --output-dir data/processed/gvpa_v1
```

当前本地数据版本为 `gvpa-recognition-c7f6f005d717`。`metadata.csv` 的冻结字段包括：

```text
sequence_id
internal_id
accession
database
source_url
retrieval_date
sequence
sequence_length
sequence_sha256
organism
description
source_file
representative_gene
representative_gvp_types
member_count
redundant_member_count
member_gvp_types
is_partial_representative
is_partial_any_member
has_type_conflict
has_gvpj_member
source_missing
organism_missing
nonstandard_residues
is_exact_duplicate
exact_duplicate_of
is_length_outlier
sequence_qc_eligible
primary_analysis_eligible
analysis_cohort
exclusion_reason
```

每条序列先删除空白并转为大写，再以UTF-8字节计算SHA-256。`internal_id` 是 `GVPA_000001` 格式的无歧义流水线ID，原始accession完整保存在 `sequence_id`；这可避免MMseqs2重写含 `|` 的UniProt ID。`sequence_qc_eligible` 用于同源聚类输入，`primary_analysis_eligible` 表示无partial、无类型冲突且通过序列质量检查。当前没有功能标签，缺失值不会转换为负例。

当前split版本为 `homology-8b9005e2d9-s42`。逐条 `split_manifest.csv` 保存在本地处理目录；公共仓库只提交不含accession的汇总和哈希。
