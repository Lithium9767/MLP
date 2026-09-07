# 数据获取与字段规范

当前没有已下载数据，采集脚本和共享存储地址 TODO。B 负责从经审核的 NCBI/UniProt 来源采集，C 核查许可与标签证据；下载时记录 URL、日期、版本和文件哈希。原始数据存入 data/raw/，只读、不覆盖；处理数据另建版本，位于 data/processed/，均不进入 Git。小型数据卡和划分 ID 可入库。

metadata.csv 的预定义字段如下，尚未创建真实样本表：

```text
sequence_id
accession
database
source_url
retrieval_date
sequence
sequence_length
sequence_hash
valid_residues
organism
taxonomy_id
lineage
pfam_start
pfam_end
label
label_definition
label_source
evidence_level
duplicate_cluster
similarity_cluster
split
exclusion_reason
```

sequence_hash 的规范化与哈希算法 TODO；标签证据等级及 Pfam 坐标基准由 B/C 冻结，不擅自填充。缺失或排除样本保留 exclusion_reason，不能把未知标签自动当负例。
