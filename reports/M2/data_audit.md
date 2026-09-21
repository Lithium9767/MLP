# M2 数据审计记录

## 冻结对象

- 数据版本：`gvpa-recognition-c7f6f005d717`
- 来源：老师提供的 `gv.zip`
- ZIP成员：`gv/recognition/data/gvp/gvpa/GvpA_sequences.json`
- 候选代表序列：2078
- 原始记录合计：17108（2078条代表 + 15030条冗余成员）

逐条序列和成员表只保存在本地 `data/processed/gvpa_v1/`。公共仓库只提交脚本、字段定义及不含序列的汇总。

## 审计结论

| 项目 | 数量 | 处理方式 |
| --- | ---: | --- |
| 唯一序列ID | 2078 | 通过 |
| 唯一序列哈希 | 2078 | 通过 |
| 非标准残基序列 | 2 | 不进入同源聚类 |
| 可进入同源聚类 | 2076 | 用于建立全体同源簇 |
| primary队列 | 1721 | 后续主分析 |
| 代表注释partial | 17 | 敏感性队列 |
| 任一成员注释partial | 50 | 敏感性队列 |
| 含非GvpA成员类型 | 335 | 敏感性队列 |
| 代表物种缺失 | 8 | 保留缺失标记 |
| 来源缺失 | 0 | 无 |

`partial` 必须匹配注释文本中的独立单词。来源文件 `GvpA_not_partial_merged.fasta` 中的字符串不能作为partial证据；这解释了早期简单字符串搜索产生的814/938误计数。

## 同源划分

在WSL Ubuntu 24.04中使用MMseqs2 15-6f452，参数为 `--min-seq-id 0.8 -c 0.8 --cov-mode 1`。第一次运行暴露出MMseqs2会重写8个含 `|` 的UniProt FASTA ID，划分校验因此主动失败。流水线改用 `GVPA_000001` 格式内部ID后重新提取和聚类，原始accession继续保存在metadata中。

最终478个簇完整覆盖2076条序列，冻结 `split_version=homology-8b9005e2d9-s42`：

| 集合 | 全部序列 | primary序列 | 同源簇 |
| --- | ---: | ---: | ---: |
| discovery | 1453 | 1202 | 335 |
| validation | 623 | 519 | 143 |

校验未发现未知ID、缺失ID、重复归属或同簇跨集合。生成命令见 `preprocessing/README.md`，运行参数及结果哈希见 `results/data_audit/`。M2的数据审计与同源划分已完成；PF00741/HMM坐标仍由C后续冻结。
