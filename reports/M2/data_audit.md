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

## 当前边界

数据审计已经可复现，但M2尚未完成。当前机器未检测到MMseqs2，WSL实例也无法启动，因此 `split_version=pending_mmseqs2`。正式ESM-2发现/验证实验必须等同源簇和split冻结后再运行。

生成命令、MMseqs2参数和划分命令见 `preprocessing/README.md`。完整计数与哈希见 `results/data_audit/gvpa_v1_audit_summary.json`。
