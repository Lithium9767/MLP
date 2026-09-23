# M2 B 角色数据复现审查

## 运行信息

- 负责人：B
- 数据版本：`gvpa-recognition-c7f6f005d717`
- 输入：`gv/gv/recognition/data/gvp/gvpa/GvpA_sequences.json`
- MMseqs2：`18.8cc5c`
- 参数：`--min-seq-id 0.8 -c 0.8 --cov-mode 1`
- split seed：`42`

## 复现结果

| 项目 | 冻结基线 | 本次复现 | 结论 |
| --- | ---: | ---: | --- |
| 候选序列 | 2078 | 2078 | 一致 |
| QC 合格 | 2076 | 2076 | 一致 |
| primary | 1721 | 1721 | 一致 |
| 同源簇 | 478 | 478 | 一致 |
| discovery | 1453 / 335 簇 | 1453 / 335 簇 | 一致 |
| validation | 623 / 143 簇 | 623 / 143 簇 | 一致 |
| cluster leakage | 0 | 0 | 一致 |

数据输出哈希与冻结数据审计一致，split manifest 哈希也一致。MMseqs2 原始 cluster TSV 哈希不同，原因是本机使用 `18.8cc5c`，冻结记录使用 `15-6f452`；该差异已记录在 receipt，不能把本次结果标记为完全同版本复现。

## 交接给 C

C 可使用 `data/processed/gvpa_v1_reproduction/` 下的 metadata、FASTA 和 split manifest。正式扫描前应核对 `metadata_sha256`、FASTA 哈希和本报告中的数据版本；PF00741 扫描输入应保持 2076 条 QC 合格序列。

## 限制

工作区只有解压后的 `gv/`，没有原始 `gv.zip`，因此无法核对原始 ZIP 文件哈希；本次使用 recognition JSON 的 SHA-256 作为数据内容身份。
