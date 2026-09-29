# M3 B 数据交接

交接基于 M2 冻结版本，不重新随机划分。

| 队列 | 数量 |
| --- | ---: |
| discovery primary | 1202 |
| validation primary | 519 |
| sensitivity | 355 |

敏感性队列分解为 partial 20、类型冲突 305、partial 与类型冲突兼有 30。每条 manifest 记录 `internal_id`、`sequence_id`、序列 SHA-256、长度、cohort、split 和同源簇。

## 本地/共享交接包

`data/processed/m3_b_handoff/m3_b_handoff_bundle.zip` 包含 metadata、聚类 FASTA、原始 split manifest 和三个 M3 manifest。包 SHA-256 见 `results/M3/data_handoff_receipt.json`。该目录被 Git 忽略，不能通过 GitHub 下载；需要将 zip 复制到团队共享存储后，把共享路径补入 receipt。

## 校验

```bash
python scripts/m3_data_handoff.py verify \
  --metadata data/processed/gvpa_v1_reproduction/metadata.csv \
  --split-manifest data/processed/gvpa_v1_reproduction/split/split_manifest.csv \
  --output-dir data/processed/m3_b_handoff
```

脚本拒绝重复或未知 ID、序列 SHA-256/长度错误、metadata 与 split 不一致、同源簇跨集合，以及数量不符。D 的 embedding manifest 可用同一脚本按 cohort 检查漏条、重复、错误哈希和错误 split。
