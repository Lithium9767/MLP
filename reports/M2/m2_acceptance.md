# M2 验收

数据版本 `gvpa-recognition-c7f6f005d717`。划分版本 `homology-8b9005e2d9-s42`。本报告汇总已经进入 `main` 的 B/C/D 结果和 E 的图表抽查。PF00741 命中不是功能标签。

| 关闭条件 | 状态 | 证据 |
| --- | --- | --- |
| 数据审计可复现 | 完成 | `results/data_audit/m2_data_reproduction.json`：2078 / 2076 / 1721，输出哈希与冻结审计一致 |
| MMseqs2 同源划分可复现 | 完成 | 478 簇。`split_manifest.csv` 哈希与冻结划分一致。MMseqs2 `15-6f452` 的原始 cluster TSV 哈希不同，已记入复现说明 |
| 簇泄漏为 0 | 完成 | 冻结 split 与 B 的复现 receipt 均为 `cluster_leakage=false` |
| PF00741 全量扫描 | 完成 | `results/bioinformatics/pf00741_scan_summary.json`：2076 条 `accepted`，discovery 1453，validation 623 |
| HMM 与 7R1C 坐标 | 完成 | `hmm_coordinate_summary.json`；7R1C 链 N，88 个提交残基，65 个建模，39 个 match state |
| 正式运行有 receipt | 完成 | 审计、划分、扫描、结构、敏感性 JSON 均含哈希或版本 |
| 台账与结果一致 | 完成 | `experiments/registry.csv` 已去掉 `split_pending`，并登记 C/D/E |
| 自动验收通过 | 完成 | `scripts/validate_m2_release.py` 对已发布 receipt 退出码为 0。测试见 `tests/test_m2_release.py` |
| 图表与数字一致 | 完成 | `reports/M2/figures/` 与 `reports/M2/manual_spot_check.md`，由 `scripts/m2_visual_summary.py` 生成 |
| 原始大文件不入库 | 完成 | 序列、PDB、HMM 和 MMseqs2 中间文件留在 `data/processed/`、`data/raw/`、`work/` |
| 敏感性不覆盖主划分 | 完成 | `results/data_audit/mmseqs_sensitivity.json`。identity 0.8 的 manifest 与冻结划分一致；0.7 与 0.9 只作对照 |

仍需人工完成、不能由这份报告代替的两项：

- 本验收提交的 PR 需要一名非作者审查。
- A 在审查通过并合并后关闭 M2。
