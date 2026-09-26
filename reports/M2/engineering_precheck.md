# M2 助手技术预检（不等于独立审查）

日期：2026-09-26。执行者：Codex 工程助手，未冒用 B/C/A 身份。系统 Windows 11，Python 3.12.11。
远端 API 核实的 main：`5891c2f7a6eb79e48dd7b75175a3af8176e410a2`，与历史 `5891c2f` 为同一提交。以下基线结果不自动覆盖未来最终待验收提交。

## 实际执行

```powershell
python -m unittest discover -s tests -p test_*.py -v
python scripts/validate_m2_release.py --manifest data/processed/gvpa_v1/split/split_manifest.csv --b-data results/data_audit/m2_data_reproduction.json --b-split results/data_audit/m2_split_reproduction.json --structure results/bioinformatics/structure_7r1c_summary.json --pf00741-scan results/bioinformatics/pf00741_scan_summary.json --hmm-coordinates results/bioinformatics/hmm_coordinate_summary.json
```

34 项测试通过；validator 退出码 0。validator 检查摘要、划分一致性，不证明原始文件可下载，也未发现下面的图表错误。

本地文件另以 Python hashlib.sha256 实测，与冻结摘要比较：6 项匹配、7 项文件缺失未验证。详见 [engineering_precheck.json](engineering_precheck.json)。实际 manifest 有 2076 行，discovery 1453、validation 623；重复 ID、跨 split ID、跨 split 簇均为 0。这只证明所查 manifest 的隔离；后续模型拟合过程仍需独立检查。

## 图表抽查：存在阻塞错误

实际查看 `figures/cohort_counts.png`，并对本地已匹配哈希的 metadata.csv 重计数：

| 图中标签 | metadata.analysis_cohort | 源数据计数 | 当前图值 |
| --- | --- | ---: | ---: |
| primary | primary | 1721 | 1721 |
| partial | sensitivity_partial | 20 | 0 |
| type conflict | sensitivity_type_conflict | 305 | 0 |
| partial+type | sensitivity_partial_and_type_conflict | 30 | 0 |
| excluded | excluded_sequence_qc | 2 | 0 |

原因：`scripts/m2_visual_summary.py` 把显示标签作为 Counter 键，缺失键默认为 0。E 应修复字段映射、添加真实队列键的回归测试、先提交代码，再从干净提交重跑图表并登记新的实验编号。不得手工改 PNG 数字。本轮保留原始图和脚本作为待修证据；其他图表的逐条原始输出核对未完成。脚本生成的 `manual_spot_check.md` 不等于真实人工复核。

## 提交来源与共享输出

- `666de2bfb0fb62e660a848fec0c6662d034c076d` 确实存在，PR #14 用它归档 E 脚本和图表；该提交同时加入代码和结果，未提供运行时 HEAD，不能据此认定它是执行代码提交。
- C 扫描 receipt 记录 HEAD `c2a88bd5854859a61c778168db045ccdf9d30253` 且 dirty=true；workflow/coordinates 文件哈希与当前 main 匹配，但不能据此宣称整个 dirty 环境可复现。归档提交是 `ae954d5a43f8ed96e24955a68b276e20a3f75866`。
- [shared_artifacts.json](shared_artifacts.json) 的 expected_sha256 来自已有结果摘要；只有 actual_sha256 非空且 status=match 的条目经过本轮实际文件核对。所有 shared_uri 未提供，访问方式未验证；原负责人须补充固定版本的共享位置、无凭据访问说明及真实下载核验。
- PR #11 未合并而关闭。其四份实现（m2_coordinates、m2_pf00741、download_7r1c、download_hmms）和结构摘要与 main 内容相同；该分支当时只有参考结构结果，完整 2076 条扫描由 PR #13 归档。不能称 #11 的所有文件均被合并。

结论：工程预检完成，存在已知阻塞项；待真实独立复核和最终验收，不建议现在关闭 M2。
