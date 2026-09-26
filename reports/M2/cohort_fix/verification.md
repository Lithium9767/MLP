# M2 队列图修复与原始输出预检

状态：**修复已完成并提交审查；M2仍待真实独立复核和A验收。** 执行者为Codex工程助手，不是B/C独立复核人。
基线main：`5891c2f7a6eb79e48dd7b75175a3af8176e410a2`。Issue #15；#9/#16仍待Review，#17保持Draft。没有修改或伪造历史审核。

## 队列图：真实重跑

运行编号 `M2-E-FIGURES-002`；实际运行代码SHA：`9c4c1f717e27d90aadbeec4cbe347de0a71e218f`。
运行前git status为空，入口也强制检查干净工作区。Python 3.13.12、Matplotlib 3.10.8、Windows 11；绘图不涉及随机采样，seed为null。

```powershell
python scripts/rerun_m2_cohort.py --metadata data/processed/gvpa_v1/metadata.csv --audit-summary results/data_audit/gvpa_v1_audit_summary.json --run-id M2-E-FIGURES-002
```

以上命令使用已安装的CPython 3.13.12解释器；本机默认MSYS Python 3.12缺matplotlib，首次测试因此报ImportError，切换已有环境后全部通过。复核人应先安装仓库requirements。复跑需新编号（例如003），入口拒绝覆盖002。

| 源字段 analysis_cohort | 源统计 | PNG实际显示 | 结果 |
| --- | ---: | ---: | --- |
| primary | 1721 | 1721 | 一致 |
| sensitivity_partial | 20 | 20 | 一致 |
| sensitivity_type_conflict | 305 | 305 | 一致 |
| sensitivity_partial_and_type_conflict | 30 | 30 | 一致 |
| excluded_sequence_qc | 2 | 2 | 一致 |

最后一类为QC排除类，不是敏感性队列。总数2078。脚本不硬编码这些计数，而是从已核对hash的metadata分组，逐类比较冻结audit，并读取matplotlib真实柱长对照源记录。

五项新增回归检查覆盖：实际柱长与源类别一致、非零源值被绘成零时失败、未知类别失败、dirty工作区拒绝运行、旧运行编号不可覆盖。原错误通过显示标签查Counter键会得到零，现在用显式字段映射修复。原完整绘图入口也调用同一修复函数。

已查看生成PNG：数字正确、标签可读。另重读导出CSV、figure_manifest逐项对照，并核对Git存储的PNG/CSV/manifest字节哈希均与receipt一致。`.gitattributes`仅对有字节证据的产物禁用换行转换，保留真实CRLF并保持其他空白检查，避免跨平台checkout改变已登记hash。

- 图：[cohort_counts.png](../figures/cohort_counts.png)
- 逐类统计：[cohort_counts.csv](../runs/M2-E-FIGURES-002/cohort_counts.csv)
- 命令、SHA、版本、输入/输出hash：[run_receipt.json](../runs/M2-E-FIGURES-002/run_receipt.json)
- 图表清单：[figure_manifest.csv](../figure_manifest.csv)；仅cohort行更新到002，另外六图未重跑。
- 图SHA256：`cafd11ab1db91fd5020060c48b838cfb5e4b54e23ea7fb246796bbef3ad12021`。

旧 `M2-E-FIGURES-001` 的git_commit继续为空；`666de2bfb0fb62e660a848fec0c6662d034c076d`仅是归档提交。新运行不补造旧运行来源，也不证明其他六图已复核。

## C：真实文件与dirty来源

```powershell
python scripts/audit_m2_artifacts.py --pdb-candidate ../github_remote_snapshots_2026-09-21/lrf/Lithium9767-MLP-d272ad7/results/structure/7R1C.pdb --out reports/M2/cohort_fix/c_artifact_audit.json
```

记录：[c_artifact_audit.json](c_artifact_audit.json)。这次检查17个receipt条目：4个match、2个mismatch、11个未取得。HMM及下载receipt按扫描/结构两个运行分别列出，因此17是条目数，并非17个不同文件。每项保留expected/actual、来源和实际检查路径；未取得项actual=null。

| 项目 | 实际状态 |
| --- | --- |
| metadata、FASTA、split manifest | 三项实际文件字节hash匹配 |
| 7R1C替代本地快照 | 实际字节hash匹配；不是receipt原路径，也不是共享访问证明 |
| 两份audit/split摘要 | 当前字节hash不匹配；转换为CRLF后与receipt完全匹配。只作换行诊断，未篡改源文件或冒充原始输入 |
| 四份PF00741输出 | failures.json、hmm_coordinate_map.csv、per_sequence_status.csv、raw_hits.jsonl均未取得 |
| HMM、下载receipt、原路径PDB与结构坐标表 | 未验证，逐项见JSON |
| 共享位置与访问方式 | 未提供；receipt只有local/shared占位描述，Issue #15无访问说明 |

扫描receipt的HEAD `c2a88bd5854859a61c778168db045ccdf9d30253`存在，但该提交本身没有m2_pf00741.py和m2_coordinates.py；运行时dirty新增/变更不可仅凭HEAD复现。目前工作区这两个文件的真实字节hash与receipt匹配；Git中的LF blob hash不同，不混淆二者。未取得原始dirty补丁、完整环境锁定或原运行工作区，因此仍标记来源未验证。

## 工程验证与关闭门槛

```powershell
python -m unittest discover -s tests -p test_*.py -v
python scripts/validate_m2_release.py --manifest data/processed/gvpa_v1/split/split_manifest.csv --b-data results/data_audit/m2_data_reproduction.json --b-split results/data_audit/m2_split_reproduction.json --structure results/bioinformatics/structure_7r1c_summary.json --pf00741-scan results/bioinformatics/pf00741_scan_summary.json --hmm-coordinates results/bioinformatics/hmm_coordinate_summary.json
git diff --check origin/main...HEAD
```

最终提交上的执行结果与完整SHA记录在本修复PR说明及验证日志中，不把receipt字段校验等同于原始文件核验。

A须先核实账号及贡献范围，再指派未编写受审改动的真实B/C。复核人必须亲自获取文件、计算hash、检查图表、运行全套检查，以自己的账号提交署名报告。助手不在independent_review.md签名，也不代发Review或A验收。

只有修复/证据和真实独立报告进入main、#9获非作者Review并合并、最终SHA测试通过、全部证据阻塞处理完毕且A在#15留下真实验收结论后，M2才能关闭；此前#17保持Draft，不执行M3正式实验。
