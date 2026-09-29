# M3 A/D 工程验证记录

执行者：Codex 工程助手。非成员独立 Review，非 A 最终验收。

- 全量自动测试：55 passed、2 subtests passed，20.37秒，退出0。包括M2真PyHMMER流程、新窗口坐标、NaN、跨split毒化、重复坐标、候选批准锁和无映射簇不捏造state1。
- M2 release validator：退出0，仍使用真实下载包中的冻结manifest与scan002/coordinate002收据。
- `git diff --check`：通过。M2 preprocessing、bioinformatics及validator实现未被覆盖。
- 六个运行目录的outputs_sha256逐一从实际本地文件重算并匹配，原始JSON/CSV字节归档时使用.gitattributes的M3 runs规则保留。
- cache smoke和readback的24条embedding实际SHA逐项一致；追加敏感性队列复用1202条精确缓存，只对额外251条推理。
- PCA/UMAP图已实际打开检查；UMAP仅展示固定3000窗口子样本，不用于聚类。

## 命令

```text
wsl.exe -d Ubuntu --exec python3 /mnt/c/Users/Lenovo/Desktop/ML/tmp/run_m2_linux_tests.py
python scripts/validate_m2_release.py --manifest data/processed/m2_handoff_received/unpacked/B_input/gvpa_v1_reproduction/split/split_manifest.csv --b-data results/data_audit/m2_data_reproduction.json --b-split results/data_audit/m2_split_reproduction.json --structure results/bioinformatics/structure_7r1c_summary.json --pf00741-scan results/bioinformatics/pf00741_scan_002_summary.json --hmm-coordinates results/bioinformatics/hmm_coordinate_002_summary.json --report ../tmp/m3_preserved_m2_validator.json --sensitivity-out ../tmp/m3_preserved_m2_sensitivity.json
git diff --check
git diff origin/main -- preprocessing bioinformatics scripts/validate_m2_release.py
```

正式运行的完整参数在每份runs/<experiment_id>/run_receipt.json的provenance.command，执行前工作树干净，实际producer SHA分别记录在收据与台账。后续归档提交不是运行代码SHA。

WSL测试环境：Python3.12.3、pytest9.1.1、PyHMMER0.12.3、numpy2.2.6。Windows实测环境：Python3.13.12、torch2.11.0+cpu、transformers4.57.6、numpy2.4.4、scikit-learn1.8.0、hdbscan0.8.44、umap-learn0.5.12、matplotlib3.10.8。测试与实验的不同环境不混为同一环境。

首次新增组成基线测试因float32约1.19e-7舍入差失败，已采用1e-6数值容差后通过；未计作首次即通过。归档助手遇到GBK读取中文文档的编码问题后改为UTF-8，不改变实验输出。

大数组、模型、全窗口表与joblib只在忽略目录；本PR只公开代码、配置、收据、小型统计与图。共享交接仍需建立，不能把本地artifact_uri称为团队已可访问的地址。
