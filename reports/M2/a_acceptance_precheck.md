# M2 原始交接技术核验（非 A 署名）

## 2026-09-29 最终验收决定（当前有效）

负责人 A 已在会话亲自确认：“我确认M2验收了”。助手按其明确指令转录此决定，未代签或新增虚假 Review。验收基线为 `567b941e9c80eda5901b094fee1baf72a7e6a730`，确认记录见 [Issue #15](https://github.com/Lithium9767/MLP/issues/15#issuecomment-5889226612)，该任务已关闭。

**M2 已验收，M3 可以推进。** PR #17 已从 Draft 转为待审，PR #7 计划可继续审查；尚未合并 M3 代码或执行 M3 正式实验。旧 dirty 来源、E001 未知运行 SHA 和复现摘要差异等限制继续保留。下方“等待 A/未验收”的文字属于确认前记录，由本节覆盖。

执行人：Codex 工程助手；日期 2026-09-29。被核验 PR #19 HEAD：`dd77543f0580377f3a349e132219cc15d2f8565c`。合并到 main 的提交：`e8b8d2d3db4d61f3449063b9aafe78a32f5e1869`，父基线为 `d64653e4c7b361ada9abf6e3cf09b16ca58fd064`。

## 实际获取与核对

实际下载 [C Release ZIP](https://github.com/Lithium9767/MLP/releases/download/c-m2-scan-002-handoff/MLP_C_M2_scan_002_handoff.zip)，1177773 字节，SHA256 `7c65d1d4ca33918e14aed1d01e0b45b002a90e6e47a55562064ffb382698691e`。包内 12 文件长度与 SHA256 全匹配；原始序列、ZIP 和大输出仅在被忽略的 `data/processed/m2_handoff_received/`，没有提交进 Git。

三份冻结输入字节一致。包内两份 B 复现摘要不同于 scan002 使用的原扫描摘要（复现 split_version 为 homology-3ce9307efa-s42，冻结版为 homology-8b9005e2d9-s42）；核心统计相同，实际 split_manifest 完全相同。原摘要 Git blob 转 CRLF 后哈希与收据一致，记录为重建，不冒称从共享包取得。

重算结果：2078 条 metadata、2076 QC 输入、1721 primary；discovery/validation 为 1453/623，primary 为 1202/519；478 簇无跨集合簇，两个集合 ID 不重叠；2076 accepted、0 failures；80219 行映射逐行残基/序列 SHA/划分一致，39 state 范围合法。PDB 7R1C 88 deposited/65 modeled/39 HMM states，实际 PDB 和结构表哈希匹配。

scan002 receipt 的干净运行提交 `1b7b249f986713a3b4c8343b898486be822fc16f` 存在；两份实现 Git blob 与记录哈希和受审代码完全一致。旧 dirty 环境来源未恢复，保留限制。

## 执行命令和结果

```bash
# 在仓库根目录；ZIP 先从上方实际下载
python scripts/verify_m2_c_handoff.py

python scripts/validate_m2_release.py --manifest data/processed/m2_handoff_received/unpacked/B_input/gvpa_v1_reproduction/split/split_manifest.csv --b-data results/data_audit/m2_data_reproduction.json --b-split results/data_audit/m2_split_reproduction.json --structure results/bioinformatics/structure_7r1c_summary.json --pf00741-scan results/bioinformatics/pf00741_scan_002_summary.json --hmm-coordinates results/bioinformatics/hmm_coordinate_002_summary.json --report ../tmp/c19_validator.json --sensitivity-out ../tmp/c19_sensitivity.json

# Linux/WSL 环境已安装离线依赖并加入 PYTHONPATH
python3 -m pytest tests -q -p no:cacheprovider
git diff --check origin/main...HEAD
```

首次核验执行入口是仓库外 `../m2_handoff_verify.py`；本轮将相同检查入库为 verify_m2_c_handoff.py，便于他人执行。输出为 c_handoff_technical_verification.json，含实测哈希和逐项区别。

完整测试实际使用 `wsl.exe -d Ubuntu --exec python3 /mnt/c/Users/Lenovo/Desktop/ML/tmp/run_m2_linux_tests.py`，启动器配置 PYTHONPATH 后调用 `pytest.main(['tests','-q','-p','no:cacheprovider'])`。Python 3.12.3、pytest 9.1.1、PyHMMER 0.12.3、matplotlib 3.11.2、numpy 2.2.6：**45 passed、2 subtests passed，19.67 秒，退出 0，无跳过**。validator 退出 0；差异检查退出 0。

环境准备曾因 WSL 临时目录被清理和仓库模块路径未加入而失败；修复环境后才得到上述完整通过结果，没有将环境失败或跳过计为成功。未重跑正式天然序列扫描；核验 C 已提供的实际输出，测试中的合成扫描不是 M3 实验。

队列图实测哈希与 E002 收据、图表清单一致。另用 Windows Python 3.13.12/matplotlib 3.10.8 对下载原始数据诊断重绘七图（未覆盖入库图和台账）：队列图字节完全一致；其他六图整体像素不一致，不能称原 E001 字节复现。原图源统计已重算，视觉抽查另附记录；旧 E001 SHA 仍未知。

## 最终验收

本报告是助手技术结果，不是 B/C 独立报告，也不是 A 的签名。最终待验收提交的测试另在 Issue #15 记录。由 A 本人确认证据及 m2_acceptance.md 限制后，才能登记 M2 关闭和启动 M3。
