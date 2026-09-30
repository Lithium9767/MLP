# M3 B 交接包：助手技术复核

复核人：Codex 工程助手。复核对象：[Draft PR #23](https://github.com/Lithium9767/MLP/pull/23) 的数据交接校验及本机收到的 `m3_b_handoff.zip`。此文件不署名为成员 B，也不是 GitHub 上真实非作者成员的 Approve Review。

## 结论

**本机 ZIP 的数据完整性与冻结 split 关系技术核对通过；B 的正式交接仍未验收。** 外层 ZIP SHA-256 为 `d4d5a53c78c10eb4bfc61008ec17eae225b9cab6305c69b7ed194de7d63307a6`。内层三个清单与外层副本字节相同，metadata、split、FASTA 的 SHA-256 与 M2 冻结记录相同。未发现 sequence hash、ID 对应、cohort、manifest 分区或同源簇泄漏错误。

我另用一次性 CSV/ZIP 读取代码独立重算，而非调用 PR 的校验函数，得到：metadata 2078、eligible/split 2076、primary 1721、敏感性 305/30/20、QC 排除 2；discovery/validation 1453/623，335/143 个同源簇，跨集重叠 0；三份 manifest 合计恰好覆盖 2076 个不重复 ID；逐条序列 SHA-256 错误数 0。该独立计算只作技术交叉核对，不能证明 B 的本人复跑。

## 代码审查与验证

首轮代码审查发现 ZIP 条目用集合比较会忽略重复的同名成员。已在 `f414e319775dfc15d88c0ec67d4cc146014881eb` 修复并增加对应失败测试，保留原 001 收据。修复后从干净提交重新运行，产生 [M3-B-HANDOFF-CHECK-002 收据](runs/M3-B-HANDOFF-CHECK-002/run_receipt.json)，其 `run_commit` 为该完整 SHA。

执行与结果：

```powershell
python -m unittest discover -s tests -p test_m3_b_handoff.py -v
# 5 tests passed
python scripts/verify_m3_b_handoff.py --bundle ..\..\m3_b_handoff.zip --receipt reports/M3/runs/M3-B-HANDOFF-CHECK-002/run_receipt.json
# passed; run_commit=f414e319775dfc15d88c0ec67d4cc146014881eb
python scripts/validate_m2_release.py
# M2 release validation passed
git diff --check
# passed
```

此前尝试全量 `python -m unittest discover -s tests -v`：有两个既有测试模块因本机缺少 `matplotlib` 和 `pytest` 而在导入阶段失败，故**全量测试未获通过结论**。不能把这些导入失败描述为代码逻辑通过或失败；需在完整依赖环境重跑。

## 尚不能代替 B 确认的事实

- ZIP 的实际制作人、制作命令、运行代码提交及其与成员 B 的账号关系未核实。包与 M2 冻结文件相同，只能证明交接字节一致，不能证明 B 独立重建了数据。
- 共享存储中的原始交接包位置、访问权限与该处文件哈希未取得；目前验证的仅是本机副本。
- C 的坐标接口审查、真实非作者 PR Review 尚无本报告可代替的证据。

因此 [Issue #24](https://github.com/Lithium9767/MLP/issues/24) 的上述待办仍应开放；A 可确认本技术预检，但不可把本报告标为 B 本人的运行或复核。没有对 validation 集运行 M3 模型推理、拟合或选参。
