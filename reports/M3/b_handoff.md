# M3 B 数据交接：工程预检与待办

状态：**交接包技术预检通过，B 本人的任务提交与非作者审查待完成**。本报告由 Codex 工程助手整理，不代表 B 签署，也不代表 A 验收。沿用冻结数据版本 `gvpa-recognition-c7f6f005d717`、split 版本 `homology-8b9005e2d9-s42`；不重新划分数据。

## 已收到的内容

本机收到 `m3_b_handoff.zip`，外层 SHA-256 为 `d4d5a53c78c10eb4bfc61008ec17eae225b9cab6305c69b7ed194de7d63307a6`。内层包 SHA-256 为 `5bb7649e1c8f6de6830e71dc7d024157e7055257cf3add28f22d78f067ada18c`。两层均可读取且 CRC 检查通过；外层三份 manifest 与内层同名文件逐字节相同。完整六文件哈希见[运行收据](runs/M3-B-HANDOFF-CHECK-001/run_receipt.json)。包仅在本机读取，未将序列或大文件提交 Git。

| 检查项 | 实测结果 |
| --- | --- |
| metadata | 2078 条，primary 1721、type conflict 305、partial + conflict 30、partial 20、QC 排除 2 |
| 可用 FASTA 与 split | 各 2076 条；ID、序列、长度、序列 SHA-256、cohort 一致 |
| split | discovery 1453 / validation 623 |
| 同源簇 | discovery 335 / validation 143；跨集合重叠 0 |
| 主队列清单 | discovery 1202 / validation 519 |
| 敏感性清单 | 355；discovery 251 / validation 104 |
| 冻结输入哈希 | metadata、split、FASTA 与 M2 已归档 SHA-256 一致 |

实际命令（在仓库根目录）：

```powershell
python scripts/verify_m3_b_handoff.py --bundle ..\..\m3_b_handoff.zip --receipt reports/M3/runs/M3-B-HANDOFF-CHECK-001/run_receipt.json
python -m unittest discover -s tests -p test_m3_b_handoff.py -v
```

正式预检从干净提交 `3167fe440e9224690e50010fafff0c5f7070fd07` 执行，收据中的 `run_commit` 指向该提交，不指向后续归档提交。脚本核对 ZIP 内容、M2 冻结哈希、逐条序列哈希、manifest 分区、同源簇隔离；没有运行 validation 模型或据此选择候选区域。四项本任务测试通过，覆盖哈希不符、manifest 错配和同源簇泄漏。

`python scripts/validate_m2_release.py` 通过。尝试执行 `python -m unittest discover -s tests -v` 时，本机 Python 环境缺少 `matplotlib` 和 `pytest`，两个既有测试模块在导入阶段失败；这不构成全量测试通过的证据。应在装有 `requirements.txt` 与 M2 C 依赖的环境中重跑并保存真实结果。

## B 仍需亲自完成

1. 确认交接包制作者、生成代码与命令、运行环境、真正的运行提交、数据来源。仅凭包名不能核实贡献身份或独立复现。
2. 在 B 自己的账号和任务分支提交或关联数据构建方法、完整性检查记录和数据/split receipt；说明这些文件是 M2 冻结输入的重打包，还是 B 独立重跑的输出，不得把助手预检冒作 B 的运行。
3. 提供团队共享存储中原始包的固定位置、访问方式和该位置文件的实际哈希；目前只核对了本机收到的 ZIP，**共享副本及其权限未验证**。不要上传未经许可的原始序列或凭据。
4. 请 C 核对交给其坐标流程的 ID/序列接口，A 核对版本与 discovery/validation 边界；按 `CONTRIBUTING.md` 由非作者真实成员审查对应 PR。B 的 M3 任务只有在证据、身份和审查齐备后才标记完成。

当前包可以供 D 使用 **discovery primary** 输入；validation primary 清单只做完整性核对，候选冻结和 A 批准前不运行验证推理或调参。
