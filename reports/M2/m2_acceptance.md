# M2 验收：负责人已确认，阶段关闭

## 2026-09-29 最终验收决定（当前有效）

负责人 A 已在会话亲自确认：“我确认M2验收了”。助手按其明确指令转录此决定，未代签或新增虚假 Review。验收基线为 `567b941e9c80eda5901b094fee1baf72a7e6a730`，确认记录见 [Issue #15](https://github.com/Lithium9767/MLP/issues/15#issuecomment-5889226612)，该任务已关闭。

**M2 已验收，M3 可以推进。** PR #17 已从 Draft 转为待审，PR #7 计划可继续审查；尚未合并 M3 代码或执行 M3 正式实验。旧 dirty 来源、E001 未知运行 SHA 和复现摘要差异等限制继续保留。下方“等待 A/未验收”的文字属于确认前记录，由本节覆盖。

2026-09-29 更新。当前有效规则：取消强制 B/C 独立复核，由 A 本人验收；技术检查不取消。依据负责人本轮明确授权，详见 [CONTRIBUTING](../../CONTRIBUTING.md)。任务为 [Issue #15](https://github.com/Lithium9767/MLP/issues/15)。本报告不是 A 署名，也不补造 Review。

| 关闭条件 | 状态 | 证据 |
| --- | --- | --- |
| C 原始交接实际下载及逐文件哈希 | 已核验 | 整包、12 个文件及三份冻结输入匹配；shared_artifacts.json |
| 干净扫描及代码来源 | 已核验 | #19 已合并；实际运行 1b7b249f986713a3b4c8343b898486be822fc16f，dirty=false；两份实现 Git blob 哈希匹配 |
| 数据与划分 | 已核验 | 2078/2076/1721，478 簇，1453/623；ID/簇交集 0 |
| 原始扫描与坐标/结构统计 | 已核验 | 2076 accepted、0 失败；80219 行逐残基核对；7R1C 88/65/39 |
| 自动测试及 validator | #19 提交已通过，最终提交另记 | 45 passed、2 subtests passed；validator 通过；a_acceptance_precheck.md |
| 队列图与图表清单 | 已核验 | 1721/20/305/30/2；E002 PNG/CSV/manifest 收据哈希匹配 |
| 其他六图 | 已诊断，限制保留 | 原始数据重绘与源统计已核对；不声称像素或旧运行字节复现，figure_download_crosscheck.json |
| E 台账准确 | 已更正并保留未知项 | E002 实际运行 SHA 已登记；E001 运行 SHA 仍为空，不以 archive_commit 冒填 |
| README 与收尾计划 | 已更新 | #9/#16/#18 已合并；本轮规则覆盖历史条款 |
| 强制独立报告、#9 非作者 Review | 不再作为门槛 | 负责人明确取消 M2 强制独立复核；真实历史无 Review 保持不变 |
| A 本人最终验收 | **已确认** | 负责人会话原话已按其指令转录到 #15；不冒称助手独立验收 |

## 必须保留的限制

1. 旧 C scan001 的完整 dirty 工作区仍未知。scan002 是干净代码的替代复现证据，四份原始输出字节与旧收据一致，不证明旧环境已恢复。
2. B 复现摘要有不同 source/split_version；实际 metadata、FASTA、split_manifest 三份字节与冻结记录一致。C 原扫描两份摘要的 CRLF 哈希可由记录提交重建，但并未从包中取得原字节文件。
3. E001 实际运行 SHA 未知；六图诊断不是原运行复现。E002 队列图有完整干净运行来源。
4. 原始上游 gv.zip 和最初 MMseqs TSV 的共享交接不在本次 C 包内；B 的已入库复现收据与冻结划分仍保留，不将它们冒充新取得的上游文件。
5. Release 已公开。上游来源与再分发条件仍由团队确认，助手不代做许可结论。
6. PF00741 是家族坐标，不能当功能标签；validation 未用于本次候选拟合，M3 仍需冻结发现/验证边界。

## 当前决定

技术结果已可交给 A 验收；M2 **尚未代为关闭**。无需等待强制 B/C 独立报告，但 A 应确认上述限制是否可接受。#17 保持 Draft，未运行 M3 正式实验。

证据：[技术核验](a_acceptance_precheck.md)、[下载清单](c_scan_002_handoff.json)、[逐文件核验](c_handoff_technical_verification.json)、[图表诊断](figure_download_crosscheck.json)。历史 engineering_precheck.md、cohort_fix/ 与独立模板保留，不改写历史结论。
