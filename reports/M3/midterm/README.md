# M3 中期汇报与项目成果索引

关联任务：[#6](https://github.com/Lithium9767/MLP/issues/6)。本目录汇报材料由 A 统稿，工程助手整理，不代替成员署名或验收。

## 汇报文件

- [PPT v2](GvpA_M3_midterm_v2.pptx)：9 页主讲，约 6 分钟；另有 3 页默认隐藏的问答备用页。
- [对应讲稿、要求核对与问答](speaker_notes_v2.md)：适配 2–4 分钟问答。
- [文件大小与 SHA-256](manifest.json)：对本目录实际文件计算。

汇报人：史晨皓。成员：A 史晨皓、B 韩涛鸿、C 燕相楠、D 李若凡、E 杨傲翔。首页不强调组长身份。

## 已有成果位置

截至本次整理，M2 已由负责人确认验收。M3 代码与结果在 [PR #22](https://github.com/Lithium9767/MLP/pull/22)，其核实后的提交为 `7bb22fcdb30363e7b0116dd0291835cb07f9b13b`，尚未合并。以下固定链接指向该证据版本：

- [ESM-2 表征入口](https://github.com/Lithium9767/MLP/blob/7bb22fcdb30363e7b0116dd0291835cb07f9b13b/features/esm2_embed.py)
- [发现流程与组成基线](https://github.com/Lithium9767/MLP/blob/7bb22fcdb30363e7b0116dd0291835cb07f9b13b/scripts/m3_discovery.py)
- [实验台账](https://github.com/Lithium9767/MLP/blob/7bb22fcdb30363e7b0116dd0291835cb07f9b13b/experiments/registry.csv)
- [运行收据、统计表及图](https://github.com/Lithium9767/MLP/tree/7bb22fcdb30363e7b0116dd0291835cb07f9b13b/reports/M3/runs)
- [候选与长度混杂诊断](https://github.com/Lithium9767/MLP/blob/7bb22fcdb30363e7b0116dd0291835cb07f9b13b/reports/M3/candidate_freeze_decision.md)
- [C 交接包技术检查](https://github.com/Lithium9767/MLP/blob/7bb22fcdb30363e7b0116dd0291835cb07f9b13b/reports/M3/c_handoff_technical_check.md)
- [E 修订包技术检查](https://github.com/Lithium9767/MLP/blob/7bb22fcdb30363e7b0116dd0291835cb07f9b13b/reports/M3/e_handoff_technical_check.md)

代码、原始实验和成员贡献保留原 PR 的历史。本次仅发布统稿材料，不代替 B/C/E 提交其原始工作，也不合并未审查分支。原始序列、论文、交接 ZIP、权重、embedding 和大数组不包含在本目录。

## 结果与边界

- 已有 1202 条主分析 discovery 序列、24581 个窗口、12 组参数实验与氨基酸组成基线。
- 原始 CSV 的组成基线簇数为 **74**，噪声比例 68.33%；ESM-2 主配置 4 簇、噪声 30.68%。聚类指标不是功能准确率。
- 簇 3 的长度诊断 AUC 为 0.9688，候选尚未冻结，519 条主分析 validation 尚未运行；M3 尚未最终验收。
- 对“研究设计与证据链”类页面，推荐写“模型与 Baseline 已运行 / 初步候选已获得 / 混杂诊断已完成”，后续写“控制混杂与确定候选 / 独立验证 / 多源证据分析”。尚无新的控制实验时，不标为“混杂已控制”或“修订已完成”。

## 验证与复现说明

这次发布未修改模型或重新运行实验，因此不将历史测试描述为本次运行。已核对 PPT 的 ZIP/XML、内部关系目标、9 页显示与 3 页隐藏设置、文件哈希及差异；已检查近似版式预览。原生 PowerPoint 渲染在制作环境不可用，放映电脑仍需确认字体与换行。

PPT 的备注与讲稿列有来源；讲稿中 `MLP/` 表示项目根目录。报告中的运行 SHA 来自实验收据，本次汇报归档提交不作为运行 SHA。
