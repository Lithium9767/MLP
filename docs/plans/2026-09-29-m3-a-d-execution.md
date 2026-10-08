# M3 A/D 执行方案

复用 Issue #5（ESM-2 入口）与 #6（M3 统筹）。M2 已验收，基线 main `355424f63d60376ca66c6cc33f92c94181c7995f`。本轮只推进 A 和暂由 A 兼任的 D，不冒充其他成员的工作或 Review。

## A 的具体交付

- 固定 configs/m3_discovery.json：输入/坐标哈希、模型 revision、随机种子、30/5 主窗口、25/35 敏感性、20/50 PCA、50/100 min_cluster_size。
- discovery primary 1202 条用于拟合；validation 519 条保持锁定；partial/type-conflict 不混入主分析。
- 参数主配置在运行前确定，所有敏感性配置都报告，不能事后挑选最美观的图。
- 实验编号：M3-D-SMOKE-001、M3-D-EMBED-001、M3-D-DISCOVERY-001。失败运行记录真实原因，运行前提交代码并保持工作区干净。
- 建立证据矩阵、台账和 reports/M3/m3_summary.md；结果进入独立 PR，M3 非作者 Review 规则保持有效。
- 教师确认、真实成员审查及最终候选冻结批准由真人完成，助手不能代签或虚构。

## D 的计算交付

1. 固定 revision 的真实 ESM-2 35M，CPU 冒烟测试 24 条；验证 token/残基、480 维、有限值及重复推理差异。
2. 全量 discovery primary 1202 条：全长残基表示和 mean pooling；缓存隔离模型 revision 与软件版本。
3. 从全长上下文残基表示生成滑窗均值（不重新编码孤立片段）；输出 raw/HMM 坐标与缺失映射比例。
4. PCA/HDBSCAN 在 discovery 窗口表示拟合；UMAP 仅展示固定子样本。
5. 报告参数 ARI、噪声比例、解释方差，长度/物种/位置/C 端延伸混杂；primary 没有 partial/type-conflict，不能声称检验了这些敏感性队列。
6. 与简单氨基酸组成聚类、长度/位置匹配随机窗口比较；随机诊断受重叠窗口影响，不报告为功能显著性。
7. 形成候选区域和冻结提案，明确为潜在区域。A 对具体提案批准后才开展 validation 的固定变换/分配。

## 边界与验收

此方案可完成 D 的 discovery 技术链；真实非作者审查、partial/type-conflict 单独队列分析、候选批准后的 validation 和教师确认仍需后续完成。任何未执行项保持未完成，不以开发代码冒充实测结果。

源自 #17 的表示入口与三项测试按文件复用；不覆盖 M2 的 HMM、结构与数据实现，不整体合并 #4。新运行入口为 scripts/m3_discovery.py，数据、权重、embedding、joblib 和大型候选表均在忽略目录。
