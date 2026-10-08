# M3 A/D 阶段报告

执行：Codex 工程助手；项目负责人 A 暂兼 D。M2 已验收。当前已完成 A 的工程统筹准备和 D 的 discovery 计算，不等于 M3 全阶段完成或真实成员 Review。

## 实测交付

| 项目 | 实际结果 |
| --- | --- |
| 输入 | discovery primary 1202 条；固定 metadata、split、HMM 坐标哈希 |
| 模型 | ESM-2 35M，revision 6fbf070e65b0b7291e7bbcd451118c216cff79d8；第12层、480维、CPU |
| token/有限值/重复 | 逐序列 token 对齐及有限值检查通过；首批重复推理最大差异0 |
| 缓存 | 24条回读实际命中缓存；文件SHA与首次输出一致 |
| 全长表示 | 1202条残基表示与 mean pooling；164.04秒 |
| 主窗口 | 30 aa、步长5，24581个窗口；显式 raw/HMM 坐标 |
| 参数矩阵 | 25/30/35窗口 ×20/50 PCA ×50/100 min_cluster_size，共12组 |
| 主聚类 | 4簇；噪声30.68%；20PC累计解释方差89.45% |
| 简单基线 | 氨基酸组成聚类74簇，噪声68.33%；与ESM聚类ARI=-0.0687，不是功能准确率 |
| 随机窗口 | 按序列长度与相对位置匹配，50次；重叠窗口不可当独立实验样本 |
| 队列敏感性 | discovery primary1202、partial7、partial+conflict17、conflict227；固定主模型迁移，不重拟合 |
| validation | 未推理、未拟合、未调参；519 primary序列保持锁定 |

## 发现与限制

1. 不能把4个簇都叫功能区域。cluster0完全无HMM映射，原始表中state1是全零argmax占位；代码已修复，原始运行表保持原样，candidate_review.csv明确更正解释。cluster2很弥散，没有多数窗口共同支持的连续区段。
2. cluster1的HMM跨度16–39，但集中度不高于匹配随机窗口的95%分位；cluster3跨度14–39，265条序列/2216窗口，集中度0.6394，高于随机诊断95%分位0.5991。后者只是探索性支持，不能宣称功能证明，也不能直接等同于参考论文的原始序列25–55。
3. 参数稳定性不统一：主配置改min_cluster_size的ARI约0.9975，改PCA维度约0.7791；窗口25/PCA20约0.1600，窗口35/PCA20约0.5170。跨窗口ARI比较相同起点，不是完全相同片段；包含noise，不能替代区域级稳健性。
4. 长度AMI0.2900，位置AMI0.1301，物种AMI0.1091，C端延伸AMI0.1509；仍未排除混杂。AMI为诊断性关联，不是因果控制。
5. 固定主模型迁移时，type-conflict噪声比例约61.51%，primary约30.69%，partial约45.04%。不能未经解释把敏感性队列混入主分析。训练标签与approximate_predict会有一个边界窗口差异，收据保留其ARI，不冒称字节等价。
6. ESM初始化日志中的pooler为随机初始化；本实验只用预训练last_hidden_state和自行mean pooling，未使用pooler或分类预测。重复检查针对首批，未声称每条序列都做双次推理。

## A 已整理的材料

- docs/plans/2026-09-29-m3-a-d-execution.md 与 configs/m3_discovery.json 固定范围、参数、发现/验证边界。
- experiments/registry.csv登记6次真实运行；收据全部绑定实际干净运行SHA，归档提交不冒充运行提交。
- 本目录 runs/ 中保存实际收据、小型统计表与图；大数组、全窗口表、模型和joblib留在 data/processed/m3/。
- reports/M3/ad_evidence_matrix.csv 列出各项交付、证据和待完成事项。

## 仍未完成

- 教师对路线的真实确认、非作者成员审查、M3最终验收不得由助手代签。
- 候选冻结提案仍为 proposed_not_A_approved。应先审阅局部稳健性/长度混杂，再由A对具体簇和参数作决定；不能在看过validation后回头选择候选。
- scripts/m3_validate_frozen.py 已实现批准锁、519条验证推理与固定PCA/HDBSCAN迁移，未实际运行。它在读取validation序列前检查真实批准记录与提案/收据/模型哈希；不得伪造批准记录来解锁。
- B/C/E的进一步坐标、生物学证据、残基扰动与最终解释不是本轮的已完成工作。

## 下一步

审查本PR与候选表；优先补强cluster3的区域级稳健性及长度混杂解释。A如果决定将当前主配置作为探索性验证方案，须明确选定簇并批准其固定提案，之后才运行validation，不重拟合模型或挑选验证集参数。
