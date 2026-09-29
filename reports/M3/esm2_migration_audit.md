# ESM-2 最小迁移审计（M3草稿，不启动实验）

任务：[Issue #5](https://github.com/Lithium9767/MLP/issues/5)。原PR：[PR #4](https://github.com/Lithium9767/MLP/pull/4)，源提交 `e87550ca662154152e1cf21bebc0e9d6b4871b92`。新分支从main `5891c2f7a6eb79e48dd7b75175a3af8176e410a2` 建立，没有整体合并旧分支。

## 本PR范围

| 文件 | 来源/必要性 |
| --- | --- |
| features/esm2_embed.py | 从#4逐文件迁移，内容不改；延迟加载torch/transformers，读取冻结metadata/split，默认discovery primary |
| requirements-m3.txt | 沿用基础requirements中的numpy等；仅新增实际推理需要的torch、transformers，不迁入pyhmmer/UMAP/HDBSCAN |
| tests/test_esm2_embed.py | 提取#4的ESM输入/窗口两项测试；增加CLI默认锁住validation的回归测试，不加载权重 |
| reports/M3/esm2_migration_audit.md | 来源、排除范围、预检边界和启动门槛 |

没有修改M2代码、数据、receipt、坐标接口、README或实验台账。没有下载权重、生成embedding或跑聚类。

## 重叠与排除审计

| #4模块/内容 | 与main关系及处理 |
| --- | --- |
| bioinformatics/hmm_coordinates.py | 与m2_coordinates/m2_pf00741的HMM坐标逻辑重叠；旧方案阈值、多命中处理和不区分split的保守性受到Review质疑；排除 |
| bioinformatics/structure_mapping.py | 与M2结构映射重叠；旧逻辑只从已建模原子构造序列，可能丢失未建模残基；排除 |
| HMM/PDB下载代码 | main已有download_hmms/download_7r1c；旧版URL校验问题不随本PR引入；排除 |
| HMM、结构相关测试 | 依赖旧实现，排除；保留main原有34项测试 |
| 原始数据、论文、旧结果、聚类/归因等扩展 | 不属于最小ESM入口，全部排除 |

#4的fengbujue777 CHANGES_REQUESTED记录保留；本草稿不表示该Review已批准。#11关闭未合并，C技术内容已由#13归档；不重新引入旧实现。B/C旧分支和本地快照保留。

## 验证命令与适用范围

```powershell
python -m unittest discover -s tests -p test_*.py -v
python features/esm2_embed.py --help
git diff --check
git diff origin/main -- bioinformatics preprocessing scripts/validate_m2_release.py
```

测试只覆盖轻量输入/窗口/CLI隔离，不是ESM模型推理验收。模型权重加载、GPU、token逐位对齐、有限值与推理可重复性尚未验证。要求最后一条diff为空以确认未覆盖M2实现。

## M2关闭后的待办

1. 先满足[Issue #15](https://github.com/Lithium9767/MLP/issues/15)真实独立复核与A验收；本PR保持Draft，不能因已有代码即宣布M3启动。
2. 由D/E真实审查入口和接口；正式运行前固定模型revision、依赖版本、Git完整SHA、设备、种子和输入/输出哈希。目前receipt未完整记录这些字段。
3. 当前缓存路径只包含模型ID和序列hash，未包含模型revision/软件版本且会重写；正式运行需采用版本隔离目录并补齐来源。
4. iter_windows只给坐标；入口仅导出全长残基表示和mean pooling，尚未导出窗口池化/PCA/HDBSCAN，不宣称完成区域发现。
5. CLI validation默认锁定仅为误用保护；直接函数调用不能替代团队候选冻结制度。正式解锁必须有冻结候选与A确认，发现阶段不使用validation统计调参。
