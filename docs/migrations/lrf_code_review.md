# lrf代码迁移审查

来源提交：`d272ad7fe5fc1c73e158e0bfee8df702c94969ec`。

本次只迁移代码思想，不迁移原始序列、论文PDF、旧图表、旧聚类结果或PDB文件。原提交已通过本地归档分支 `archive/lrf-d272ad7` 和离线ZIP快照保留。

## 已迁移

| 来源 | 新模块 | 调整 |
| --- | --- | --- |
| `features/esm2_embed.py` | `features/esm2_embed.py` | 改用internal ID和冻结split；validation默认锁定；缓存绑定序列哈希与模型；不把attention当归因 |
| `bioinformatics/download_hmms.py` | `bioinformatics/download_hmms.py` | 默认只下载PF00741；增加响应校验、时间、URL和SHA-256收据 |
| `bioinformatics/hmm_conservation.py` | `bioinformatics/hmm_coordinates.py` | 输出逐残基原始坐标、HMM match state和插入标记；支持保守性统计；使用当前metadata字段 |
| `experiments/recluster_and_structure.py`中的PDB解析 | `bioinformatics/structure_mapping.py` | 拆除全量重聚类；只保留PDB链、二级结构和HMM坐标映射；PDB文件不入库 |

## 未迁移

- `data/GvpA_sequences.json`、`gv/**`、FASTA、MSA和论文PDF：属于原始/外部资料，不进入公共仓库。
- `results/**`旧图表与结论：没有绑定当前 `homology-8b9005e2d9-s42`，不能作为正式结果。
- `models/unsupervised_analysis.py`：依赖缺失时静默改用层次聚类，会改变预注册方法。
- `experiments/recluster_and_structure.py`的聚类部分：在全量数据上重新拟合并读取旧结果，不满足discovery/validation隔离。
- PF00741/PF01132自然序列二分类逻辑：家族命中不是功能差异标签。

## 使用边界

这些模块目前是经过单元测试的运行入口，不代表ESM-2、PF00741坐标或7R1C映射已经在新数据版本上完成。正式运行必须记录可选依赖版本、输入哈希、split版本、Git提交和失败记录。
