# 数据划分规范
- 同一 similarity_cluster 不得跨 train、val、test。
- 测试集不得用于调参、选择归因阈值或反复筛选模型。
- 划分文件必须版本化，例如 split-v1/train_ids.txt；示例名称不代表已有数据。
- 旧划分不能直接覆盖；更改产生新版本并记录原因。
- 必须记录工具版本、相似度阈值、随机种子和划分比例。
- 保存聚类映射、数据版本及标签/物种分布统计，验证 ID 和簇均无交叉。

候选实现为 `preprocessing/mmseqs_cluster_split.py`。默认候选参数为 80% identity、80% coverage、
80/10/10 和 seed 42；参数未经过小组评审前不得称为冻结 split。当前环境安装 MMseqs2 时遇到网络不可达，
因此尚未生成真实划分文件。B 主责、D 复核。
