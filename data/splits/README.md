# 发现/验证划分规范

- 先按序列相似性建立 `homology_cluster`，同簇不得跨discovery和validation。
- discovery用于区域探索、聚类参数和阈值选择；validation只执行冻结规则。
- 16–46、25–55若来自既有论文或recognition工作，登记为预注册候选，不再用同一数据宣称重新发现。
- partial、Gvp类型冲突和主要物种组作为敏感性队列分别报告。
- 记录工具、版本、identity/coverage阈值、随机种子、数据版本和输入哈希。
- 保存cluster映射并验证ID、完全重复序列和同源簇均无交叉。
- 更改规则时创建新版本，不覆盖旧划分。

当前尚未提交正式同源划分文件，由B主责，C/D复核。
