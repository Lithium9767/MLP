# 无监督表示与区域分析

主责：D；复核：E。

当前主线不训练天然GvpA功能分类器。计划流程：

1. 在discovery集生成全长和窗口ESM-2表示；
2. 在原始表示或PCA空间运行HDBSCAN；
3. UMAP只用于二维展示；
4. 将窗口簇映射回原始位置和PF00741 HMM列；
5. 冻结候选区域与参数；
6. 在同源隔离的validation集复现；
7. 检查窗口、PCA、HDBSCAN参数和随机种子敏感性；
8. 检查长度、物种、partial和C端延伸等混杂。

旧LR/SVM/RF配置仅作为教师要求监督任务时的备用方案，不按当前天然序列数据启动。

谱聚类 K=12 的独立比较入口为 `spectral_regions.py` 和
`scripts/m3_spectral.py`。它只接收已核验的窗口特征，分开实现谱图拟合与
冻结 PCA 空间的邻居标签迁移；不等同于已训练天然功能分类器。
输入、运行命令及解释边界见 [中文说明](../reports/M3/spectral_k12.md)。

当前固定K12已完成519条保留验证序列的迁移，具体结果见
[验证报告](../reports/M3/spectral_validation.md)。`spectral_validation.py`
只统计固定标签的区域谱、组外距离及PCA残差，不拟合或改变簇标签；
原始序列到窗口特征由 `features/spectral_windows.py` 提供，正式验证
由 `scripts/m3_spectral_validate.py` 绑定冻结协议后编排。
