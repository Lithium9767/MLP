# 生物信息分析

主责：C；复核：B/E。

下一阶段需要：

1. 独立复跑PF00741全量扫描并记录HMMER/InterProScan、Pfam版本、GA阈值、命令和输入哈希；
2. 建立原始残基、窗口、PF00741 HMM列、文献编号和PDB 7R1C残基的显式映射；
3. 在HMM/MSA列上计算覆盖、熵、保守性和partial/类型冲突敏感性；
4. 逐条审计实验突变的亲本、位置、替换、表型、宿主和来源；
5. 用PDB 7R1C作为主要结构证据，ESMFold仅为可选补充。

PF00741可用于身份、同源坐标和候选区域验证，不作为当前天然GvpA的监督标签。MAFFT未实际运行前不得写成已完成结果。

本地试扫描工具见 [PILOT_PFAM.md](PILOT_PFAM.md)。

M2 的干净扫描 002 及公开原始交接包见 [M2 C 运行说明](M2_PF00741.md)。M3 启动前的窗口坐标接口已在 [m3_window_coordinates.py](m3_window_coordinates.py)；批量校验/映射入口为 `python -m scripts.m3_coordinate_handoff`，明确保留未对齐、插入、多域、低覆盖和 7R1C 未建模状态。使用方法、来源哈希和剩余的冻结后交接见 [M3 C 坐标交接](../reports/M3/coordinate_handoff.md)。这不是候选发现或 validation 结果。
