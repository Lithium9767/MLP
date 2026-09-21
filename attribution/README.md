# 残基证据与一致性

主责：E；复核：C/D。

attention只作为sequence-dependency evidence，不能单独称为关键位点归因。主线至少实现一种残基级直接扰动：

- mask；
- alanine substitution；
- 短窗口遮挡或删除。

每种方法必须说明标量目标、参考序列、坐标映射和聚合规则。输出保留原始残基位置、PF00741 HMM列、文献编号和PDB编号。比较扰动、attention、保守性、实验突变和结构证据，并加入等长度随机区域/随机位点基线。

候选位点只在validation集和敏感性队列中按冻结规则评估。无湿实验时输出“候选关键残基”和证据等级。
