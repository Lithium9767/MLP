# 配置
baseline.yaml 是计划接口，不是已经实现的训练程序。
target_label、dataset_version 和 split_version 尚未冻结，正式训练前必须填写；未来入口应拒绝这些字段为 null 的正式运行。
当前特征为 AAC；kmer_k: 2 为后续 k-mer 实验预留，不表示当前同时使用 k-mer。
deterministic 等字段表达复现要求，尚无实现保证。指标的二分类/多分类适用性、AUC 汇总与缺失类别处理由 A/D/E 冻结，TODO。
配置变更通过 PR，正式实验保存实际使用的完整配置快照。
