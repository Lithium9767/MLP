# 验证计划
已实现 C 的独立 MSA 分析测试：`python -m pytest tests/test_msa_analysis.py -q`。
覆盖已知信息熵、缺口/未知残基、坐标往返、原始序列核对、错误 FASTA、CLI 端到端输出、文件哈希及拒绝覆盖。数据全部为合成计算示例，不是生物学实验结果。
TODO：B/D 实现 ID 与 similarity_cluster 不跨集合、元数据字段和标签审计；C/E 实现残基/MSA/token 坐标映射校验；D/E 实现小样本训练评估冒烟检查及指标输出一致性检查。
真实测试添加后由对应 PR 给出运行命令与结果。最终由非主责成员从空环境复跑核心流程。
