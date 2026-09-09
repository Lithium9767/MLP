# 验证记录

C 分支完整验证命令：

```powershell
python -m pytest tests -q -p no:cacheprovider
python -m bioinformatics.demo --out-dir outputs/c_demo_v1
```

2026-09-09：152 项测试通过，环境为 Windows、Python 3.14.3、pytest 9.1.1、matplotlib 3.10.8。示例目录须不存在，复跑请使用新路径。示例全部是合成计算材料，没有真实生物学结论。

| 测试文件 | 数量 | 主要验证 |
| --- | ---: | --- |
| `test_msa_analysis.py` | 23 | 已知熵值、缺口/未知残基、坐标恢复、输入核对、图与哈希、拒绝覆盖 |
| `test_a3m.py` | 31 | 插入/匹配状态、辅助记录、原始序列一致性、严格格式与输出 |
| `test_data_audit.py` | 33 | ID/序列/哈希、排除记录、标签缺失与冲突、划分覆盖和跨集泄漏 |
| `test_pfam.py` | 46 | TSV 与原始序列对应、区间合法性、重叠并集、未命中的不确定状态 |
| `test_regions.py` | 17 | 连续位置、缺失/低支持断点、显式阈值、输入校验和输出 |
| `test_workflow.py` | 2 | FASTA/A3M 完整流程、插入位 Pfam 连接、完整残基重建、哈希与禁止覆盖 |

另运行 `demo`，验证数据审计、FASTA/Pfam、A3M 三条路径并人工检查两张保守性 PNG 图。首次沙箱内 pytest 遇到 Windows 临时目录 ACL 限制，随后通过自动批准的正常沙箱外执行通过；未修改测试断言。

目前未覆盖真实数据验收、模型 token—残基映射、模型训练与归因统计，因为相应数据或实现未交付。B 分支测试不在 C 当前代码树内，本表不宣称运行了 B 全套测试。最终核心流程仍需非主责成员从空环境复现。
