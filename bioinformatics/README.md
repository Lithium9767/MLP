# C：已完成的独立工作

本次完成两项可独立验收的成果：数据库与背景事实核查、已有蛋白质 MSA 的描述性分析工具。工具不依赖 B 的数据采集程序、D 的模型或 E 的归因实现。

**编号核查：仓库原文中的 PF01132 是延伸因子 P 的 OB 结构域，不是 GvpA 家族。PF00741 是 Gas vesicle protein 家族，但也包含 GvpJ 等成员，不能把家族命中直接当作 GvpA 功能标签。** 证据及项目影响见 [事实核查报告](../reports/M1/C_fact_check.md)。本分支不替团队冻结新的数据定义。

## 输入与输出

运行入口：`python -m bioinformatics.msa_analysis`。输入为**已有、等长的蛋白质 FASTA 比对**，不是未比对的序列。可选提供对应原始 FASTA，程序会核对 ID 集合，并逐条检查去掉比对缺口后是否准确还原原始序列。

FASTA 标题首个空白前的字符串作为 ID；ID 不得重复。大小写统一为大写，`.` 转为缺口 `-`。接受 20 种标准氨基酸及 `B/J/O/U/X/Z`；后六种仅参与坐标映射和非标准残基计数，不进入本工具的 20 字母熵估计。拒绝空序列、全缺口序列、不等长比对、停止符、数字及内部空白。全缺口的**列**允许保留。

| 文件 | 内容 |
| --- | --- |
| `coordinate_map.csv` | 每条序列、每个比对位置一行；含 `sequence_id`、`msa_position`、`residue_position`、`residue` |
| `column_scores.csv` | 各列残基支持数、缺口比例、标准残基比例、众数、信息熵、保守性及支持标记 |
| `conservation.png` / `.svg` | 信息熵、保守性与覆盖比例三联图；低支持位置用橙色叉号标出 |
| `manifest.json` | 输入/输出 SHA-256、运行编号、参数、来源工具与版本、Python/绘图库版本、Git 状态及实现文件哈希 |

坐标均为 **1-based**；缺口没有原始残基坐标，CSV 留空。未知残基依然占一个原始序列位置。并列众数全部输出并按字母排序，例如 `AC`。空值不会被写成零。

## 计算定义与适用范围

每列只对 20 种标准氨基酸计数，使用等权序列频率：

`H = -sum(p[a] * log2(p[a]))`；`conservation = 1 - H / log2(20)`。

缺口和非标准残基不进入频率分母，但其数量和比例单独报告。没有标准氨基酸的列，信息熵和保守性都为空。只有一个标准残基的列可以得到数学上的零熵，但默认不会被标成充分支持。

`sufficient_support` 默认要求至少 2 个标准残基、标准残基占全部序列至少 50%。这只是工具的描述性支持标记，可通过参数设置；它不是统计显著性、保守区域判定或团队已冻结的评价规则。没有按物种/同源簇加权，近缘重复序列仍可能主导结果。跨物种保守性及功能结论需要额外分析。

## 完整示例

在仓库根目录运行。核心计算仅用 Python 标准库，绘图需要 matplotlib，测试需要 pytest（已包含于根目录依赖清单）。独立使用时可执行：

```powershell
python -m pip install matplotlib pytest
python -m bioinformatics.msa_analysis --alignment bioinformatics/examples/synthetic_alignment.fasta --original bioinformatics/examples/synthetic_original.fasta --out-dir outputs/c_synthetic_v1 --run-id c-synthetic-v1 --dataset-version synthetic-v1 --purpose synthetic --alignment-tool hand-constructed --alignment-tool-version 1 --title "Synthetic fixture - not a biological result"
python -m pytest tests/test_msa_analysis.py -q
```

四条示例序列是人工构造的计算测试材料，**不是 GvpA 序列、真实训练数据或生物学实验结果**。示例含全缺口列、只有一个有效残基的列和未知残基。首次运行后再次执行同一输出路径会被拒绝；复跑请换用 `outputs/c_synthetic_v2` 等新目录。

工具不修改输入，也不覆盖输出目录。失败目录中的 `FAILED.txt` 表示产物不完整，不应使用；成功标志为 `manifest.json`。`--purpose` 仅接受 `synthetic` 或 `exploratory`，不声称替代正式实验台账或项目冻结流程。

## 验收边界

本工具从已有 MSA 输入到五个输出文件的流程已完整实现并测试；没有执行真正的多序列比对，也没有实现/声称完成 Pfam 扫描、正式保守区域发现、结构映射或模型归因验证。它可单独使用，后续真实项目可按已验证的输入格式调用。
