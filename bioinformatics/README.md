# C：生物信息分析工具与复现

已完成数据交接审计、A3M 转换与核对、MSA 坐标和信息熵、保守区域提取、InterProScan/Pfam 结果解析，以及按完整序列残基坐标连接这些结果的流程。各工具均可独立运行；真实项目完成情况见 [C 任务清单](../reports/C_status.md)。

**PF01132 是延伸因子 P 的 OB 结构域；气囊蛋白家族 PF00741 也包含 GvpJ 等成员。** 家族命中不等于 GvpA 功能标签，正式对照对象仍须团队确认。见 [事实核查](../reports/M1/C_fact_check.md)和 [标签可行性表 v0.1](../reports/M1/C_label_feasibility_v0.1.md)。

## 一条命令验收

在仓库根目录运行，需 Python、matplotlib；测试另需 pytest。核心计算使用标准库，不需要 MMseqs2、比对软件或数据库服务。

```powershell
python -m bioinformatics.demo --out-dir outputs/c_demo_v1
python -m pytest tests -q -p no:cacheprovider
```

`demo` 完整运行无标签数据交接审计、FASTA/Pfam 注释连接、含插入的 A3M 分析。输出目录必须不存在；复跑请换路径。所有示例序列和 Pfam 区间均为人工构造，**不是 GvpA 数据、数据库扫描结果或生物学发现**。阈值 0.9、最短 2 列仅用于演示，不是正式评价规则。

| demo 输出 | 内容 |
| --- | --- |
| `data_audit/audit.json` | 序列/元数据一致性、来源缺失、标签/物种/长度分布、划分状态 |
| `fasta_pfam/` | 4 条合成序列的 MSA、区域、Pfam 和 69 个完整残基的连接表 |
| `a3m/` | 3 条合成序列的 A3M 转换、MSA、区域和 18 个完整残基的连接表 |
| `manifest.json` | 执行命令、Git 状态、产物 SHA-256 |

图表为 PNG/SVG。`FAILED.txt` 表示流程未完整成功，不应使用部分产物；成功流程生成 `manifest.json`。数据审计另以 `audit.json` 的状态与退出码判定。所有工具拒绝覆盖旧目录，不修改输入。

## 模块与输入

入口为 `python -m bioinformatics.<模块>`，各模块均支持 `--help`。

| 模块 | 输入与工作 |
| --- | --- |
| `data_audit` | 原始 FASTA、metadata、可选 split；审核 ID、序列/长度/哈希、证据缺失、相同序列或簇跨集及分布 |
| `a3m` | 已有 A3M、可选原始 FASTA；区分匹配位和插入，转换等长匹配位比对并保留完整坐标 |
| `msa_analysis` | 已有等长 FASTA 比对、可选原始 FASTA；逐列统计、位置映射、熵和保守性曲线 |
| `regions` | `column_scores.csv` 和显式阈值；提取满足支持条件的连续保守区间 |
| `pfam` | 原始 FASTA、已有 InterProScan 5 TSV、指定 accession 与来源版本；核对区间、覆盖和逐残基注释 |
| `workflow` | FASTA/A3M 比对，及可选完整 Pfam 输入；串联分析，输出按完整原始位置连接的注释表 |

数据审计采用当前 B 交接字段约定，详见 [B 实现复核](../reports/M1/C_B_review.md)，不替代团队冻结数据定义。未知标签始终保留为未知，格式通过不会自动批准监督训练。

## 完整分析入口

FASTA 合成示例：

```powershell
python -m bioinformatics.workflow --alignment bioinformatics/examples/synthetic_alignment.fasta --format fasta --original bioinformatics/examples/synthetic_original.fasta --out-dir outputs/c_fasta_v1 --run-id C-fasta-demo --dataset-version synthetic-v1 --purpose synthetic --alignment-tool hand-constructed --alignment-tool-version 1 --min-conservation 0.9 --min-length 2 --title "Synthetic fixture - not a biological result"
```

A3M 合成示例：

```powershell
python -m bioinformatics.workflow --alignment bioinformatics/examples/synthetic.a3m --format a3m --original bioinformatics/examples/synthetic_a3m_original.fasta --out-dir outputs/c_a3m_v1 --run-id C-a3m-demo --dataset-version synthetic-v1 --purpose synthetic --alignment-tool hand-constructed --alignment-tool-version 1 --min-conservation 0.9 --min-length 2 --title "Synthetic fixture - not a biological result"
```

连接 Pfam 时须提供原始 FASTA，并同时补充 `--interproscan-tsv`、`--pfam-accession`、`--interproscan-version`、`--pfam-version`、`--source-description`。完整合成参数可见 [demo.py](demo.py)。工具不执行数据库扫描；版本参数是运行者填写的来源记录，程序不能独立证明其真实性。

## 坐标与 A3M

所有生物学位置均为 **1-based，区间两端包含**。FASTA 标题首个空白前的内容为唯一 ID。给定原始 FASTA 时要求 ID 集合及每条完整序列精确一致，不默默截短或猜测 ID。

普通 FASTA 比对统一大小写，将 `.` 作为缺口。A3M 必须走 `a3m` 或 `workflow --format a3m`，不能直接转大写：大写字母和 `-` 为匹配状态，小写为插入，`.` 为插入区填充。程序排除并记录 `ss_`、`sa_`、`aa_` 辅助记录，核查匹配列数和字符。

**完整流程根目录的 `residue_annotations.csv` 是向 D/E 交接的坐标表：**

- 每个完整序列残基一行，以 `sequence_id + residue_position` 连接，并保留残基身份。
- `alignment_position` 是 FASTA 比对列或 A3M 匹配列；缺口没有原始残基行。
- A3M 插入保留完整原始位置，但匹配列、熵、保守性、区域归属留空。`insertion_anchor` 是插入前经过的匹配列数，`insertion_rank` 不代表跨序列同源关系。
- 插入残基仍可能位于 Pfam 的原始序列区间内，Pfam 按完整原始坐标连接。
- A3M 流程的 `msa/coordinate_map.csv` 使用去插入序列的残基索引，**不能作为完整原始坐标**；使用根目录连接表。`a3m/coordinate_map.csv` 同时保留两种索引。
- 未提供原始文件时，只能恢复输入比对中的序列，不能证明它与未截短的数据库序列一致；manifest 标记 `original_verified=false`。

接收 `B/J/O/U/X/Z` 并保留坐标，但不纳入 20 字母熵计算。拒绝非 ASCII、重复 ID、空/全缺口序列、非法残基、内部空白和不等长匹配比对；全缺口列可保留。

## 熵与保守区域定义

每列仅统计 20 种标准氨基酸，以等权序列频率计算：

`H = -sum(p[a] * log2(p[a]))`；`conservation = 1 - H / log2(20)`。

缺口和非标准残基不进入频率分母，但单独报告数量/比例。没有标准残基的列，熵和保守性留空。默认充分支持要求至少 2 个标准残基、占全部序列至少 50%，可用 `--min-canonical`、`--min-fraction` 调整。这不是显著性检验。

区域提取要求显式 `--min-conservation`、`--min-length`，只连接连续、充分支持且达到阈值的列，不跨缺失或低支持列；输出 `regions/regions.csv` 与逐列归属表。这是描述性区间，不能仅凭该表宣称功能关键区域。

未做物种/同源簇加权，近缘序列可能主导结果。跨物种分析还需物种信息、冗余控制和组间复核。正式区域选择与阈值应遵守划分及评价规则，不能利用测试集挑选最有利参数。

## Pfam 结果与缺失值

解析标准 InterProScan 5 TSV，核查 ID、MD5、完整序列长度、区间、日期/状态与 accession；其他分析行也需通过结构/序列一致性核查。未知 ID、错位或越界直接失败。

输出 `hits.csv`、`residue_annotations.csv`、`per_sequence_coverage.csv`、审计及来源清单。覆盖长度按同一序列命中区间并集计算，不重复累计重叠部分。区间起止、相对位置和逐残基覆盖可用于后续位置分布分析。

TSV 通常只列命中，不能用缺行证明扫描完成且不存在家族。无目标命中的序列保留不确定状态和空覆盖值；已有目标命中的序列才区分命中区间内外。`all_sequences_reported_hit_fraction` 只是输入中有报告命中的比例，不是扫描阴性率。**Pfam 未注释不等于没有功能。**

## 验收边界

152 项测试覆盖数学边界、输入损坏、A3M 恢复、数据交接、Pfam 坐标/覆盖和端到端连接。合成示例已运行并检查图表，见 [验收与剩余依赖](../reports/C_status.md)。

未执行真实 GvpA 的比对、数据库扫描、模型训练、归因检验或结构映射。`--purpose` 仅接受 `synthetic` / `exploratory`；正式项目还需实际输入、版本、台账及非作者复核。
