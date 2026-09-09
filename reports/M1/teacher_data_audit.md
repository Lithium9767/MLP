# 教师提供数据初审（成员 B）

状态：候选数据，未冻结。审计对象为本地 `gv/gv/`；原始目录已由根 `.gitignore` 排除。

## 选择结论

根据实验要求，B 的主数据候选使用 `design/data/GvpA.fasta`，只负责 GvpA 序列、元数据、清洗、
去冗余和 similarity-aware split。`GvpA_msa.a3m` 交由 C 负责，不纳入 B 的处理结论。

`recognition` 是八种 GVP 亚型识别数据，与本项目要求的 GvpA 功能属性预测不同，仅作为补充来源审计。
`circuit` 面向基因线路与调控元件设计，其中存在 `ATG...TAA` 一类占位序列，不作为主数据。

## GvpA 主候选集审计

| 项目 | 结果 |
| --- | ---: |
| FASTA 原始记录 | 856 |
| 格式/残基检查后记录 | 856 |
| 唯一 accession | 856 |
| 唯一序列 | 856 |
| 缺失 organism | 0 |
| 长度范围 | 60–100 aa |

候选集可通过 `python preprocessing/prepare_gvpa_dataset.py` 复现。监督功能标签尚未提供，
因此 `label` 保持为空；当前数据只能作为无标签 GvpA 序列集，不能交给 D 正式训练监督模型。

## recognition 数据概况

| 标签 | 条数 | 最短 | 中位长度 | 最长 |
| --- | ---: | ---: | ---: | ---: |
| GvpA | 2,078 | 59 | 143 | 177 |
| GvpC | 233 | 18 | 212 | 520 |
| GvpG | 3,576 | 77 | 85 | 134 |
| GvpJ | 3,135 | 32 | 136 | 274 |
| GvpK | 377 | 26 | 106 | 184 |
| GvpN | 758 | 245 | 326 | 490 |
| GvpO | 1,802 | 28 | 110 | 287 |
| GvpP | 220 | 135 | 171.5 | 623 |

合计 12,179 条；当前 JSON 内序列哈希和 accession 均唯一，未发现跨标签的完全相同序列；
30 条记录缺少 organism，需补查或明确缺失原因。
最大类与最小类数量比约 16.25，后续不能只报告 accuracy，应使用 Macro-F1、每类召回率、
PR-AUC/MCC 等，并采用类权重或仅在训练集内进行采样。

## 三套 GvpA 的关系

- `design` 与 `circuit`：分别 856/862 条，完全相同序列交集 821 条，二者高度重合。
- `recognition` 与 `design`：交集 425 条。
- `recognition` 与 `circuit`：交集 440 条。

因此不能把三套来源直接合并后随机切分。若后续使用 `design` 的 MSA 或结果作外部验证，必须先按
序列哈希及相似度检查其与训练/验证/测试集的重合，并明确它是内部参考还是独立外部证据。

## 当前阻塞项

1. A/教师需提供或确认逐序列功能标签；B 负责合并和质量审计，但不负责定义生物学标签。
2. 老师资料只有处理后结果，仍需追溯原始数据库查询、检索日期、版本、许可及 60–100 aa 筛选依据。
3. 当前运行环境未安装 MMseqs2；通过 conda 安装时网络不可达。
4. 分簇与按簇划分脚本已经完成，待 MMseqs2 可用后运行并审计泄漏。
5. 在功能标签和 similarity-aware split 完成前，不向 D 交付正式监督训练集。

## 可复现入口

运行 `python preprocessing/prepare_gvpa_dataset.py`，会在被 Git 忽略的
`data/processed/gvpa_v0.1_candidate/` 生成候选 `metadata.csv`、统一 FASTA 和审计 JSON。
相似度聚类与分组划分入口为 `preprocessing/mmseqs_cluster_split.py`。
