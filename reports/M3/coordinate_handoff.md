# M3 C：窗口、PF00741 与 7R1C 坐标交接

状态：**C 的窗口接口和冻结来源已就绪，正在交付完整 discovery primary 窗口网格的坐标；候选冻结后生物学核对及 validation 映射待 A 冻结。** 本文件是 C 的交接，不代替 A 对 M3 的验收或组员的人工独立 Review。M2 已由 A 验收，PR #19 已合并；旧运行 `dirty=true` 的历史记录不改写，本交接复用干净 run 002。

## 可取得的冻结来源

- [PR #19](https://github.com/Lithium9767/MLP/pull/19) 的 C 原始扫描包位于 [公开 Release](https://github.com/Lithium9767/MLP/releases/download/c-m2-scan-002-handoff/MLP_C_M2_scan_002_handoff.zip)；完整 ZIP SHA-256 是 `7c65d1d4ca33918e14aed1d01e0b45b002a90e6e47a55562064ffb382698691e`，逐文件清单在 [M2 交接清单](../M2/c_scan_002_handoff.json)。自动化下载与 13 个 ZIP 成员、内层 B 输入三文件的校验见 [核验收据](pr19_handoff_verification.json)。该收据标注 `assistant_automated`，**不是非作者组员签字**；M2 验收决定另见 [Issue #15](https://github.com/Lithium9767/MLP/issues/15)。
- 本机 B 冻结输入实际位于 `data/processed/gvpa_v1_reproduction/`；不使用旧文档中的 `gvpa_v1/` 路径。M2 scan 002 的 `per_sequence_status.csv` 和 `hmm_coordinate_map.csv` 位于被 Git 忽略的 `data/processed/m2_c/pf00741_scan_002/`。逐位 7R1C 表位于 `data/processed/m2_c/7r1c_002/7r1c_residue_hmm_map.csv`。从 Release 获取后必须按清单复核 SHA-256，不在 Git 提交原始序列或大坐标表。
- 五份被引用的 M2 小型 JSON 收据原先记录 CRLF 字节哈希，但 Git 的自动换行可能让检出副本变成 LF。此分支通过 `.gitattributes -text` 将其**原有工作副本字节**按记录哈希存入 Git；只修复跨平台字节一致性，不改 JSON 字段、运行结果或历史审核状态。
- `python -m scripts.m3_verify_pr19_handoff --out <新收据路径>` 可重新下载核验；不要把自动核验写成 M2 正式验收。引用的冻结数据版本为 `gvpa-recognition-c7f6f005d717`，划分版本为 `homology-8b9005e2d9-s42`。

## 新收到的 B M3 数据包

B 交付的 `m3_b_handoff.zip` 已作本地只读核验；2026-10-05 用户明确提供微信收到的文件并确认它是 B 的交付，实际 ZIP 与先前收到的包字节一致。完整 ZIP SHA-256 为 `d4d5a53c78c10eb4bfc61008ec17eae225b9cab6305c69b7ed194de7d63307a6`，内层 `m3_b_handoff_bundle.zip` 为 `5bb7649e1c8f6de6830e71dc7d024157e7055257cf3add28f22d78f067ada18c`。外层三份 manifest 与内层副本逐字节相同。自动核验细节和每份输入哈希见 [B 包核验收据](b_handoff_verification.json)；复核入口为 `python -m scripts.m3_verify_b_handoff --zip <B 包路径> --out <新收据路径>`。用户确认交付方不等于确认 B 独立重跑的命令、环境或运行 SHA；未知来源不补造。

三份清单包含 [discovery primary](../../data/splits/m3/discovery_primary_manifest.csv) 1,202、[validation primary](../../data/splits/m3/validation_primary_manifest.csv) 519、[sensitivity](../../data/splits/m3/sensitivity_manifest.csv) 355，合计 2,076 个唯一 ID；每条记录的序列哈希、长度、队列、划分及同源簇与 M2 冻结数据一致，478 个同源簇不跨 discovery/validation。本分支只收录这三份**不含原始序列**的字节级副本，供 D 按收据 SHA 取得清单。内层 `metadata.csv`、`sequences_for_clustering.fasta`、`split_manifest.csv` 的 SHA-256 与上面 C 已核验的 M2 冻结输入完全相同，因此 **C 坐标索引无需迁移或重扫**。这属于 B→C 输入兼容性检查，不是用 validation 挑选窗口、调整参数或验证候选。

新包没有 B 自带的运行收据、完整性检查脚本或固定团队共享 URI；目前也没有其他组员实际下载该新包的核验记录。现在可用本分支的三份 manifest 配合 PR #19 Release 中同 SHA 的 M2 原始三文件重建输入，但这不能冒充 B 原包的独立下载签收或 A 的正式验收。包内有 validation 全序列；正式 validation 窗口核对仍须等待 A 冻结候选和评价规则。

## 已核实坐标与覆盖

`PF00741.24` 有 39 个 HMM match state。天然序列 `raw_position`、窗口边界采用 **1-based 闭区间**；HMM match state 是另一套同源坐标。M2 原始扫描有 2,076 条 QC 合格序列、80,219 行接受域逐残基映射；primary 队列 discovery 1,202、validation 519。入口同时锁定 B 冻结文件、三份 M2 收据的字节 SHA、收据中的 HMM/PDB 来源 SHA 及逐表 SHA；不能只同步修改某个本地收据来替换来源。完整输入、原始表、结构表哈希与数量的启动前复核见 [coordinate_summary.json](../../results/M3/coordinate_summary.json)。这些是 M2 坐标完整性事实，不是 M3 的发现或验证结果。

7R1C 使用 PDB **作者链 N**（RCSB label 链 A）。88 个 SEQRES 提交位点中 65 个有 CA 实验建模，23 个未建模；未建模位点为提交位置 1 和 67–88。HMM 状态 1–39 在该参考中对应提交位点／作者编号 11–49，均已建模。50–66 虽建模却不在该局部 HMM 对齐内。天然序列通过 HMM 状态连接 7R1C 只表示**同源位置**，不能把两个残基称为相同残基，也不能由 PF00741 命中推断功能。未建模位置保留 `modeled=false`、PDB 作者编号为空，不补造结构坐标。

## 给 D/E 的窗口接口

使用 [m3_window_coordinates.py](../../bioinformatics/m3_window_coordinates.py) 的 `load_coordinate_index(records, status_path, coordinate_path, structure_path)`；`records` 应由 M2 `validate_b_handoff(...)` 提供，先核对原始 SHA 再加载。随后调用 `index.map_window(internal_id, start, end, min_match_fraction=0.8)`。`start/end` 是天然序列 **1-based 且两端包含**；如果 D 的 Python 窗口为 `[start0, end0)`，必须显式转换为 `start=start0+1, end=end0`，并检查所得长度，不能按列名猜坐标系。

批量入口 `python -m scripts.m3_coordinate_handoff` 在不提供窗口时只校验 M2 来源并生成小型技术摘要，不进行候选分析。拿到 D 的窗口文件后，新建输出目录运行：

```powershell
python -m scripts.m3_coordinate_handoff `
  --windows data/processed/m3_d/discovery_windows.csv `
  --split discovery --cohort primary `
  --out-dir data/processed/m3_c/discovery_window_map_001
```

窗口 CSV 必须有 `internal_id,start_1based,end_1based_inclusive`，可有唯一 `window_id`。入口拒绝未知 ID、重复窗口、越界、队列/划分错用，输出 `window_coordinate_map.csv`、`residue_coordinate_map.csv` 和含源/输出哈希的 `mapping_receipt.json`；大表留在共享存储。`--cohort sensitivity` 用于非 primary 的单独交接，不将 partial／类型冲突混入主分析。对 validation 窗口还必须提供 `--split validation --validation-freeze <A 已冻结的文件>`；程序核对下面的 Markdown front matter 并记录文件 SHA。**示例是格式说明，不是实际冻结：**

```text
---
status: frozen
frozen_at_utc: 2026-09-29T10:00:00Z
candidate_manifest_sha256: <64 位小写 SHA-256>
evaluation_rules_sha256: <64 位小写 SHA-256>
---
```

实际值须由 A 在打开 validation 前写入并提交。程序只检查格式与文件哈希，冻结时间是否真的早于验证、候选和规则的内容是否充分仍由 A 审核。目前没有实际冻结文件，因此不能运行该阶段。

已用 discovery primary 的 `GVPA_000005` 两个预定位置 30 aa 窗口做**纯接口冒烟测试**：生成 2 行窗口表、60 行逐残基表与输出哈希，分别检查域外残基留空及域内状态连接。最终代码对应的文件在被忽略的 `data/processed/m3_c/smoke_discovery_map_003/`，不作为候选发现、参数选择或 validation 证据。

每个残基标记为 `match`、`insertion`、`unaligned` 或 `ambiguous_multi_domain`。插入没有 HMM state；域外/未对齐不强制填状态；重叠多域保留候选坐标但不给单一答案。窗口输出列出映射比例、HMM state 集合、仅在连续对齐内部可判定的缺失 state，以及 `mapped`／`low_match_coverage`／`unmapped`／`multi_domain_window`／`ambiguous_multi_domain`。默认 0.8 仅是**技术覆盖标记**，不是已冻结的发现阈值；A 应在查看 validation 前冻结候选与评价规则。7R1C 连接另标 `reference_modeled`、`reference_unmodeled`、`not_mapped_in_7r1c` 或 `no_match_state`。

## 尚待后续交接

C 可先交付完整 discovery 网格，D 无需为坐标接口重新运行 ESM-2。实际簇标签、入选窗口和候选定义仍由 D 提供；A 的具体候选及评价规则冻结记录尚未取得，所以不对 validation 作候选核对。冻结后再按同一接口处理 validation 与生物学证据。任何保守性用于发现时仅从 discovery 计算，不能偷看 validation。结构、PF00741 及同源位置都不是独立功能标签。

## 完整 discovery 网格交接

入口为 `python -m scripts.m3_discovery_coordinate_grid --window-length 30 --step 5 --out-dir <新目录>`；25/35 aa 敏感性网格用同一入口单独输出。枚举规则与 D 的 PR #22 提交 `79009c83c191a77144181346693b5c2f82e514fc` 中 `features.esm2_embed.iter_windows` 一致：`range(0, length-width+1, 5)`，转换为 1-based 闭区间，不追加不满长的末端窗口。仅处理 1202 条 discovery primary；不进行候选筛选、保守性计算或 validation 窗口映射。

输出 `windows.csv`、`window_coordinate_map.csv`、`residue_coordinate_map.csv` 与 `mapping_receipt.json`。公开逐残基表省略天然/结构参照的氨基酸字母，保留位置、状态、覆盖和编号；有序列需求时复用已核验的 B 输入。D 的 `raw_start/raw_end` 与 C 的 `start_1based/end_1based_inclusive` 均为 1-based 闭区间，以 `(internal_id, start, end)` 连接，要求唯一匹配，并核对序列 SHA。不要按表格行号或聚类编号连接。

0.8 覆盖只标记低覆盖，所有窗口均保留；空 HMM state 不补成 state 1。未建模的全部 23 个 7R1C 提交位点保留在原始结构表和摘要中，即使不属于 39 个对齐 HMM 状态，也不丢弃。D/E 可使用窗口/残基表做映射与覆盖抽查；候选选择和冻结后的生物学结论仍由对应负责人记录。
