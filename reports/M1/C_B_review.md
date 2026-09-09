# C 对 B 数据交接实现的复核

复核日期：2026-09-09。B 代码版本：`82ea140d6fc811b5bb37167da927e221d369c609`。
状态：**完成代码边界审查与 C 侧交接审计工具；未完成真实数据审计。**

## 审查对象与证据边界

- [B 的候选数据处理脚本](https://github.com/Lithium9767/MLP/blob/82ea140d6fc811b5bb37167da927e221d369c609/preprocessing/prepare_gvpa_dataset.py)。
- [B 的分簇划分脚本](https://github.com/Lithium9767/MLP/blob/82ea140d6fc811b5bb37167da927e221d369c609/preprocessing/mmseqs_cluster_split.py)。
- [B 的教师数据初审报告](https://github.com/Lithium9767/MLP/blob/82ea140d6fc811b5bb37167da927e221d369c609/reports/M1/teacher_data_audit.md)。

B 报告记载 856 条候选 GvpA 序列、长度 60–100 aa、标签为空。这是 **B 报告的结果**，
不能当作 C 已独立复现的计数。当前仓库没有该批原始 FASTA、A3M、逐序列 metadata 或正式 split，
因此 C 无法核查其完整性、物种分布、实际标签证据及真实泄漏情况。

本次问题全部由**人工构造的边界数据**复现，不证明真实 856 条序列存在这些错误。
只读审查 B 代码，没有合并 B 分支，也没有改动 B 的实现。

## 可复现的问题

| 问题 | 人工输入与实际行为 | 对交接的影响 | 修复建议，供 B 处理 |
| --- | --- | --- | --- |
| 有效输出中 ID 不唯一 | 两条记录同为 `WP_1.1`，序列分别为 `ACD`、`EFG`，均被保留；FASTA 和 metadata 出现同名有效记录 | 以 ID 连接比对、标签、划分或归因结果时产生歧义或静默覆盖 | 检查 accession/ID 与完整序列的对应；冲突进入审计并明确排除或解决，不凭行号默默改成新样本 |
| 无效记录占用了去重代表 | 首条记录没有 accession，序列 `ACD`；第二条合法 `WP_2.1` 同为 `ACD`。前者被排除，但先写入去重表，导致后者也以重复理由排除，最终有效记录为 0 | 数据保留结果受无效记录顺序影响，合法样本可能丢失 | 将合法性检查与有效代表选择分开；去重代表必须是被保留的有效记录 |
| 同一 member 可出现在两个簇 | `r1 → {r1, a}`、`r2 → {r2, a}` 被接受；按 0.5/0.5/0 划分时 `a` 可同时写入 `train` 与 `val` | 异常聚类输入可能产生跨集合重复与泄漏 | 读入时检查每个 member 唯一归属；写出前检查 ID 唯一性、完整覆盖、簇与相同序列是否跨集 |

第三例直接调用 `assign_clusters()` 验证内部函数边界；0.5/0.5/0 是用于复现的函数参数，
不是推荐的正式划分比例，也不意味着 CLI 接受零比例。正式工具参数与评估规则仍由相关成员确认。

以下 Python 片段从指定 Git 提交读入两个模块，在临时目录构造上述三例。需要本地已有该提交、Python 与 Git，
不需要教师数据、MMseqs2 或联网，也不写入 B 的源码文件：

```python
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

commit = "82ea140d6fc811b5bb37167da927e221d369c609"

def module_at_commit(path):
    code = subprocess.run(["git", "show", f"{commit}:{path}"],
                          check=True, capture_output=True).stdout.decode("utf-8")
    namespace = {"__name__": "review_only"}
    exec(compile(code, path, "exec"), namespace)
    return namespace

prepare = module_at_commit("preprocessing/prepare_gvpa_dataset.py")
split = module_at_commit("preprocessing/mmseqs_cluster_split.py")
with TemporaryDirectory() as directory:
    root = Path(directory)
    duplicate = root / "duplicate.fasta"
    duplicate.write_text(">WP_1.1\nACD\n>WP_1.1\nEFG\n", encoding="utf-8")
    rows, _ = prepare["build_rows"](duplicate)
    print("included IDs:", [r["sequence_id"] for r in rows if not r["exclusion_reason"]])

    invalid_first = root / "invalid_first.fasta"
    invalid_first.write_text(">\nACD\n>WP_2.1\nACD\n", encoding="utf-8")
    rows, _ = prepare["build_rows"](invalid_first)
    print("exclusion reasons:", [r["exclusion_reason"] for r in rows])

    clusters_file = root / "clusters.tsv"
    clusters_file.write_text("r1\tr1\nr1\ta\nr2\tr2\nr2\ta\n", encoding="utf-8")
    clusters = split["read_clusters"](clusters_file)
    assignments = split["assign_clusters"](clusters, (0.5, 0.5, 0.0), 42)
    split["write_split"](clusters, assignments, root / "split.csv")
    print((root / "split.csv").read_text(encoding="utf-8"))
```

## 已完成的 C 侧交接审计工具

实现：[bioinformatics/data_audit.py](../../bioinformatics/data_audit.py)。只用 Python 标准库。
它核查 B 提供的 `sequences.fasta`、`metadata.csv` 和可选 `split.csv`，输出带输入文件 SHA256、
执行参数、实现 SHA256、时间及 Python 版本的 `audit.json`。

审计规则：

1. `exclusion_reason` 为空的 metadata 行才是有效记录，必须与 FASTA 严格一一对应，ID 唯一且完整序列一致。
   被排除记录保留在审计统计中，但不要求在 FASTA 中存在；仅被排除的 ID 出现在 FASTA 中则报错。
2. 核对序列长度、有效残基标志与 SHA256。先拒绝非 ASCII 字符，再转大写计算哈希；允许 B 使用的标准/歧义氨基酸字母，
   拒绝缺口、终止符和序列内部空白。这是对 B 当前实现的检查约定，尚不代表 A/B/C 已冻结正式数据定义。
3. 统计 accession、database、source_url、retrieval_date、organism、taxonomy_id、lineage 缺失情况；
   `database` 标成 inferred 时单独提示。字段存在不能替代来源真实性与许可核查。
4. 统计有效样本的逐长度频数、物种文本频数和标签频数；未知标签单列计数，不变成 0、阴性或新类别。
   物种字段只是原始文本分组，不声称完成 taxonomy 校正。
5. 报告已填标签缺证据字段、多个标签定义、相同序列冲突标签等情况；不推断标签、不决定证据等级。
   同一序列的冲突标签也可能涉及不同实验背景，需要人工核对。
6. 若提供 split，核查 ID 唯一性与精确覆盖、仅接受 `train/val/test`、要求簇 ID、检查簇和完全相同序列跨集情况，
   并检查外部 split 与 metadata 已填字段是否冲突。未单独提供 split 时，metadata 中已有的划分字段也会被审计。
   工具不执行 MMseqs2，不证明不同簇之间没有高相似度序列。

全无标签且文件有效时返回成功，明确标记 `labels.status=unlabelled`、`supervised_training_approved=false`，
可用于序列的无监督描述。即使所有标签字段齐全，工具也不自动批准监督训练。标签/证据不全是准备状态，
不会伪装成 FASTA 格式错误。数据完整性或 split 错误会写出报告并以退出码 2 结束。

## 复现与验证

提交的示例是合成序列，**不是 GvpA 数据或生物学实验结果**：

```powershell
python -m bioinformatics.data_audit --fasta bioinformatics/examples/synthetic_original.fasta --metadata bioinformatics/examples/synthetic_metadata.csv --out-dir outputs/C_synthetic_handoff_audit
python -m pytest tests/test_data_audit.py -q
```

`--out-dir` 必须是尚不存在的新目录，重复运行应换目录名；工具拒绝覆盖已有结果。
真实交接时，将上述两个输入替换为 B 实际交付文件，并在具备正式划分时增加 `--split 路径`。

本次 **33 项测试通过**，覆盖 ID 冲突、排除规则、完整序列/长度/哈希一致性、非法残基与 Unicode、
CSV 结构错误、标签缺失与冲突、划分覆盖/跨集泄漏、输入文件缺失、CLI 错误报告和禁止覆盖。
首次常规测试因 Windows 临时目录 ACL 被拒绝；随后在自动批准的正常沙箱外运行中通过，未修改或绕过 pytest 断言。

## 仍不能完成的部分

- **真实数据的生物学复核、实际分布统计与真实泄漏审计**：需要 B 交付可读取的版本化 FASTA、metadata，
  正式划分可用时再交付 split 与聚类参数/版本。当前仅有报告文字，不能据此制造逐序列结果。
- **来源与许可的逐项确认**：需要 B 的原始查询/获取记录、日期、版本和可追溯来源；不能把 accession 推断当成完整来源证明。
- **监督标签可行性结论**：需要 A 组织确定任务目标、B 提供逐序列证据记录、C 核查生物学含义、D 复核建模可行性。
  本工具提供完整核查能力，但不替代这些证据与决定。

这些限制不影响已完成工具和人工边界复核的验收。真实分析前仍需取得实际输入并执行审计。
