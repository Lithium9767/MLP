# M2 数据地基设计

## 目标

把老师提供的 `gv.zip` 中 recognition GvpA JSON 转换为可追溯的候选数据集，并为按同源簇隔离的 discovery/validation 划分建立固定接口。该阶段不构造功能标签，也不提前运行候选区域发现。

## 方案选择

考虑过三种推进方式：

1. 直接复用开发分支中的 ESM-2、UMAP 和 HDBSCAN 结果。速度快，但数据版本、同源隔离和候选冻结均不完整，容易发生循环验证。
2. 只整理方法文档。风险小，但不能给后续成员提供可运行输入。
3. 先建立数据审计与同源划分流水线，再接入表示学习。增加一个前置里程碑，但能使后续结果绑定同一输入、质量队列和 split。

采用方案3。它先解决新路线中最靠前、也最容易影响所有后续结论的数据边界问题。

## 数据流

```text
gv.zip
  └─ recognition GvpA JSON
       ├─ metadata.csv          代表序列、质量标记、分析队列
       ├─ members.csv           原始成员与代表序列的关系
       ├─ sequences_for_clustering.fasta 通过序列质量检查的聚类输入
       ├─ audit_summary.json    计数、阈值、输入与输出哈希
       └─ dataset_manifest.json 数据版本和生成规则

sequences_for_clustering.fasta
  └─ MMseqs2 easy-cluster（外部步骤）
       └─ *_cluster.tsv
            ├─ split_manifest.csv
            └─ split_summary.json
```

原始序列和逐条成员表保存在 `data/processed/`，受 `.gitignore` 保护。仓库只提交脚本、测试、配置和不含序列的小型审计摘要。

## 质量规则

- 规范序列为大写并删除空白；主分析只接受20种标准氨基酸。
- recognition JSON 的完全重复状态重新核验；重复记录保留审计信息，但不进入主分析 FASTA。
- `partial` 同时检查代表注释和所有冗余成员注释。
- 类型冲突检查所有成员的 `gvp_types`，任何非 `GvpA` 类型单独标记。
- 来源、物种和 accession 缺失分别记录，不把未知值推断成阴性或某个物种。
- 长度异常使用本次数据中标准序列的均值±3个总体标准差，只作为敏感性标记，不自动排除。
- 分析队列区分 `primary`、partial、类型冲突、二者并存和非标准残基排除组。

## 同源划分

MMseqs2 建议先以 `min-seq-id=0.8`、覆盖度 `0.8` 生成簇。划分脚本只接受完整且无重复归属的 cluster TSV；主分析中的每个 sequence ID 必须恰好出现一次。簇按大小优先、种子控制的稳定顺序分配至 discovery/validation，目标比例为70/30。同一簇绝不跨集合。

由于当前 Windows 环境未检测到 MMseqs2，代码只冻结命令参数和消费 cluster TSV 的逻辑。没有真实 cluster TSV 时，split 状态必须保持 pending。

## 验收

- 对真实 `gv.zip` 能生成2078条候选记录及确定的数据版本。
- 审计摘要能复现 partial、类型冲突、非标准残基等计数。
- 输入ZIP、JSON成员、主分析FASTA和元数据均记录SHA-256。
- 合成测试证明完全重复、质量队列、缺失簇、重复簇归属和同簇隔离规则有效。
- 所有既有 PF00741 试扫描测试继续通过。
