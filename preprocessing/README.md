# 数据预处理

主责：B；复核：C。

## GvpA 主候选集

教师提供的完整资料本地放在仓库根目录 `gv/`，该目录已被 Git 忽略。根据实验要求，B 只整理
GvpA 序列数据；MSA、Pfam 和保守性分析由 C 负责。主输入为
`gv/gv/design/data/GvpA.fasta`。运行：

```bash
python preprocessing/prepare_gvpa_dataset.py
```

默认输出到被 Git 忽略的
`data/processed/gvpa_v0.1_candidate/`：

- `metadata.csv`：符合项目最小字段规范的逐序列元数据；
- `sequences.fasta`：通过检查和精确去重后、用于后续聚类的统一 FASTA；
- `audit_summary.json`：类别、长度、重复与冲突审计。

该输出只是候选集，不是冻结数据。老师资料没有逐序列监督标签，脚本会保留空 `label` 并标记
`pending_prediction_target_confirmation`，不会把“属于 GvpA”伪装成待预测的功能标签。

`preprocessing/prepare_teacher_data.py` 仅用于审计 `recognition` 中的八分类资料，不作为当前主任务入口。

## 相似度聚类与划分

安装 MMseqs2 后运行以下候选配置：

```bash
python preprocessing/mmseqs_cluster_split.py \
  --fasta data/processed/gvpa_v0.1_candidate/sequences.fasta \
  --work-dir data/cache/mmseqs-gvpa-0.8 \
  --output data/processed/gvpa_v0.1_candidate/split.csv \
  --min-seq-id 0.8 --coverage 0.8 --ratios 0.8,0.1,0.1 --seed 42
```

0.8 identity/coverage 是待评审的候选参数，不是已冻结结论。脚本保证同一 MMseqs2 cluster
整体进入一个集合。正式冻结前需记录 MMseqs2 版本、簇数量、最大簇和三个集合的样本统计。

## 下一步

1. A/教师提供或确认逐序列预测标签；B 将标签合并到 `metadata.csv` 并审计缺失率。
2. B 用 MMseqs2 建立 similarity cluster，记录阈值、覆盖率、版本和命令。
3. B 按 cluster 分组生成 train/validation/test；同一簇不得跨集合。
4. B 输出各集合的标签、物种和长度分布，再申请冻结 D1/D2。
