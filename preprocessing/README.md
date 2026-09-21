# 数据预处理

主责：B；复核：C/D。本目录负责数据审计、质量队列、同源聚类接口和发现/验证划分。

从老师提供的 `gv.zip` 生成本地候选数据：

```powershell
python preprocessing/prepare_gvpa_dataset.py `
  --input-zip C:\path\to\gv.zip `
  --output-dir data/processed/gvpa_v1
```

输出包括 `metadata.csv`、`members.csv`、同源聚类输入 `sequences_for_clustering.fasta`、`audit_summary.json` 和 `dataset_manifest.json`。`data/processed/` 不进入Git。`sequence_qc_eligible` 表示序列可进入同源聚类，`primary_analysis_eligible` 才表示无partial和类型冲突的主分析队列。

安装MMseqs2后，在主分析FASTA上生成同源簇：

```bash
mmseqs easy-cluster data/processed/gvpa_v1/sequences_for_clustering.fasta work/gvpa work/mmseqs_tmp \
  --min-seq-id 0.8 -c 0.8 --cov-mode 1
```

再生成按簇隔离的 discovery/validation 清单：

```powershell
python preprocessing/homology_split.py `
  --metadata data/processed/gvpa_v1/metadata.csv `
  --clusters-tsv work/gvpa_cluster.tsv `
  --output-dir data/processed/gvpa_v1/split
```

划分脚本要求每个主分析序列恰好出现在一个簇中；缺失、未知或重复归属都会终止运行。处理脚本不把未知值补成负例，也不把2078称为质量控制后的主分析规模。
