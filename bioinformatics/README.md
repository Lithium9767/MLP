# 生物信息分析

主责：C；复核：B/E。

下一阶段需要：

1. 独立复跑PF00741全量扫描并记录HMMER/InterProScan、Pfam版本、GA阈值、命令和输入哈希；
2. 建立原始残基、窗口、PF00741 HMM列、文献编号和PDB 7R1C残基的显式映射；
3. 在HMM/MSA列上计算覆盖、熵、保守性和partial/类型冲突敏感性；
4. 逐条审计实验突变的亲本、位置、替换、表型、宿主和来源；
5. 用PDB 7R1C作为主要结构证据，ESMFold仅为可选补充。

PF00741可用于身份、同源坐标和候选区域验证，不作为当前天然GvpA的监督标签。MAFFT未实际运行前不得写成已完成结果。

本地试扫描工具见 [PILOT_PFAM.md](PILOT_PFAM.md)。

已从归档分支审查并迁移三个入口：

- `download_hmms.py`：下载PF00741 HMM并记录URL、时间和哈希；
- `hmm_coordinates.py`：输出原始残基到HMM match state的逐残基映射、插入标记和保守性；
- `structure_mapping.py`：将PDB 7R1C指定链映射到PF00741状态，不在仓库内保存PDB文件。

```powershell
python bioinformatics/download_hmms.py --output-dir data/raw/pfam
python bioinformatics/hmm_coordinates.py `
  --metadata data/processed/gvpa_v1/metadata.csv `
  --hmm data/raw/pfam/PF00741.hmm `
  --output-dir data/processed/gvpa_v1/pf00741
python bioinformatics/structure_mapping.py `
  --pdb data/raw/structure/7R1C.pdb `
  --chain N `
  --hmm data/raw/pfam/PF00741.hmm `
  --output-dir data/processed/gvpa_v1/structure
```

这些入口尚未在当前数据版本上完成正式运行，运行后必须把receipt和小型汇总登记到实验台账。
