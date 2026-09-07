# GvpA 功能预测与关键位点归因

## 项目简介
五人机器学习课程项目，仓库为 [Lithium9767/MLP](https://github.com/Lithium9767/MLP)。构建“预测 → 归因 → 生物学验证 → 候选区域发现”的证据链。

## 研究背景
本项目计划从 GvpA 氨基酸序列出发，探索功能属性预测与残基级解释，将模型结果与 Pfam01132、MSA 保守性及可选结构信息比较。背景文献及数据库证据由 C 在 M1 核查，TODO。

## 研究目标
建立可追溯标签与无明显泄漏的数据划分；完成传统基线及可归因模型；比较至少两种归因方法，并报告候选区域及限制。

> 第一优先级是确定具有独立证据来源的预测标签。需要避免使用 Pfam 生成监督标签后，又使用同一个 Pfam 标注证明归因正确，从而形成循环论证。

## 三个研究问题
- RQ1：模型关注的关键位点是否与 Pfam01132 和 MSA 高保守位点一致？
- RQ2：不同归因方法得到的关键位点是否一致？
- RQ3：是否存在模型共同关注但 Pfam 尚未标注的跨物种保守区域？

## 技术路线
1. 审计预测目标、独立标签来源与证据等级。
2. 获取序列与元数据，清洗、去冗余，冻结 similarity-aware split。
3. 比较 AAC/k-mer 等手工特征与 Logistic Regression、SVM、Random Forest。
4. TODO：确认环境后评估 ESM2 frozen embedding 和残基级轻量分类头。
5. TODO：以 IG、Saliency 等不同归因方法进行一致性分析；Attention 仅作补充。
6. 对照 Pfam/MSA，加入随机位点或置换基线；有余力时做结构映射。
7. 在模型质量达标后报告 RQ 证据，候选区域不等同于实验验证。

## 五人分工
| 成员 | 主角色 | 主要职责 | 互审对象 |
| --- | --- | --- | --- |
| A | 组长与实验统筹 | RQ 拆解、排期、实验矩阵、报告和答辩整合 | 复核模型结论是否对应 RQ |
| B | 数据工程负责人 | NCBI/UniProt 采集、清洗、去冗余、相似度划分 | 与 C 互审元数据 |
| C | 生物信息负责人 | MSA、信息熵、Pfam01132、保守区域和结构验证 | 复核 D/E 的生物学合理性 |
| D | 模型负责人 | 手工特征、传统 ML、ESM2、训练和评估 | 与 E 互审归因接口 |
| E | 归因与可视化负责人 | IG、Saliency、SHAP、一致性分析和图表 | 复核 D 的评估结果 |

姓名与 GitHub 账号映射：TODO，由成员自行确认。

## M1—M5 里程碑
| 阶段 | 建议时间 | 交付与验收 |
| --- | --- | --- |
| M1 | 第 1 周 | 预测目标、标签可行性、实验设计与仓库 |
| M2 | 第 2 周 | 数据卡、去冗余、固定划分、EDA 与 MSA/Pfam |
| M3 | 第 3 周 | 三种基线、统一指标与初步归因 |
| M4 | 第 4—5 周 | 主模型、多归因、一致性与 RQ 证据 |
| M5 | 第 6 周 | 报告、PPT、结果冻结和交叉复现 |

实际日期 TODO，按课程截止时间调整；GitHub Milestones/看板尚待创建。

## 仓库结构
```text
MLP/
├── README.md
├── CONTRIBUTING.md
├── requirements.txt
├── .gitignore
├── .github/
│   ├── ISSUE_TEMPLATE/
│   │   ├── task.md
│   │   ├── experiment.md
│   │   └── bug.md
│   └── pull_request_template.md
├── configs/
│   ├── README.md
│   └── baseline.yaml
├── data/
│   ├── README.md
│   ├── dataset_card.md
│   └── splits/README.md
├── preprocessing/README.md
├── bioinformatics/README.md
├── features/README.md
├── models/README.md
├── attribution/README.md
├── experiments/
│   ├── README.md
│   └── registry.csv
├── results/README.md
├── reports/
│   ├── README.md
│   ├── meetings/README.md
│   └── M1/problem_definition.md
└── tests/README.md
```

## 环境安装
建议使用独立 Python 虚拟环境，具体支持版本与依赖锁定待 D 验证。
```bash
python -m venv .venv
```
Windows PowerShell 激活：`.venv\Scripts\Activate.ps1`；Linux/macOS：`source .venv/bin/activate`。
```bash
python -m pip install -r requirements.txt
```
当前依赖未锁版本，尚未验证完整运行环境。torch、transformers、fair-esm 等待路线确认后单独 PR 添加。

## 数据获取
TODO：B 实现下载脚本，C 核查来源和许可。当前没有数据集或可执行下载入口。字段见 [数据说明](data/README.md)，标签与许可见 [数据卡](data/dataset_card.md)。原始数据只读保存于本地 data/raw/，不提交 Git；共享存储位置和校验哈希 TODO。

## Baseline 运行计划
当前仅有 [配置模板](configs/baseline.yaml)，训练脚本尚未实现，不能据此启动训练。
正式训练前必须填写 target_label、dataset_version、split_version，并冻结指标规则。D 后续通过 PR 实现训练与评估入口，给出经过验证的命令；比较 LR、SVM、RF，报告 Macro-F1、ROC-AUC、PR-AUC、MCC 及跨种子稳定性。禁止在测试集调参。

## 实验记录规范
使用 [实验台账](experiments/registry.csv)，每次正式实验一行，包括失败实验。运行前先提交代码并确认工作区干净；绑定配置、种子、数据版本、划分版本和完整 Git commit。图表追溯到实验编号；指标由脚本生成，禁止手工修改最终指标。详见 [实验规范](experiments/README.md)。

## 项目当前状态
M1 初始化骨架：仅文档、模板、配置和空表头。TODO：独立标签、真实数据、训练入口、归因实现、测试与实验结果。暂无数据规模、模型指标或生物学发现。

## 风险与限制
标签独立性和可得性尚未确认；同源泄漏、类别/物种偏差、坐标错位、归因不稳定和算力不足均需检查。预测质量不达标时优先修正数据和模型；不将低质量模型的归因解释为功能证据。Pfam 未标注不等于无功能；无湿实验时仅报告候选区域。

## TODO 清单
- [ ] 首要 M1 Issue：确定独立预测标签并通过可行性审计（A 主责，B/C/D 复核；待创建 Issue）。
- [ ] 建立成员账号映射、M1—M5 和看板。
- [ ] 实现小样本下载与元数据审计。
- [ ] 冻结标签定义、数据和相似度划分。
- [ ] 实现三种基线与评估入口。
- [ ] 实现归因坐标接口与生物学对照。
- [ ] 由非主责成员复跑核心流程。
