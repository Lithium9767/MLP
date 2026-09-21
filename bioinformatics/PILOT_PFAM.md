# PF00741 试标注：可复现操作说明

## 当前状态

Windows PATH和Ubuntu WSL的PATH及常见安装目录未找到InterProScan/HMMER/Pfam数据库。本轮停在 **scan_not_run**。仅生成真实输入、待执行状态和解析代码；raw目录没有伪造结果。60条标签均为空。没有调用远程服务，`run --execute`也只运行本地程序并强制`-dp`关闭远程预计算查询。

## 输入与实际选样

主输入为用户gv.zip内`gv/recognition/data/gvp/gvpa/GvpA_sequences.json`（2078条）。已经核实结构为list：unique_sequence_id、sequence、representative_annotation、cluster_statistics、redundant_members_details。注释含id、description、organism、gene、source_file、full_header、gvp_types。没有独立partial/已审核状态字段，partial从描述/标题提取，conflict从代表及成员类型提取。字段缺失不编造。

背景显式来自同一recognition集合的GvpC/G/J/K/N/O/P JSON，分别233/3576/3135/377/758/1802/220条。它们不是GvpA的另一正式版本，单独用作背景组，类型含GvpA或多类型者不入背景。

实际选择30条清晰GvpA、15条冲突、15条背景。各组organism字符串分别30/15/15个，不代表经过taxonomy核实的独立物种。背景含GvpC 3条，其余六类各2条。长度范围分别59–177、72–157、49–520 aa。

seed=20260908。选择前按来源和记录序号排序，再用固定种子打破并列；优先新增organism、类型、25aa长度区间、来源。完全重复按完整序列哈希去重。无可信的既有聚类版本，使用difflib matching-block近似筛查：较短/较长长度≥0.8，匹配字符数/较短长度≥0.9则不同时选入。本轮选择中10条因该规则排除；这不是完整全库聚类或严格生物序列比对。

40–800 aa是本批操作性质的宽松长度质控阈值，非各家族完整性的生物学定义。明显无效字母、任何代表/成员partial、缺少必要元数据、单一残基比例>60%等隔离，不当阴性。无partial文本不证明完整；最终阴性还需人工核查。

## 本地复现选样

在仓库根目录执行。`archive`改为实际压缩包路径，Python脚本本身没有个人路径；ZIP内容仅读取不解包、不运行。

```powershell
$archive = 'C:/path/to/gv.zip'
$selectionArgs = @('scripts/pilot_pfam.py', 'select', '--input-json', "${archive}::gv/recognition/data/gvp/gvpa/GvpA_sequences.json", '--output-dir', 'data/pilot', '--seed', '20260908', '--targets', '30', '15', '15')
foreach ($family in @('C','G','J','K','N','O','P')) {
  $folder = 'gvp' + $family.ToLower()
  $selectionArgs += @('--background-json', "${archive}::gv/recognition/data/gvp/$folder/Gvp${family}_sequences.json")
}
python @selectionArgs
```

已有data/pilot时必须使用新的输出目录。任何子命令都拒绝覆盖非空输出目录；重新运行可能的人工记录不得覆盖。manifest为UTF-8 CSV，不含序列本体，FASTA只含内部ID和原始完整序列。哈希按原始大写标准氨基酸字符串计算；需要规范化的异常输入隔离，不擅自修正。

## 当前准备命令（不扫描）

```bash
python scripts/pilot_pfam.py run --pilot-dir data/pilot --rules configs/pilot_pfam_rules.json --output-dir results/pilot_pfam
python scripts/pilot_pfam.py summarize --parsed-dir results/pilot_pfam/parsed --output-dir results/pilot_pfam/reports
```

这些命令本轮已经执行。selected=60、submitted=0、success=0、positive=0、not_detected=0、not_submitted=60。这里的0是没有已完成观察，不是60条全部阴性。

## 在已安装InterProScan 5的Linux环境实际执行

本机尚无安装，以下为准确的包装器调用模板，需将路径替换为真实安装位置。先把完整仓库代码和data/pilot放到该Linux环境；涉及其他机器传输由用户自行安排。

```bash
python3 scripts/pilot_pfam.py run \
  --pilot-dir data/pilot \
  --rules configs/pilot_pfam_rules.json \
  --interproscan /path/to/interproscan.sh \
  --execute \
  --output-dir results/pilot_pfam/run_001
```

包装器先运行`--help`与`--version`，要求可识别InterProScan5版本、Pfam版本和参数；失败则不提交。仅在通过后执行如下参数结构（实际绝对路径和退出码写入receipt）：

```text
interproscan.sh -i <pilot_sequences.fasta> -appl Pfam -f TSV,XML -b <run_001/raw/pilot> -dp
```

本轮没有实际扫描命令、扫描工具版本或数据库版本。后续真实值来自安装帮助/版本输出，并与XML中版本交叉核对；不以“最新版”替代版本。

该包装器当前只支持InterProScan5本地XML；不支持直接导入远程结果、InterProScan6或HMMER输出。若安装帮助格式不同，将拒绝运行而非猜参数。

## 真实结果解析与统计

```bash
python3 scripts/pilot_pfam.py parse \
  --pilot-dir data/pilot \
  --run-dir results/pilot_pfam/run_001 \
  --output-dir results/pilot_pfam/run_001/parsed_verified
python3 scripts/pilot_pfam.py summarize \
  --parsed-dir results/pilot_pfam/run_001/parsed_verified \
  --output-dir results/pilot_pfam/run_001/reports_verified
```

`run`创建的parsed是待解析状态；真实解析写新目录，避免覆盖原始状态。只有本地运行receipt确认退出码0、配置/FASTA/manifest/原始XML哈希一致、版本一致、XML存在该ID及原序列和matches容器，才允许判定0/1。TSV只输出有命中的序列，不能靠TSV缺行判阴性。若XML省略无命中序列，则保持result_missing，不猜阴性。

## 固定PF00741规则

版本`pf00741-pilot-v1`，扫描前已保存configs/pilot_pfam_rules.json：

- 必须同时为Pfam和PF00741，名称相似或只有IPR000638不算；保留全部位置。
- 原生InterProScan/Pfam筛选后，对序列E-value和位置domain E-value均要求严格小于1e-5；等于阈值进入边界复核。
- 命中≥30aa；HMM覆盖≥0.5；这些是保守的pilot复核规则，不声称是官方GA阈值。
- 缺少上述显著性或HMM覆盖字段则needs_review；保留空值，不推算分数。sequence coverage只报告，不按全长占比排除多域蛋白。
- 任意目标命中待复核时，整个序列标签空；只有全部目标命中可接受才标1。没有目标命中，且完整运行与输入QC均通过，才标0。
- 0解释为“在当前工具、数据库版本和固定计算规则下未检出合格PF00741命中”，不代表确定不是GvpA或无功能。

统计success定义为扫描及标签判断均完成。完成计算但需复核者单列needs_review，并纳入已提交pending；未提交另列not_submitted，因此满足 submitted=success+pending，success=positive+not_detected。

## 文件与验证

单脚本包含select/run/parse/summarize子命令，复用共享的哈希、CSV及验证逻辑，因此未拆成四个重复包装脚本。

data/pilot含FASTA、manifest、excluded_or_qc_sequences.csv、selection_summary.json、selection.log。results/pilot_pfam含rules.json、run_receipt.json、logs、parsed三张表、reports统计及报告。raw为空是正确的scan_not_run状态。

```bash
python -m unittest discover -s tests -p test_pilot_pfam.py -v
```

测试中的合成XML只存在临时测试目录，用于边界/失败规则验证，绝不当真实扫描结果。全部生成输入和本地结果默认.gitignore，不自动发布用户序列。

官方说明：
- https://interproscan-docs.readthedocs.io/en/v5/HowToRun.html
- https://interproscan-docs.readthedocs.io/en/v5/OutputFormats.html

下一阶段：准备带Pfam数据库的本地InterProScan5安装或由用户在已有Linux服务器执行；取得真实结果后才评估标注流程稳定性。不得将本次准备表述为“已经完成扫描”。
