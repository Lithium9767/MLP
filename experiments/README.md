# 实验台账
registry.csv 当前只有表头，没有实验结果。每次正式实验记录一行，experiment_id 唯一；复跑分配新的运行编号。

- 实验运行前先提交代码并确认工作区干净，保存运行时完整 Git commit。
- 每次实验绑定配置、种子、数据版本、划分版本和 Git commit。
- status 建议使用 running、completed、failed；失败实验也应保留必要记录及错误原因。
- rq 填 RQ1/RQ2/RQ3 或“模型质量保障”；owner 使用实际运行者，当前成员映射 TODO。
- metrics_path 指向脚本生成的小型指标文件；artifact_uri 指向共享存储的固定版本，存储位置 TODO，不能包含凭据。
- conclusion 和 limitations 由实际证据填写；未完成时留空，不虚构结果。
- 最终图表必须能够追溯到实验编号、配置和数据。
- 禁止手工修改最终指标；发现错误应修复脚本、重新运行并保留变更记录。
- CSV 含逗号或换行的内容按标准 CSV 引号转义。

字段：experiment_id、owner、status、rq、dataset_version、split_version、config_path、seed、git_commit、metrics_path、artifact_uri、conclusion、limitations。
实验入口及自动汇总脚本 TODO，由 D/E 实现，A 复核 RQ 对应关系。
