# M2 验收：待独立复核与最终验收

核实日期：2026-09-26。基线 `5891c2f7a6eb79e48dd7b75175a3af8176e410a2`，即历史 5891c2f。数据 `gvpa-recognition-c7f6f005d717`；划分 `homology-8b9005e2d9-s42`。

PR #13/#14 已把技术结果归档进 main，不等于满足工程与协作验收。PF00741 命中表示家族坐标，不是功能标签。任务：[M2 post-merge independent audit #15](https://github.com/Lithium9767/MLP/issues/15)。

| 关闭条件 | 当前状态 | 证据或剩余动作 |
| --- | --- | --- |
| 真实独立复核报告进入 main | 待完成 | [independent_review.md](independent_review.md) 仅为模板；B/C身份、独立性待A核实 |
| 最终待验收提交全部测试通过 | 基线预检通过，最终待复核 | 34 tests、validator退出0；每次后续变更须绑定最终完整SHA重新验证 |
| README状态准确 | 本收尾分支已修正，待合并 | 区分全量扫描、60条试样及独立复核 |
| E实验运行commit准确补齐 | 部分解决，待复核 | PR #18新增运行002有真实SHA；旧001仍未知，不冒填666de2b |
| 验收报告更新 | 本收尾分支已更新，待合并 | 本文件不代替A验收 |
| PR #9非作者审查并合并/认可的替代PR | 待完成 | #9仍open，无Review；与main对齐后等待真实审查 |
| 原始输出共享位置、访问与哈希可追溯 | 阻塞 | [shared_artifacts.json](shared_artifacts.json)：位置待补，7份本地文件缺失未验证 |
| 图表与源数据一致 | 队列图已修复，待合并/独立复核 | PR #18的002显示1721/20/305/30/2；其余六图未重跑 |
| 正式扫描代码来源可复现 | 待补证据 | C receipt dirty=true；已核对两个实现文件hash，仍需完整变更清单或干净重跑 |
| A在Issue留下真实最终验收评论 | 待完成 | 助手不代签、不提交冒名Review |

## 技术结果与预检范围

- 候选2078、QC可用2076、primary1721；478簇，discovery1453/335簇、validation623/143簇。
- 本地六份输入/划分文件hash匹配；实际manifest跨集合ID和簇交集为0。原始cluster TSV复现差异已有记录，但不能仅凭版本差异断言全部原因。
- PF00741摘要报告2076 accepted；结构摘要报告7R1C N链88个提交残基、65个已建模残基和39个match states。原始扫描/坐标文件本轮缺失，数字未从原始输出独立重算。
- D敏感性结果已归档；不替换冻结主划分。E七张图已归档，但发现一张确定错误，不能维持“全部数字一致”的旧结论。
- 命令、实测结果、局限详见 [engineering_precheck.md](engineering_precheck.md)，该记录不是B/C独立复核。

## PR核实

| PR | 作者 | 状态/Review | 本轮处理原则 |
| --- | --- | --- | --- |
| #4 | ngocnamd93-spec | open；fengbujue777请求修改 | 不整体合并；另准备ESM-2最小迁移草稿；保留原分支 |
| #7 | ngocnamd93-spec | open；无Review | M2关闭前待启动，不能据计划宣称已做实验 |
| #9 | ngocnamd93-spec | open；无Review | 对齐main、保留验收要求，等真实非作者审查 |
| #11 | fengbujue777 | closed、未合并 | C实现及结构摘要由#13整合，不重复合并 |
| #13 | Lithium9767 | merged；无Review | B/C/D归档；需合并后独立复核 |
| #14 | Lithium9767 | merged；无Review | E归档；需图表修复和来源补证 |

作者账号不等于角色映射或真实运行者；README的身份映射仍待确认。不把新的合并后复核伪装成历史合并前Review。B/C旧分支保留，不删除、不重写。

M3可准备代码草稿和计划；上述关闭门槛全部满足、A亲自验收前，不正式启动M3实验。

## 工程交接链接

- [收尾PR #16](https://github.com/Lithium9767/MLP/pull/16)：`docs/15-close-m2`，本报告、README、台账、预检和独立复核模板。
- [M2计划PR #9](https://github.com/Lithium9767/MLP/pull/9)：`docs/8-m2-team-workplan`，已merge当前main并更新收尾任务，等待真实非作者Review；未合并。
- [M3计划PR #7](https://github.com/Lithium9767/MLP/pull/7)：`docs/6-m3-team-workplan`，已对齐main，明确待启动；未合并。
- [ESM-2最小迁移Draft PR #17](https://github.com/Lithium9767/MLP/pull/17)：`feature/5-m3-esm2-foundation`，仅入口、依赖、ESM测试与迁移审计；37项轻量测试通过，未执行模型推理或M3实验。
- 原PR #4与源分支保留并标注替代关系；没有覆盖M2实现。
- B/C本地快照：`archive/B-82ea140` = `82ea140d6fc811b5bb37167da927e221d369c609`；`archive/C-80892e2` = `80892e2e2fbbbf6c3ef40cbb91d3c0774135770f`。远端旧分支未删除或重写。

后续更新：[修复PR #18](https://github.com/Lithium9767/MLP/pull/18)在main基线上新增运行002，干净运行SHA为 `9c4c1f717e27d90aadbeec4cbe347de0a71e218f`。其待审提交 `2d776fa5341816b95de2c362099a064c23d6f7d3`通过39项测试、validator和差异检查；详细证据在该PR的reports/M2/cohort_fix/目录。C审计17个receipt条目中4匹配、2摘要字节不匹配（CRLF差异已定位）、11未取得；共享访问和dirty来源仍未验证。17是核验条目数，包含两次运行的重复路径，不能和此前13项文件清单直接比较。历史预检保持原样，不把新结果冒充旧检查。先合并真实审查通过的#18，再合并本收尾文档并核对新旧台账行都保留；独立报告和A验收仍为关闭门槛。

工程准备完成，等待独立复核/最终验收；同时须解决上表图表与证据阻塞。以上PR准备和助手测试均不构成成员签字或M2正式关闭。
