# 协作规范

## 2026-09-29 M2 验收规则调整（当前有效）

项目负责人在本轮明确确认：B/C 独立复核不是老师硬性要求，取消其作为 M2 关闭的必要条件，改由 A 本人验收。本轮 M2 收尾 PR 可按负责人授权、实际技术检查后合并；不补造历史非作者 Review，也不将助手核验署名为成员复核。其他阶段默认协作规则不变。

仍须实际取得原始交接包，核对输入/输出哈希、划分泄漏、图表、代码来源、全部测试及 M2 validator。保留历史 dirty 收据、未知运行 SHA 和复现差异。A 在 Issue #15 亲自确认最终提交、证据和限制后才能关闭 M2。助手可以确认技术检查结果，不能代签 A。

本节优先于本文下方保留的历史 M2 独立复核与非作者 Review 门槛；独立报告改为可选，不再是阻塞项。PR #17 保持 Draft，A 验收前不运行 M3 正式实验。

## 工作流
- main 始终保持可运行；当前骨架阶段尚无训练入口，不宣称模型可运行。
- 首次共享骨架按初始化授权提交 main；后续正式功能禁止直接提交 main。
- 每项任务必须关联 Issue，每个分支从最新 main 创建。
- 合并前至少一名非作者审核；文档约定不等同于已配置 GitHub 分支保护。
- 数据定义变更由 A、B、C 共同确认；评估规则变更由 A、D、E 共同确认。
- 每人使用自己的 GitHub 账号提交，禁止组长代替成员集中提交。
- 不覆盖旧划分；更改冻结评估规则必须重跑可比较的基线。
- 大文件不得进入普通 Git；原始数据、模型、缓存保存于共享存储并记录版本和哈希。
- Word、PPT 指定统稿人，普通协作内容优先使用 Markdown。
- 不提交密钥、个人隐私或未获许可的数据，不改写共享历史。

## 分支与提交
常规分支：docs/<issue-id>-<description>、feature/<issue-id>-<description>、analysis/<issue-id>-<description>、model/<issue-id>-<description>、attribution/<issue-id>-<description>、fix/<issue-id>-<description>。

提交前缀：feat:、data:、model:、analysis:、docs:、test:、fix:、chore:。

首批启动分支按项目约定使用以下名称，仍需在 PR 中关联实际 Issue；由各成员自行创建，不代建提交：
- A：docs/m1-scope
- B：feature/data-ingestion
- C：analysis/label-feasibility
- D：model/baseline-design
- E：attribution/method-interface

```bash
git switch main
git pull origin main
git switch -c <自己的分支名>
```
完成工作并检查暂存内容后：
```bash
git add .
git diff --cached --check
git diff --cached
git commit -m "<符合规范的提交信息>"
git push -u origin <自己的分支名>
```
创建 PR，填写验证、复现命令、版本和非作者复核人。无适用实验的纯文档修改写“不适用”并说明原因。

## 实验与冻结
正式运行前提交代码，检查工作区干净，登记完整 commit 和配置；失败实验保留原因。数据定义、划分、评估、最终结果分别冻结，变更通过 Issue/PR 留痕；禁止手工修改最终指标。结果必须由非主责成员复核。
