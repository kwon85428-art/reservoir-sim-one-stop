# RELEASE-NOTES · reservoir-sim-one-stop

格式：语义化版本。四件套对齐口径：SKILL.md frontmatter `version` = README 元数据表「当前版本」= 本文件最新条目 = git tag。

## v1.1.1 — 2026-10-01

- 页面按 skill-page-design-template 八节结构重构：§0 三秒路由、符号约定、判据体系表、硬边界路由到具体 skill
- 仓库一致性补齐：`.gitattributes`（LF 归一 + 二进制白名单）、`RELEASE-NOTES.md`、`assets/workbuddy-invite-poster.png`
- 元数据表补齐仓库/底本/知识规模/测试集字段；SKILL↔README 版本对齐
- 无功能变更（sim_pipeline.py 与 v1.1.0 相同）

## v1.1.0 — 2026-10-01

- 新增第 ⑤ 步**置信度节点**：
  - χ² 似然 `N·(RMSE/σ)²` → 后验权重 `exp(−Δχ²/2)`
  - log10(c) 加权 CDF 反插值出参数 1σ / 95% 区间
  - 后验加权候选模型集合 → 含水率 P10/P50/P90 预测带（含制度变更段）
  - `RMSE_min/σ` 判据：≤1.2 数据榨干 / 1.2–2 轻度结构失配 / >2 结构性失配
  - ESS<2 时自动标注"区间为离散候选近似，无直接意义"
- 新增图件：`posterior.png`（Δχ² 曲线 + 阈值线 + 95% 区间带）、`confidence_predict.png`（预测带）
- 修复：`ratio` 误用 `sqrt(χ²_min)`（混入 √N）；候选跳过后 pairs 元组解包错位
- 触发词扩充：置信区间 / 后验 / P10P50P90 / 不确定性量化
- 实测验证：快速档（16×16/dx15/1200d）54s 全五步通过，RMSE_min/σ=0.90 判"已到噪声底"，ESS=1.0 如实标注

## v1.0.0 — 2026-09-30

- 首版四步闭环：真值生成（lognormal 渗透率场 + 高渗条带 + 加噪月度采样）→ 盲跑（均质）→ 历史拟合（扫条带对比度 c，水率 RMSE U 形）→ 制度变更预测（注水量 ×1.5，真值/拟合/盲三线）
- 自研 2D 油水两相 IMPES 模拟器（场单位制，BHP 生产井 + Peaceman 井指数）
- 三个实测实现坑入档：奇异系统 NaN 静默卡死 / 井项符号反压力倒挂 / CFL 双相口径
- 核心判据：不可压缩 + 定率注采下渗透率整体尺度不可辨识，可辨识的是非均质结构
