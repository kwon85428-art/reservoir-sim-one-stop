---
name: reservoir-sim-one-stop
display_name: 油气田油藏模拟一条龙
display_name_en: Reservoir Simulation One-Stop
description: |
  油藏数值模拟教学/验证一条龙：合成矿场数据（真值）→ 自研 2D 油水两相 IMPES 模拟器 → 历史拟合（参数扫描）→ 制度变更外推预测，产出含水率拟合图、RMSE 曲线、饱和度快照、预测对比图 + CSV/JSON。
  触发词：油藏模拟、数值模拟、储层模拟、合成数据跑模拟、历史拟合、水驱、见水、含水率、IMPES、黑油模型、注采制度、可辨识性、置信区间、后验、P10P50P90、不确定性量化、跑一轮数据。
  不适用：替代商业模拟器（Eclipse/CMG/tNavigator/HiSimPro）；真实矿场数据的历史拟合生产服务；组分/热采/压裂/全隐式模拟。
agent_created: true
description_zh: |
  油藏数值模拟教学与验证一条龙：用合成矿场数据跑通「真值生成 → 2D 油水两相 IMPES 模拟 → 历史拟合参数扫描 → 注采制度变更外推」完整闭环，产出含水率拟合图、RMSE 曲线、饱和度快照、预测对比图及 CSV/JSON 数据件。内含三个实测物理实现坑（奇异系统、井项符号、CFL 口径）与「渗透率整体尺度不可辨识」判据，可复用于能源 AI 产品拆解验证与油藏工程教学。
description_en: |
  One-stop reservoir simulation teaching and validation pipeline: synthetic field data generation -> self-contained 2D oil-water IMPES simulator -> history matching via parameter sweep -> injection-schedule-change forecast. Outputs water-cut matching plots, RMSE curves, saturation snapshots, prediction comparison charts plus CSV/JSON. Includes three measured implementation pitfalls (singular system, well-term sign, CFL scope) and the permeability-scale non-identifiability criterion, reusable for energy-AI product verification and petroleum engineering education.
category: industry
version: 1.1.0
author: 老桂
tags:
  - 油藏模拟
  - 数值模拟
  - 历史拟合
  - 水驱
  - IMPES
  - 合成数据
  - 石油工程
---

# 油气田油藏模拟一条龙

## 触发场景

- 「用模拟数据跑一轮储层模拟」「跑一轮油藏数值模拟」
- 「演示历史拟合怎么工作」「含水率拟合」「见水时间」
- 「验证这个油藏 AI 产品的拟合/预测逻辑」（拆解判据侧）
- 油藏工程教学：水驱前缘、相渗、井指数、CFL、可辨识性

## 执行步骤

1. **确认参数**（有默认值，可直接跑）：
   网格 `--nx/--ny`（默认 24×24）、井距 `--dx`（m，默认 10）、真值条带对比度 `--c-true`（默认 10）、
   注水量 `--q-inj`（stb/d，默认 3000）、模拟时长 `--t-end`（天，默认 2400）、观测噪声 `--noise`（默认 0.02）、
   随机种子 `--seed`（默认 42）、输出目录 `--out`。
2. **跑五步闭环**（默认全跑）：
   - ① 真值：含高渗条带渗透率场 + 加噪月度采样 → `truth_monthly.csv`
   - ② 盲跑：均质初值模型
   - ③ 历史拟合：扫条带对比度 c，水率 RMSE U 形曲线 → 定 `c_best`
   - ④ 预测：day 2400 起注水量 ×1.5，真值/拟合/盲模型三线对比
   - ⑤ **置信度节点**：χ² 似然 → 参数后验（1σ/95% 区间、ESS）→ 后验加权集合预测带（P10/P50/P90）→ `RMSE_min/σ` 判据
3. **读结果**：`fit_watercut.png`（拟合）、`rmse_curve.png`（可辨识性）、
   `saturation_snapshots.png`（前缘）、`predict.png`（外推）、
   `posterior.png`（参数后验）、`confidence_predict.png`（预测区间）、
   `run_meta.json`（全部参数与结论）。
4. **解读判据**（见 `references/method-and-pitfalls.md`）：
   - RMSE U 形有唯一低点 → 参数可辨识；
   - **不可压缩 + 定率注采下，渗透率整体尺度不可辨识**——历史拟合产品若主打"调 k 倍数"，理论上是在拟合平方向；
   - `RMSE_min/σ ≤ 1.2` → 数据榨干；`> 2` → 结构性失配，点估计不可信；
   - 预测带（P10–P90）宽度才是"敢不敢拿去做方案"的依据，单线预测一律打对折。

## 运行方式

```bash
# 全流程（24×24 网格约 6-8 分钟，置信度节点占一半）
python scripts/sim_pipeline.py --out outputs_sim
# 快速验证（小网格短时长，约 1 分钟）
python scripts/sim_pipeline.py --nx 16 --ny 16 --dx 15 --t-end 1200 --out outputs_quick
```

依赖：numpy / scipy / matplotlib（中文标注用 Microsoft YaHei，非 Windows 需改字体）。

## 输出要求

- 图件 6 张 PNG + 数据件 CSV/JSON，全部落在 `--out` 目录；
- `run_meta.json` 必须含 `C_TRUE / c_best / rmses / posterior（ci1sigma、ci95、ess、rmse_min_over_sigma、verdict）/ identifiability` 字段；
- 对话内总结必须给：真值 c、拟合 c_best、后验 95% 区间、`RMSE_min/σ` 比值与判定四个数。

## 硬边界

- **不替代商业模拟器**：本技能是教学/验证级 toy（不可压缩、无毛管压力、无重力、显式饱和度），
  真实方案设计用 Eclipse / CMG / tNavigator / HiSimPro；
- **不做真实矿场数据拟合服务**：真实数据涉及计量口径清洗、措施作业拆分，需另行处理；
- 不做组分 / 热采 / 压裂模拟（黑油两相而已）；
- 数值结论仅在本合成体系内成立，外推到具体油藏需重跑真实数据。
