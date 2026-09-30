# 油气田油藏模拟一条龙 · Reservoir Simulation One-Stop

用**合成矿场数据**跑通油藏数值模拟的完整闭环，每一步都可核、可复现：

```
① 真值生成          ② 盲跑            ③ 历史拟合              ④ 制度变更预测        ⑤ 置信度节点
渗透率场(含高渗条带)  均质初值模型      扫条带对比度 c            注水量 ×1.5 外推      χ² 似然 → 参数后验
+ 加噪月度采样  →    不知道有条带  →   水率 RMSE U 形曲线   →   真值/拟合/盲 三线  →   P10/P50/P90 预测带
truth_monthly.csv                    c_best                    拟合收敛才可外推       RMSE_min/σ 判据
```

## 快速开始

```bash
pip install numpy scipy matplotlib
# 快速验证（约 1 分钟）
python scripts/sim_pipeline.py --nx 16 --ny 16 --dx 15 --t-end 1200 --out outputs_quick
# 全流程（24×24，约 6-8 分钟）
python scripts/sim_pipeline.py --out outputs_sim
```

产出 6 张图 + CSV/JSON：含水率拟合、RMSE 曲线、饱和度前缘快照、预测对比、**参数后验 Δχ² 曲线**、**P10–P90 预测带**。

## 模型与物理

- 2D 油水两相 IMPES（隐式压力 + 显式饱和度），场单位制（mD/ft/cp/psi/stb·day）
- 不可压缩、无毛管压力、无重力——**教学/验证级 toy，不替代 Eclipse/CMG/tNavigator/HiSimPro**
- Brooks-Corey 相渗，Peaceman 井指数，BHP 生产井 + 定率注水井

## 两个核心判据（拆解油藏 AI 产品时直接用）

**1. 可辨识性**：不可压缩 + 定率注采下，渗透率**整体尺度**不可辨识（k→k/m 时压力场按 m 缩放、饱和度路径不变），可辨识的只有**非均质结构**。历史拟合产品若主打"调 k 倍数"，理论上是在拟合一个平方向。

**2. 置信度三查**：
- `RMSE_min/σ ≤ 1.2` → 数据榨干；`> 2` → 结构性失配，点估计不可信
- 拟合结果带参数置信区间吗？（后验 ESS 小于 2 时区间是离散近似）
- 预测输出是单线还是 P10–P50–P90 带？单线 = 缺 UQ，打对折

## 实测坑（三个，每个都能让第一次跑全错）

1. **奇异系统**：双定率井 → 压力零空间 → `spsolve` 返 NaN → 时间循环空转不报错。处置：RuntimeWarning 升级为 error。
2. **井项符号反**：装配约定 `A·p = −q_ext` 时注入源 b 取负、BHP 井对角 `-= wi·mobt`。症状：压力场倒挂水堵注入井 / 全场漂负。
3. **CFL 口径**：显式饱和度步长限制来自水相+油相面通量之和，只算水相会提前失稳。

详见 `references/method-and-pitfalls.md`。

## 边界

不做真实矿场数据的历史拟合服务（计量口径清洗、措施拆分需另行处理）；不做组分/热采/压裂；数值结论仅在本合成体系内成立。

---
author: 老桂（@kwon85428-art）｜ WorkBuddy Skill：`reservoir-sim-one-stop` v1.1.0
