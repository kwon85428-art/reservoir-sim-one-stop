# -*- coding: utf-8 -*-
"""油气田油藏模拟一条龙：合成数据 -> IMPES 模拟 -> 历史拟合 -> 预测。

模型：2D 油水两相 IMPES（隐式压力 + 显式饱和度），场单位制
（k mD / ft / cp / psi / stb·day）。不可压缩、无毛管压力、无重力。
井控：注水井定率注水；生产井 BHP 约束（Peaceman 井指数）。

可辨识性要点：不可压缩 + 定率注水下，渗透率**整体尺度**不可辨识
（k->k/m 时压力场按 m 缩放、饱和度路径不变）；可辨识的是**非均质
结构**——所以历史拟合参数选"高渗条带对比度 c"（相控建模经典旋钮）。

用法：
  python sim_pipeline.py                          # 默认 24x24 / 2400 天全流程
  python sim_pipeline.py --nx 16 --ny 16 --t-end 1200 --out outputs_quick
"""
import argparse, json, os, time, warnings
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("error", category=RuntimeWarning)  # 奇异系统直接报错，不静默

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

# ---------------- 常量（流体/岩石，与网格无关） ----------------
FT = 3.28084                # m -> ft
MU_W, MU_O = 0.5, 5.0       # cp
SWI, SOR = 0.15, 0.20
PHI = 0.25
RW = 0.328                  # ft，井半径

def kr(S):
    s = np.clip((S - SWI) / (1 - SWI - SOR), 0.0, 1.0)
    return s * s, (1 - s) * (1 - s)

class Model:
    """网格 + 渗透率场 + 面索引。"""
    def __init__(self, nx, ny, dx, c_streak, seed, sigma, base, streak=True):
        self.nx, self.ny, self.dx = nx, ny, dx
        self.dy = self.h = dx
        self.dxf = self.dyf = self.hf = dx * FT
        self.vp_cell = PHI * dx * dx * dx * 6.2898      # rb
        self.y0 = ny // 2
        rng = np.random.default_rng(seed)
        k = base * np.exp(sigma * rng.standard_normal((nx, ny)))
        if streak and c_streak > 1:
            k[:, self.y0-1:self.y0+2] *= c_streak
        self.k = k
        kx = 2.0 * k[:-1, :] * k[1:, :] / (k[:-1, :] + k[1:, :])
        ky = 2.0 * k[:, :-1] * k[:, 1:] / (k[:, :-1] + k[:, 1:])
        self.tx = 0.006328 * kx * (self.dyf * self.hf) / self.dxf
        self.ty = 0.006328 * ky * (self.dxf * self.hf) / self.dyf
        self.re = 0.208 * self.dxf
        self.wi = 0.006328 * 2*np.pi * k[nx-1, self.y0] * self.hf / np.log(self.re / RW)
        I, J = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
        self.fx_i = I[:-1, :].ravel() * ny + J[:-1, :].ravel()
        self.fx_j = I[1:, :].ravel() * ny + J[1:, :].ravel()
        self.fy_i = I[:, :-1].ravel() * ny + J[:, :-1].ravel()
        self.fy_j = I[:, 1:].ravel() * ny + J[:, 1:].ravel()
        self.i_inj = 0 * ny + self.y0
        self.i_prd = (nx - 1) * ny + self.y0

def step(m, S, p, q_inj, dt_cap=30.0):
    krw, kro = kr(S)
    mobw = krw / MU_W
    mobo = kro / (MU_O * 1.0)
    mobt = mobw + mobo

    dp_x = p[:-1, :] - p[1:, :]
    upx = dp_x > 0
    w_x = m.tx * np.where(upx, mobt[:-1, :], mobt[1:, :])
    dp_y = p[:, :-1] - p[:, 1:]
    upy = dp_y > 0
    w_y = m.ty * np.where(upy, mobt[:, :-1], mobt[:, 1:])

    wx, wy = w_x.ravel(), w_y.ravel()
    rows = np.concatenate([m.fx_i, m.fx_i, m.fx_j, m.fx_j, m.fy_i, m.fy_i, m.fy_j, m.fy_j])
    cols = np.concatenate([m.fx_i, m.fx_j, m.fx_j, m.fx_i, m.fy_i, m.fy_j, m.fy_j, m.fy_i])
    vals = np.concatenate([-wx, wx, -wx, wx, -wy, wy, -wy, wy])
    A = coo_matrix((vals, (rows, cols)), shape=(m.nx*m.ny, m.nx*m.ny)).tocsr()
    b = np.zeros(m.nx*m.ny)
    # 约定：A·p|_i = Σ w(p_j-p_i)，质量守恒 Σw(p_i-p_j)=q_ext → A·p = -q_ext
    b[m.i_inj] -= q_inj
    # 生产井 BHP 锚定：Σw(pj-pi) = wi*mobt*(p-pbh) → 对角 -= wi*mobt，b -= wi*mobt*pbh
    pbh = 1500.0
    mt_p = mobt[m.nx-1, m.y0]
    A = A.tolil()
    A[m.i_prd, m.i_prd] -= m.wi * mt_p
    b[m.i_prd] -= m.wi * mt_p * pbh
    p = spsolve(A.tocsr(), b).reshape(m.nx, m.ny)
    if not np.all(np.isfinite(p)):
        raise RuntimeError("压力方程发散（检查源汇符号/奇异系统）")

    mw_x = np.where(upx, mobw[:-1, :], mobw[1:, :])
    qw_x = m.tx * mw_x * dp_x
    mw_y = np.where(upy, mobw[:, :-1], mobw[:, 1:])
    qw_y = m.ty * mw_y * dp_y
    mo_x = np.where(upx, mobo[:-1, :], mobo[1:, :])
    qo_x = m.tx * mo_x * dp_x
    mo_y = np.where(upy, mobo[:, :-1], mobo[:, 1:])
    qo_y = m.ty * mo_y * dp_y

    dS = np.zeros((m.nx, m.ny))
    dS[:-1, :] -= qw_x; dS[1:, :] += qw_x
    dS[:, :-1] -= qw_y; dS[:, 1:] += qw_y
    dS[0, m.y0] += q_inj
    qw_prod = m.wi * mobw[m.nx-1, m.y0] * (p[m.nx-1, m.y0] - pbh)
    dS[m.nx-1, m.y0] -= max(qw_prod, 0.0)

    flux_out = np.zeros((m.nx, m.ny))
    for qf in (qw_x, qo_x):
        flux_out[:-1, :] += np.maximum(qf, 0); flux_out[1:, :] += np.maximum(-qf, 0)
    for qf in (qw_y, qo_y):
        flux_out[:, :-1] += np.maximum(qf, 0); flux_out[:, 1:] += np.maximum(-qf, 0)
    dt = min(dt_cap, 0.6 * m.vp_cell / max(flux_out.max(), 1e-6))
    S = np.clip(S + dt * dS / m.vp_cell, SWI, 1 - SOR)

    krwp, krop = kr(S[m.nx-1, m.y0])
    mw = krwp / MU_W; mo = krop / MU_O
    qo = m.wi * mo * max(p[m.nx-1, m.y0] - pbh, 0)
    qw = max(qw_prod, 0.0)
    return S, p, dt, qo, qw

def run(m, t_end=2400.0, sample_every=30.0, noise=0.0, seed=7, sat_days=(),
        q_inj=3000.0, ramp=None, max_steps=12000):
    S = np.full((m.nx, m.ny), SWI)
    p = np.full((m.nx, m.ny), 3500.0)
    t, rec, snaps, next_s, steps = 0.0, [], {}, 0.0, 0
    rng = np.random.default_rng(seed)
    truncated = False
    while t < t_end - 1e-9:
        if steps >= max_steps:
            truncated = True
            break
        q = q_inj if ramp is None else (q_inj if t < ramp[0] else q_inj * ramp[1])
        S, p, dt, qo, qw = step(m, S, p, q)
        t += dt; steps += 1
        for sd in sat_days:
            if sd not in snaps and t >= sd:
                snaps[sd] = S.copy()
        if t >= next_s:
            wc = qw / max(qo + qw, 1e-9)
            if noise > 0:
                wc = float(np.clip(wc + rng.normal(0, noise), 0, 1))
            rec.append({"day": round(t, 1), "oil_stb_d": round(qo, 1),
                        "water_stb_d": round(qw, 1), "wc": round(wc, 4)})
            next_s += sample_every
    return rec, snaps, steps, truncated

def main():
    ap = argparse.ArgumentParser(description="油藏模拟一条龙：真值→盲跑→拟合→预测")
    ap.add_argument("--nx", type=int, default=24)
    ap.add_argument("--ny", type=int, default=24)
    ap.add_argument("--dx", type=float, default=10.0, help="网格步长 m")
    ap.add_argument("--t-end", type=float, default=2400.0)
    ap.add_argument("--q-inj", type=float, default=3000.0)
    ap.add_argument("--c-true", type=float, default=10.0, help="真值条带对比度")
    ap.add_argument("--noise", type=float, default=0.02, help="观测噪声 std")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="outputs_sim")
    args = ap.parse_args()
    OUT = os.path.abspath(args.out)
    os.makedirs(OUT, exist_ok=True)
    T0 = time.time()

    def fmt_model(c, streak):
        return Model(args.nx, args.ny, args.dx, c, args.seed, 0.5, 80.0, streak=streak)

    print(f"[1/5] 真值模型 c={args.c_true:g} ...")
    m_true = fmt_model(args.c_true, streak=True)
    rec_true, sat_true, n1, tr1 = run(m_true, t_end=args.t_end, noise=args.noise,
                                      sat_days=(args.t_end*0.125, args.t_end*0.3, args.t_end),
                                      q_inj=args.q_inj)
    with open(os.path.join(OUT, "truth_monthly.csv"), "w", encoding="utf-8") as f:
        f.write("day,oil_stb_d,water_stb_d,wc_noisy\n")
        for r in rec_true:
            f.write(f"{r['day']},{r['oil_stb_d']},{r['water_stb_d']},{r['wc']}\n")
    print(f"    {n1} 步 / {len(rec_true)} 采样点 / 截断={tr1} / {time.time()-T0:.0f}s")

    print("[2/5] 盲跑均质 c=1 ...")
    m_blind = fmt_model(1.0, streak=False)
    rec_blind, _, _, _ = run(m_blind, t_end=args.t_end, q_inj=args.q_inj)

    print("[3/5] 拟合：扫条带对比度 ...")
    cands = [1.0, 2.0, 5.0, 10.0, 20.0, 40.0]
    obs = {int(r["day"]): r["wc"] for r in rec_true}
    pairs, recs_c, steps_c = [], {}, {}
    for c in cands:
        rec_c, _, ns, trc = run(fmt_model(c, streak=True), t_end=args.t_end, q_inj=args.q_inj)
        wc_c = {int(r["day"]): r["wc"] for r in rec_c}
        common = sorted(set(obs) & set(wc_c))
        if len(common) < 5:
            print(f"    c={c:g}: 采样重叠不足，跳过"); continue
        e = np.array([wc_c[d] - obs[d] for d in common])
        rmse = float(np.sqrt((e ** 2).mean()))
        pairs.append((c, rmse, len(common))); recs_c[c] = rec_c; steps_c[c] = (ns, trc)
        print(f"    c={c:>5}: RMSE={rmse:.4f}  ({ns} 步, 截断={trc}, 累计 {time.time()-T0:.0f}s)")
    if not pairs:
        raise SystemExit("无有效拟合样本")
    c_best = pairs[int(np.argmin([r for _, r, _ in pairs]))][0]
    print(f"    -> 最优 c={c_best:g}（真值 {args.c_true:g}）")

    print("[4/5] 预测：注水量 x1.5 制度变更外推 ...")
    horizon = args.t_end + 1200.0
    ramp = (args.t_end, 1.5)
    rec_pred_true = run(m_true, t_end=horizon, ramp=ramp, q_inj=args.q_inj)[0]
    rec_pred_fit = run(fmt_model(c_best, streak=True), t_end=horizon, ramp=ramp, q_inj=args.q_inj)[0]
    rec_pred_blind = run(m_blind, t_end=horizon, ramp=ramp, q_inj=args.q_inj)[0]

    # ---- [5/5] 置信度节点：似然 → 参数后验 → 集合预测带 ----
    print("[5/5] 置信度节点：后验权重 + 集合预测带 ...")
    sigma = max(args.noise, 1e-6)
    chi2 = np.array([n * (r / sigma) ** 2 for _, r, n in pairs])
    dchi2 = chi2 - chi2.min()
    w = np.exp(-dchi2 / 2.0); w /= w.sum()
    ess = float(1.0 / np.sum(w ** 2))
    cs = np.array([c for c, _, _ in pairs])

    def wq(q):
        xs = np.log10(cs); order = np.argsort(xs)
        xs_s, w_s = xs[order], w[order]
        return float(10.0 ** np.interp(q, np.cumsum(w_s), xs_s))

    c_lo1, c_hi1 = wq(0.1587), wq(0.8413)      # 1σ 区间
    c_lo2, c_hi2 = wq(0.025), wq(0.975)        # 95% 区间
    c_mean = float(np.sum(w * cs))
    ratio = float(min(r for _, r, _ in pairs) / sigma)   # = RMSE_min/σ，勿用 sqrt(χ²min)（混入√N）
    verdict = ("RMSE 已到噪声底：数据榨干，后验可信"
               if ratio <= 1.2 else
               ("轻度结构失配：预测区间偏窄，勿直接外推" if ratio <= 2.0
                else "结构性失配：模型与数据矛盾，点估计不可信"))

    # 后验加权集合预测带（各候选模型带制度变更跑到 horizon）
    grid_days = np.arange(30.0, horizon, 30.0)
    M = np.zeros((len(pairs), len(grid_days)))
    for i, (c, _, _) in enumerate(pairs):
        rec_m = run(fmt_model(c, streak=True), t_end=horizon, ramp=ramp, q_inj=args.q_inj)[0]
        M[i] = np.interp(grid_days, [r["day"] for r in rec_m], [r["wc"] for r in rec_m])
    order = np.argsort(M, axis=0)
    cw = np.cumsum(w[order], axis=0)
    p10 = np.array([np.interp(0.1, cw[:, j], M[order[:, j], j]) for j in range(M.shape[1])])
    p50 = np.array([np.interp(0.5, cw[:, j], M[order[:, j], j]) for j in range(M.shape[1])])
    p90 = np.array([np.interp(0.9, cw[:, j], M[order[:, j], j]) for j in range(M.shape[1])])
    ci_note = "" if ess >= 2.0 else "（ESS≈1：后验塌缩到单点，区间只是离散候选近似，无直接意义）"
    print(f"    RMSE_min/σ={ratio:.2f} → {verdict}")
    print(f"    c 后验：均值={c_mean:.1f}，1σ=[{c_lo1:.1f},{c_hi1:.1f}]，95%=[{c_lo2:.1f},{c_hi2:.1f}]，ESS={ess:.1f}{ci_note}")

    # ---- 图件 ----
    fig, ax = plt.subplots(figsize=(8, 4.6), dpi=130)
    ax.plot([r["day"] for r in rec_true], [r["wc"] for r in rec_true], "o", ms=2.5,
            color="#888", label="合成观测（真值+噪声）")
    for c, col, lbl in [(1.0, "#E24B4A", "盲模型 c=1"), (40.0, "#BA7517", "c=40（过强）"),
                        (c_best, "#185FA5", f"拟合模型 c={c_best:g}")]:
        if c in recs_c:
            ax.plot([r["day"] for r in recs_c[c]], [r["wc"] for r in recs_c[c]], "-",
                    lw=1.4, color=col, label=lbl)
    ax.set_xlabel("天"); ax.set_ylabel("含水率")
    ax.set_title("历史拟合 = 参数估计：扫高渗条带对比度 c")
    ax.legend(loc="lower right", fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fit_watercut.png")); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4), dpi=130)
    ax.semilogx([c for c, _, _ in pairs], [r for _, r, _ in pairs], "o-", color="#0F6E56")
    ax.axvline(args.c_true, color="#A32D2D", ls="--", lw=1, label=f"真值 c={args.c_true:g}")
    ax.axvline(c_best, color="#185FA5", ls=":", lw=1, label=f"拟合最优 c={c_best:g}")
    ax.set_xlabel("条带对比度 c（对数轴）"); ax.set_ylabel("水率 RMSE")
    ax.set_title("目标函数有唯一低点 → 参数可辨识"); ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "rmse_curve.png")); plt.close(fig)

    if sat_true:
        fig, axs = plt.subplots(1, max(1, len(sat_true)), figsize=(3.7*max(1,len(sat_true)), 3.6), dpi=130)
        axs = np.atleast_1d(axs)
        for ax, sd in zip(axs, sorted(sat_true)):
            im = ax.imshow(sat_true[sd].T, origin="lower", cmap="cividis",
                           vmin=SWI, vmax=1-SOR, extent=[0, m_true.nx*m_true.dx, 0, m_true.ny*m_true.dy],
                           aspect="equal")
            ax.set_title(f"真值含水饱和度 day {sd:g}", fontsize=9)
            ax.plot(0.5*m_true.dx, (m_true.y0+0.5)*m_true.dy, "r^", ms=6)
            ax.plot((m_true.nx-0.5)*m_true.dx, (m_true.y0+0.5)*m_true.dy, "kv", ms=6)
        fig.colorbar(im, ax=axs, shrink=0.85, label="Sw")
        fig.savefig(os.path.join(OUT, "saturation_snapshots.png")); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.4), dpi=130)
    for rec, col, lbl in [(rec_pred_true, "#444", "真值模型"),
                          (rec_pred_fit, "#185FA5", f"拟合模型 c={c_best:g}"),
                          (rec_pred_blind, "#E24B4A", "盲模型 c=1")]:
        ax.plot([r["day"] for r in rec], [r["wc"] for r in rec], "-", lw=1.4, color=col, label=lbl)
    ax.axvline(args.t_end, color="#BA7517", ls="--", lw=1)
    ax.text(args.t_end + 30, 0.05, "注水量 x1.5\n（制度变更）", fontsize=8, color="#854F0B")
    ax.set_xlabel("天"); ax.set_ylabel("含水率")
    ax.set_title("预测：拟合模型贴真值，盲模型在制度变更后跑偏")
    ax.legend(fontsize=8, loc="lower right"); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "predict.png")); plt.close(fig)

    # 后验 Δχ² 曲线
    fig, ax = plt.subplots(figsize=(6, 4), dpi=130)
    ax.semilogx(cs, dchi2, "o-", color="#0F6E56")
    ax.axhline(1.0, color="#888", ls=":", lw=1)
    ax.axhline(3.84, color="#A32D2D", ls="--", lw=1)
    ax.text(cs[0] * 1.1, 1.15, "Δχ²=1（1σ）", fontsize=7, color="#666")
    ax.text(cs[0] * 1.1, 4.1, "Δχ²=3.84（95%）", fontsize=7, color="#A32D2D")
    ax.axvline(args.c_true, color="#A32D2D", ls="--", lw=1, label=f"真值 {args.c_true:g}")
    ax.axvline(c_best, color="#185FA5", ls=":", lw=1, label=f"点估计 {c_best:g}")
    ax.axvspan(c_lo2, c_hi2, color="#B5D4F4", alpha=0.3, label=f"95% 区间 [{c_lo2:.1f},{c_hi2:.1f}]")
    ax.set_xlabel("条带对比度 c（对数轴）"); ax.set_ylabel("Δχ²")
    ax.set_title("参数后验：阈值线与区间宽度即置信度")
    ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "posterior.png")); plt.close(fig)

    # 置信度预测带
    fig, ax = plt.subplots(figsize=(8, 4.4), dpi=130)
    ax.fill_between(grid_days, p10, p90, color="#185FA5", alpha=0.22,
                    label="P10–P90（后验集合）")
    ax.plot(grid_days, p50, "-", lw=1.6, color="#185FA5", label="P50")
    ax.plot([r["day"] for r in rec_pred_true], [r["wc"] for r in rec_pred_true],
            "-", lw=1.3, color="#444", label="真值模型")
    ax.plot([r["day"] for r in rec_pred_blind], [r["wc"] for r in rec_pred_blind],
            "-", lw=1.3, color="#E24B4A", label="盲模型 c=1")
    ax.axvline(args.t_end, color="#BA7517", ls="--", lw=1)
    ax.text(args.t_end + 30, 0.05, "注水量 x1.5", fontsize=8, color="#854F0B")
    ax.set_xlabel("天"); ax.set_ylabel("含水率")
    ax.set_title(f"预测区间：RMSE_min/σ={ratio:.2f}，{verdict}")
    ax.legend(fontsize=8, loc="lower right"); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "confidence_predict.png")); plt.close(fig)

    meta = {"C_TRUE": args.c_true, "c_best": c_best,
            "rmses": {str(c): r for c, r, _ in pairs},
            "grid": f"{args.nx}x{args.ny}", "dx_m": args.dx,
            "q_inj_stb_d": args.q_inj, "T_END_days": args.t_end,
            "noise_std": args.noise, "seed": args.seed,
            "steps": {"truth": n1, "truncated": tr1},
            "physics": "2D oil-water IMPES, incompressible, no Pc/gravity, producer BHP",
            "identifiability": "定率注采+不可压缩：k 整体尺度不可辨识，非均质结构可辨识",
            "posterior": {"sigma": sigma, "ess": ess,
                          "c_mean": c_mean, "ci1sigma": [c_lo1, c_hi1],
                          "ci95": [c_lo2, c_hi2],
                          "rmse_min_over_sigma": ratio, "verdict": verdict}}
    with open(os.path.join(OUT, "run_meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=1)
    print(f"完成 {time.time()-T0:.0f}s -> {OUT}")

if __name__ == "__main__":
    main()
