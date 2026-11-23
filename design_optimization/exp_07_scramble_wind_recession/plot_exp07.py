import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import ngsolve as ngs

from solvers import CircuitParams, StaticSolverSession, solve_transient_pulse
from plotting import plot_field_map, plot_transient_summary
from run_exp07 import build_recessed_design
import config as cfg

ngs.SetNumThreads(1)

out_dir = Path(__file__).resolve().parent
s_csv = out_dir / "static_recession_grid.csv"
t_csv = out_dir / "transient_recession_results.csv"

if not s_csv.exists() or not t_csv.exists():
    raise FileNotFoundError("Run run_exp07.py first!")

s_rows = []
with s_csv.open("r", newline="") as f:
    for r in csv.DictReader(f):
        s_rows.append({k: float(v) for k, v in r.items()})

t_rows = []
with t_csv.open("r", newline="") as f:
    for r in csv.DictReader(f):
        t_rows.append({k: (v if k == "series" else float(v)) for k, v in r.items()})


# --- Plot 1: Static Sensitivity to Bottom Recession z_coil_bot & Scramble Fill Factor ---

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 5.2), dpi=150)

styles = [
    (8.5, 1.8, 0.50, "#1b9e77", "-",  "L=8.5mm, w=1.8mm, Fill=0.50 (Semi-Neat)"),
    (8.5, 1.8, 0.42, "#377eb8", "-",  "L=8.5mm, w=1.8mm, Fill=0.42 (Typical Blitz)"),
    (8.5, 1.8, 0.35, "#e41a1c", "--", "L=8.5mm, w=1.8mm, Fill=0.35 (Very Messy)"),
    (8.5, 2.3, 0.35, "#984ea3", "-",  "L=8.5mm, w=2.3mm, Fill=0.35 (Messy + Wider Bobbin)"),
    (7.5, 1.8, 0.42, "#ff7f00", ":",  "L=7.5mm, w=1.8mm, Fill=0.42 (Shorter Core)"),
]

for cl, cw, fl, col, ls, lbl in styles:
    sub = [
        r for r in s_rows
        if abs(r["core_len_mm"] - cl) < 1e-4
        and abs(r["coil_w_mm"] - cw) < 1e-4
        and abs(r["fill_factor"] - fl) < 1e-4
    ]
    sub.sort(key=lambda r: r["z_coil_bot_mm"])
    if not sub:
        continue
    zb = [r["z_coil_bot_mm"] for r in sub]
    ax1.plot(zb, [-r["f_18v_gf"] for r in sub], marker="o", lw=2.0, color=col, ls=ls, label=lbl)
    ax2.plot(zb, [r["impulse_18v_gf_ms"] for r in sub], marker="s", lw=2.0, color=col, ls=ls, label=lbl)

ax1.axhline(0.0, color="black", ls="--", lw=1.0)
ax1.set_xlabel("Bottom Coil Recession z_coil_bot (mm)")
ax1.set_ylabel("Repulsive Shove @ 18V (-F_z in gf)")
ax1.set_title("How Bottom Coil Recession Affects Repulsive Force")
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=7.8, loc="best")

ax2.set_xlabel("Bottom Coil Recession z_coil_bot (mm)")
ax2.set_ylabel("18V Repulsive Impulse Proxy (-F_18V * tau_RC, gf*ms)")
ax2.set_title("How Scramble Fill & Bobbin Width Affect Total Push Impulse")
ax2.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(out_dir / "fig_07_recession_and_fill_sensitivity.png")
plt.close(fig)


# --- Plot 2: Moving-Gap Transient Voltage Launch Window Across Recession & Scramble Fill ---

fig2, (ax_sep, ax_f) = plt.subplots(2, 1, figsize=(9.0, 7.5), dpi=150, sharex=True)

tr_meta = [
    ("Recess_0p5mm_Fill0p42_w1p8",    "#1b9e77", "-",  "z_bot=0.5mm, Fill=0.42, w=1.8mm"),
    ("Recess_0p9mm_Fill0p42_w1p8",    "#377eb8", "-",  "z_bot=0.9mm, Fill=0.42, w=1.8mm"),
    ("Recess_1p3mm_Fill0p42_w1p8",    "#ff7f00", "-",  "z_bot=1.3mm, Fill=0.42, w=1.8mm"),
    ("Recess_1p7mm_Fill0p42_w1p8",    "#a65628", "--", "z_bot=1.7mm, Fill=0.42, w=1.8mm"),
    ("MessyFill_0p35_Recess1p3_w1p8", "#e41a1c", ":",  "z_bot=1.3mm, Messy Fill=0.35, w=1.8mm"),
    ("MessyFill_0p35_Recess1p3_w2p3", "#984ea3", "-.", "z_bot=1.3mm, Messy Fill=0.35, w=2.3mm (Wider Bobbin)"),
]

for key, col, ls, desc in tr_meta:
    sub = [r for r in t_rows if r["series"] == key]
    if not sub:
        continue
    sub.sort(key=lambda r: r["v0_v"])
    n_pass = sum(int(r["dropped"]) for r in sub)
    full_lbl = f"{desc} (N={int(sub[0]['n_turns'])}t, m={sub[0]['m_dropper_g']:.2f}g) -> {n_pass}/{len(sub)} Launch"
    ax_sep.plot([r["v0_v"] for r in sub], [r["sep_42ms_mm"] for r in sub], marker="o", lw=2.0, color=col, ls=ls, label=full_lbl)
    ax_f.plot([r["v0_v"] for r in sub], [r["min_f_gf"] for r in sub], marker="s", lw=1.8, color=col, ls=ls, label=desc)

ax_sep.axhline(1.5, color="gray", ls=":", lw=1.1)
ax_sep.set_ylabel("Separation @ 42 ms (mm)")
ax_sep.set_title("Exp 07 Transient Proof: Scramble Hand-Winding + Bottom Lead/Bobbin Recession!")
ax_sep.grid(True, alpha=0.3)
ax_sep.legend(fontsize=7.6, loc="upper left")

ax_f.axhline(0.0, color="black", ls="--", lw=1.1)
ax_f.set_ylabel("Peak Repulsive Force (gf)")
ax_f.set_xlabel("Trimpot Capacitor Voltage V0 (V)")
ax_f.grid(True, alpha=0.3)

fig2.tight_layout()
fig2.savefig(out_dir / "fig_07_transient_recession_margin.png")
plt.close(fig2)


# --- Render 2D Field Map & 20V Launch Plot for Practical Shop Design (z_bot = 1.3mm, Scramble Fill = 0.42) ---

mdl_shop = build_recessed_design(core_len_mm=8.5, coil_w_mm=1.8, z_coil_bot_mm=1.30, fill_factor=0.42)
sess = StaticSolverSession(mdl_shop, fe_order=2)
r_coil = mdl_shop.parts[1].coil_res_ohm
m_shop = mdl_shop.parts[0].mass_g + mdl_shop.parts[1].mass_g

res_hold = sess.solve(0.0, compute_diagnostics=True)
i_20v = 0.84 * (20.0 / (r_coil + 0.152))
res_rep = sess.solve(i_20v, compute_diagnostics=True)

plot_field_map(
    mdl_shop,
    res_hold,
    str(out_dir / "fig_07_shop_build_hold_0A.png"),
    r_max_plot_mm=8.0,
    z_min_plot_mm=-5.5,
    z_max_plot_mm=10.5,
    title=f"Shop-Ready Build (Flush Top, 1.3mm Bottom Recess, Fill=0.42, m={m_shop:.2f}g) — HOLD: {res_hold.force_z_n/9.80665*1e3:+.0f} gf",
)
plot_field_map(
    mdl_shop,
    res_rep,
    str(out_dir / "fig_07_shop_build_repel_20V.png"),
    r_max_plot_mm=8.0,
    z_min_plot_mm=-5.5,
    z_max_plot_mm=10.5,
    title=f"Shop-Ready Build — 20V REPULSION w/ 1.3mm Bare Nose (I={i_20v:.2f}A): {res_rep.force_z_n/9.80665*1e3:+.1f} gf",
)

ckt_20 = CircuitParams(v0_cap_v=20.0, c1_f=cfg.CAP_C1_F, pulse_ms=cfg.PULSE_MS)
tr_shop = solve_transient_pulse(mdl_shop, ckt_20, t_end_ms=42.0, dt_ms=0.2, zippy_mass_g=27.0)
plot_transient_summary(
    tr_shop,
    str(out_dir / "fig_07_shop_build_transient_20V.png"),
    title=f"Shop-Ready Build (z_bot=1.3mm, Scramble Fill=0.42, 34 AWG, m={m_shop:.2f}g @ 20V): Launch!",
)

print(f"All Exp 07 plots saved to {out_dir}")
