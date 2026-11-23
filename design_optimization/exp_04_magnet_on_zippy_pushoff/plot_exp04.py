import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import ngsolve as ngs

from solvers import CircuitParams, StaticSolverSession, solve_transient_pulse
from plotting import plot_field_map, plot_transient_summary
from run_exp04 import build_architecture
import config as cfg

ngs.SetNumThreads(1)

out_dir = Path(__file__).resolve().parent
static_csv = out_dir / "static_curves.csv"
tr_csv = out_dir / "transient_margin.csv"

if not static_csv.exists() or not tr_csv.exists():
    raise FileNotFoundError("Run run_exp04.py first to generate static_curves.csv and transient_margin.csv!")


# --- Load CSV Data ---

s_rows = []
with static_csv.open("r", newline="") as f:
    for r in csv.DictReader(f):
        s_rows.append({k: (v if k == "arch" else float(v)) for k, v in r.items()})

t_rows = []
with tr_csv.open("r", newline="") as f:
    for r in csv.DictReader(f):
        t_rows.append({k: (v if k == "arch" else float(v)) for k, v in r.items()})


# --- Plot 1: Static Force F_z(I) Across All 7 Magnet-on-Zippy Architectures ---

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.2, 8.6), dpi=150, sharex=True)

arch_labels = {
    "inv_n52_touch_15um":     ("Bare N52 on Zippy @ 18um Touch (Huge Mirror Hold!)", "#e41a1c", "-"),
    "inv_n52_spacer_0p55mm":  ("Inverted: N52 on Zippy + 0.55mm Spacer (3.15g Dropper!)", "#1b9e77", "-"),
    "inv_n52_iron_behind":    ("Inverted: N52 + 1mm Back-Iron on Zippy + 0.65mm Spacer", "#377eb8", "--"),
    "inv_iron_cap_on_n52":    ("Iron Cap IN FRONT of Zippy N52 @ 15um (Re-Attracts!)", "#ff7f00", "-."),
    "inv_weaker_n35_0p40mm":  ("Inverted: Weaker N35 on Zippy + 0.40mm Spacer", "#984ea3", "-"),
    "inv_ferrite_touch_15um": ("Inverted: Weak Ferrite C8 (Br=0.39T) @ 15um Touch", "#a65628", ":"),
    "dual_n52_spacer_0p85mm": ("Dual N52 (Dropper + Zippy) + 0.85mm Spacer", "#666666", "--"),
}

for arch, (lbl, col, ls) in arch_labels.items():
    sub = [r for r in s_rows if r["arch"] == arch]
    if not sub:
        continue
    cur = [r["current_a"] for r in sub]
    fgf = [r["force_z_gf"] for r in sub]
    dmg = [r["pm_demag_margin_ka_m"] for r in sub]
    ax1.plot(cur, fgf, label=lbl, color=col, ls=ls, lw=2.0)
    ax2.plot(cur, dmg, label=lbl, color=col, ls=ls, lw=1.8)

ax1.axhline(27.0, color="orange", ls=":", lw=1.2, label="+27 gf (Gravity Release Threshold)")
ax1.axhline(0.0, color="black", ls="--", lw=1.2, label="0 gf (Below 0 = ACTIVE MAGNETIC REPULSION!)")
ax1.set_ylim(-450.0, 600.0)
ax1.set_ylabel("Force on Zippy F_z (gf)\n[+ = Attract, - = Active Push-Off!]")
ax1.set_title("Exp 04: Magnet-on-Zippy Static Force vs 32 AWG Coil Current (Negative = Active Launch!)")
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=7.6, loc="upper right")

ax2.axhline(0.0, color="red", ls="--", lw=1.2, label="0 kA/m (80°C N52 Knee; 25°C N52 Safe Down to ~ -100 kA/m)")
ax2.set_ylabel("PM Demag Margin (kA/m)")
ax2.set_xlabel("32 AWG Coil Current I (A)  [~197 turns, R_coil = 2.95 Ohm]")
ax2.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(out_dir / "fig_04_static_pushoff_comparison.png")
plt.close(fig)


# --- Plot 2: Trimpot Voltage Margin Width (Separation @ 40 ms vs Trimpot V0 with Existing 32 AWG + 1000 uF!) ---

fig2, (ax_sep, ax_minf) = plt.subplots(2, 1, figsize=(8.2, 6.8), dpi=150, sharex=True)

for arch in ["inv_n52_spacer_0p55mm", "inv_weaker_n35_0p40mm", "inv_iron_cap_on_n52", "dual_n52_spacer_0p85mm"]:
    sub = [r for r in t_rows if r["arch"] == arch]
    if not sub:
        continue
    lbl, col, ls = arch_labels[arch]
    v_vals = [r["v0_v"] for r in sub]
    ax_sep.plot(v_vals, [r["sep_40ms_mm"] for r in sub], marker="o", lw=2.0, color=col, ls=ls, label=lbl)
    ax_minf.plot(v_vals, [r["min_force_gf"] for r in sub], marker="s", lw=2.0, color=col, ls=ls, label=lbl)

ax_sep.axhline(1.5, color="gray", ls=":", lw=1.1, label="Clean Escape Threshold")
ax_sep.set_ylabel("Separation @ 40 ms (mm)")
ax_sep.set_title("Trimpot Voltage Tolerance Window (Existing 32 AWG Wire + 1000 uF Cap!)")
ax_sep.grid(True, alpha=0.3)
ax_sep.legend(fontsize=8.0, loc="best")

ax_minf.axhline(0.0, color="black", ls="--", lw=1.1, label="0 gf (Negative = Active Downward Kick!)")
ax_minf.set_ylabel("Peak Repulsive Kick (gf)")
ax_minf.set_xlabel("Trimpot Capacitor Voltage V0 (V)")
ax_minf.grid(True, alpha=0.3)

fig2.tight_layout()
fig2.savefig(out_dir / "fig_04_trimpot_voltage_margin.png")
plt.close(fig2)


# --- Render 2D Field Maps & Transient Waveform of the Actual Winning Design (inv_iron_cap_on_n52) ---

best_arch = "inv_iron_cap_on_n52"
mdl_best = build_architecture(best_arch)
sess_best = StaticSolverSession(mdl_best, fe_order=2)

res_hold = sess_best.solve(0.0, compute_diagnostics=True)
res_repel = sess_best.solve(5.6, compute_diagnostics=True)  # ~22V peak current in 32 AWG!

plot_field_map(
    mdl_best,
    res_hold,
    str(out_dir / "fig_04_winner_iron_cap_hold_0A.png"),
    r_max_plot_mm=8.5,
    z_min_plot_mm=-6.0,
    z_max_plot_mm=9.5,
    title=(
        f"WINNER: Iron-Cap + N52 on Zippy (m_dropper=3.15g!) — HOLD (I=0A): "
        f"F_z = {res_hold.force_z_n/9.80665*1e3:+.1f} gf"
    ),
)

plot_field_map(
    mdl_best,
    res_repel,
    str(out_dir / "fig_04_winner_iron_cap_repel_5p6A.png"),
    r_max_plot_mm=8.5,
    z_min_plot_mm=-6.0,
    z_max_plot_mm=9.5,
    title=(
        f"WINNER: Iron-Cap + N52 on Zippy — ACTIVE REPULSION (I=5.6A): "
        f"F_z = {res_repel.force_z_n/9.80665*1e3:+.1f} gf (SHOOTS DOWNWARD!)"
    ),
)

# Detailed 22V transient launch waveform for the winning Iron-Cap + N52 combo
ckt_22v = CircuitParams(v0_cap_v=22.0, c1_f=cfg.CAP_C1_F, pulse_ms=cfg.PULSE_MS)
tr_best = solve_transient_pulse(mdl_best, ckt_22v, t_end_ms=40.0, dt_ms=0.2, zippy_mass_g=27.0)
plot_transient_summary(
    tr_best,
    str(out_dir / "fig_04_winner_iron_cap_transient_22V.png"),
    title="WINNER: Iron-Cap + N52 on Zippy (32 AWG, 1000 uF @ 22V): Active Launch Waveform!",
)

print(f"All Exp 04 plots saved to {out_dir}")
