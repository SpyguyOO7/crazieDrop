import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import ngsolve as ngs

from solvers import CircuitParams, StaticSolverSession, solve_transient_pulse
from plotting import plot_field_map, plot_transient_summary
from run_exp05 import build_design
import config as cfg

ngs.SetNumThreads(1)

out_dir = Path(__file__).resolve().parent
csv_static = out_dir / "sweep_static.csv"
csv_tr = out_dir / "transient_champions.csv"

if not csv_static.exists() or not csv_tr.exists():
    raise FileNotFoundError("Run run_exp05.py first!")

s_rows = []
with csv_static.open("r", newline="") as f:
    for r in csv.DictReader(f):
        s_rows.append({k: float(v) for k, v in r.items()})

t_rows = []
with csv_tr.open("r", newline="") as f:
    for r in csv.DictReader(f):
        t_rows.append({k: (v if k == "design_label" else float(v)) for k, v in r.items()})


# --- Plot 1: Pareto Front (Dropper Mass vs Repulsive Launch Impulse, Colored by Passive Hold) ---

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.5), dpi=150)

valid = [r for r in s_rows if r["f_hold_gf"] >= 250.0 and r["m_zippy_g"] <= 1.10]
sc = ax1.scatter(
    [r["m_dropper_g"] for r in valid],
    [-r["f_20v_gf"] for r in valid],
    c=[r["f_hold_gf"] for r in valid],
    cmap="plasma",
    vmin=250.0,
    vmax=1000.0,
    s=42,
    edgecolors="black",
    linewidths=0.35,
    alpha=0.88,
)
cb = fig.colorbar(sc, ax=ax1, pad=0.02)
cb.set_label("Passive Hold Force @ I=0A (gf)")
ax1.axhline(0.0, color="black", ls="--", lw=1.0)
ax1.set_xlabel("Drone Dropper Active Mass: Iron Core + 32 AWG Coil (g)")
ax1.set_ylabel("Active Downward Repulsive Shove @ 20V (-F_z in gf)")
ax1.set_title("360-Design Sweep: Dropper Mass vs Active Repulsion")
ax1.grid(True, alpha=0.3)

# Effect of Zippy Iron-Cap Thickness cap_h_mm
cap_hs = sorted(set(r["cap_h_mm"] for r in s_rows))
for ch_val in cap_hs:
    sub = [r for r in s_rows if abs(r["cap_h_mm"] - ch_val) < 1e-4 and abs(r["coil_w_mm"] - 1.7) < 1e-4 and abs(r["r_tip_mm"] - 1.55) < 1e-4]
    sub.sort(key=lambda r: r["m_dropper_g"])
    if sub:
        ax2.scatter(
            [r["f_hold_gf"] for r in sub],
            [-r["f_20v_gf"] for r in sub],
            s=48,
            label=f"Zippy Iron Cap = {ch_val:.2f} mm (m_zip={sub[0]['m_zippy_g']:.2f}g)",
        )

ax2.axhline(0.0, color="black", ls="--", lw=1.0)
ax2.axvline(290.0, color="red", ls=":", lw=1.1, label="290 gf Hold Target")
ax2.set_xlabel("Passive Hold Force @ I=0A (gf)")
ax2.set_ylabel("Active Repulsive Shove @ 20V (-F_z in gf)")
ax2.set_title("How Zippy Iron-Cap Thickness Tunes Hold vs Push-Off")
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=8.0, loc="best")

fig.tight_layout()
fig.savefig(out_dir / "fig_05_pareto_and_cap_thickness.png")
plt.close(fig)


# --- Plot 2: Transient Trimpot Voltage Window Width (12V .. 25V) of the Top 3 Champions ---

fig2, (ax_s, ax_f) = plt.subplots(2, 1, figsize=(8.4, 7.0), dpi=150, sharex=True)
labels = sorted(set(r["design_label"] for r in t_rows))
palette = ["#1b9e77", "#d95f02", "#7570b3"]

best_overall_id = None
best_overall_score = -1.0

for lbl, col in zip(labels, palette):
    sub = [r for r in t_rows if r["design_label"] == lbl]
    sub.sort(key=lambda r: r["v0_v"])
    v_arr = [r["v0_v"] for r in sub]
    sep_arr = [r["sep_40ms_mm"] for r in sub]
    minf_arr = [r["min_f_gf"] for r in sub]
    n_pass = sum(int(r["dropped"]) for r in sub)

    cand_info = f"{lbl} (m_drop={sub[0]['m_dropper_g']:.2f}g, Hold={sub[0]['f_hold_gf']:.0f}gf) — {n_pass}/{len(sub)} Voltages Launch!"
    ax_s.plot(v_arr, sep_arr, marker="o", lw=2.1, color=col, label=cand_info)
    ax_f.plot(v_arr, minf_arr, marker="s", lw=2.0, color=col, label=lbl)

    # Pick overall champion with widest voltage launch window & lightest mass
    score = n_pass * 100.0 + max(sep_arr) - 5.0 * sub[0]["m_dropper_g"]
    if score > best_overall_score:
        best_overall_score = score
        best_overall_id = int(sub[0]["id"])

ax_s.axhline(1.5, color="gray", ls=":", lw=1.1, label="Clean Escape Threshold")
ax_s.set_ylabel("Separation @ 40 ms (mm)")
ax_s.set_title("Exp 05 Champions: Trimpot Voltage Launch Window (12 V -> 25 V, 32 AWG + 1000 uF!)")
ax_s.grid(True, alpha=0.3)
ax_s.legend(fontsize=8.0, loc="upper left")

ax_f.axhline(0.0, color="black", ls="--", lw=1.1, label="0 gf (Below 0 = Active Downward Launch!)")
ax_f.set_ylabel("Peak Repulsive Force (gf)")
ax_f.set_xlabel("Trimpot Capacitor Voltage V0 (V)")
ax_f.grid(True, alpha=0.3)

fig2.tight_layout()
fig2.savefig(out_dir / "fig_05_champions_voltage_window.png")
plt.close(fig2)


# --- Render 2D Field Maps & Transient Launch Plot for the Best Overall Champion ---

if best_overall_id is not None:
    champ = [r for r in s_rows if int(r["id"]) == best_overall_id][0]
    print("\n=== Rendering 2D Field Maps & Launch Waveform for Overall Best Champion ===")
    print(
        f"  ID={int(champ['id'])}: r_shank={champ['r_shank_mm']}mm (D={2*champ['r_shank_mm']:.2f}mm), "
        f"L_core={champ['core_len_mm']}mm, r_tip={champ['r_tip_mm']}mm (d_tip={2*champ['r_tip_mm']:.2f}mm), "
        f"w_coil={champ['coil_w_mm']}mm (N={int(champ['n_turns_32awg'])} turns 32 AWG, OD={champ['coil_od_mm']}mm), "
        f"Zippy cap_h={champ['cap_h_mm']}mm | "
        f"m_dropper={champ['m_dropper_g']:.2f}g, m_zippy={champ['m_zippy_g']:.2f}g, Hold={champ['f_hold_gf']:.1f}gf"
    )

    mdl_c, _ = build_design(
        champ["r_shank_mm"], champ["core_len_mm"], champ["r_tip_mm"], champ["coil_w_mm"], champ["cap_h_mm"]
    )
    sess_c = StaticSolverSession(mdl_c, fe_order=2)
    res_h = sess_c.solve(0.0, compute_diagnostics=True)

    r_loop = champ["r_coil_ohm"] + 0.152
    i_peak_20v = 0.80 * (20.0 / r_loop)
    res_p = sess_c.solve(i_peak_20v, compute_diagnostics=True)

    plot_field_map(
        mdl_c,
        res_h,
        str(out_dir / "fig_05_champion_hold_0A.png"),
        r_max_plot_mm=8.0,
        z_min_plot_mm=-5.5,
        z_max_plot_mm=champ["core_len_mm"] + 2.0,
        title=f"Exp 05 Champion (m_drop={champ['m_dropper_g']:.2f}g) — PASSIVE HOLD: {res_h.force_z_n/9.80665*1e3:+.1f} gf",
    )
    plot_field_map(
        mdl_c,
        res_p,
        str(out_dir / "fig_05_champion_repel_20V.png"),
        r_max_plot_mm=8.0,
        z_min_plot_mm=-5.5,
        z_max_plot_mm=champ["core_len_mm"] + 2.0,
        title=f"Exp 05 Champion — 20V ACTIVE REPULSION (I={i_peak_20v:.2f}A): {res_p.force_z_n/9.80665*1e3:+.1f} gf",
    )

    ckt = CircuitParams(v0_cap_v=20.0, c1_f=cfg.CAP_C1_F, pulse_ms=cfg.PULSE_MS)
    tr_c = solve_transient_pulse(mdl_c, ckt, t_end_ms=40.0, dt_ms=0.2, zippy_mass_g=27.0)
    plot_transient_summary(
        tr_c,
        str(out_dir / "fig_05_champion_transient_20V.png"),
        title=f"Exp 05 Champion (m_drop={champ['m_dropper_g']:.2f}g, 32 AWG, 1000 uF @ 20V): Active Launch!",
    )

print(f"All Exp 05 plots saved to {out_dir}")
