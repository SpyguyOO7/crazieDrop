import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import ngsolve as ngs

from solvers import CircuitParams, StaticSolverSession, solve_transient_pulse
from plotting import plot_field_map, plot_transient_summary
from run_exp06 import build_design
import config as cfg

ngs.SetNumThreads(1)

out_dir = Path(__file__).resolve().parent
csv_static = out_dir / "sweep_static.csv"
csv_tr = out_dir / "transient_champions.csv"

if not csv_static.exists() or not csv_tr.exists():
    raise FileNotFoundError("Run run_exp06.py first!")

s_rows = []
with csv_static.open("r", newline="") as f:
    for r in csv.DictReader(f):
        s_rows.append({k: float(v) for k, v in r.items()})

t_rows = []
with csv_tr.open("r", newline="") as f:
    for r in csv.DictReader(f):
        t_rows.append({k: (v if k == "design_label" else float(v)) for k, v in r.items()})


# --- Plot 1: Trimpot Voltage Window (12V .. 25V) Across the 3 Bare-Contact 34-AWG Champions + 470uF Check ---

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.0, 7.4), dpi=150, sharex=True)

main_labels = [
    ("Featherweight_Sub2p2g", "#e7298a", "-"),
    ("Balanced_Sub2p8g",      "#1b9e77", "-"),
    ("HeavyDuty_Sub3p6g",     "#7570b3", "-"),
    ("Balanced_LightCap_470uF", "#d95f02", "--"),
]

best_champ_id = None
best_pass_count = -1

for lbl, col, ls in main_labels:
    sub = [r for r in t_rows if r["design_label"] == lbl]
    if not sub:
        continue
    sub.sort(key=lambda r: r["v0_v"])
    v_arr = [r["v0_v"] for r in sub]
    sep_arr = [r["sep_42ms_mm"] for r in sub]
    minf_arr = [r["min_f_gf"] for r in sub]
    n_pass = sum(int(r["dropped"]) for r in sub)

    desc = (
        f"{lbl} (m={sub[0]['m_dropper_g']:.2f}g, C1={int(sub[0]['cap_uf'])}uF, "
        f"Hold={sub[0]['hold_gf']:.0f}gf) -> {n_pass}/{len(sub)} Launch!"
    )
    ax1.plot(v_arr, sep_arr, marker="o", lw=2.1, ls=ls, color=col, label=desc)
    ax2.plot(v_arr, minf_arr, marker="s", lw=1.9, ls=ls, color=col, label=lbl)

    if "LightCap" not in lbl and n_pass > best_pass_count:
        best_pass_count = n_pass
        best_champ_id = int(sub[0]["id"])

ax1.axhline(1.5, color="gray", ls=":", lw=1.1, label="Clean Escape Threshold")
ax1.set_ylabel("Separation @ 42 ms (mm)")
ax1.set_title("Exp 06 (34 AWG, Bare 15um Lathe Contact): Trimpot Launch Window Width!")
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=8.0, loc="upper left")

ax2.axhline(0.0, color="black", ls="--", lw=1.1)
ax2.set_ylabel("Peak Repulsive Force (gf)")
ax2.set_xlabel("Trimpot Capacitor Voltage V0 (V)")
ax2.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(out_dir / "fig_06_trimpot_launch_window.png")
plt.close(fig)


# --- Plot 2: Verification of Optional Bench Kapton Shimming (15um Bare vs 65um vs 140um) ---

fig_s, ax_s = plt.subplots(figsize=(7.8, 4.8), dpi=150)
shim_series = [
    ("Balanced_Sub2p8g",    "#1b9e77", "-",  "Bare Lathe Contact (15 um)"),
    ("Balanced_Shim_65um",  "#377eb8", "--", "+ 2 Layers Kapton Shim (65 um total)"),
    ("Balanced_Shim_140um", "#ff7f00", "-.", "+ 5 Layers Kapton Shim (140 um total)"),
]
for lbl, col, ls, title_lbl in shim_series:
    sub = [r for r in t_rows if r["design_label"] == lbl]
    if not sub:
        continue
    sub.sort(key=lambda r: r["v0_v"])
    ax_s.plot(
        [r["v0_v"] for r in sub],
        [r["sep_42ms_mm"] for r in sub],
        marker="o",
        lw=2.0,
        color=col,
        ls=ls,
        label=f"{title_lbl} — Passive Hold = {sub[0]['hold_gf']:.0f} gf",
    )

ax_s.axhline(1.5, color="gray", ls=":", lw=1.1)
ax_s.set_xlabel("Trimpot Capacitor Voltage V0 (V)")
ax_s.set_ylabel("Separation @ 42 ms (mm)")
ax_s.set_title("Bench Shimming Proof: Bare Metal Works & Adding Kapton Widens Margin Further!")
ax_s.grid(True, alpha=0.3)
ax_s.legend(fontsize=8.5, loc="best")
fig_s.tight_layout()
fig_s.savefig(out_dir / "fig_06_kapton_shim_verification.png")
plt.close(fig_s)


# --- Plot 3: Render 2D Field Maps & 18V Transient Waveform of Winning 34-AWG Bare-Metal Champion ---

if best_champ_id is not None:
    ch = [r for r in s_rows if int(r["id"]) == best_champ_id][0]
    print("\n=== Winning 34-AWG Bare-Metal Design Specs ===")
    print(
        f"  ID={int(ch['id'])} | Dropper Shank D={2*ch['r_shank_mm']:.2f}mm (r={ch['r_shank_mm']}mm), "
        f"Length={ch['core_len_mm']}mm, Tip D={2*ch['r_tip_mm']:.2f}mm (r={ch['r_tip_mm']}mm)\n"
        f"  34 AWG Coil: w={ch['coil_w_mm']}mm (OD={ch['coil_od_mm']}mm), N={int(ch['n_turns_34awg'])} turns, "
        f"R_dc={ch['r_coil_ohm']:.2f} Ω (tau_RC={ch['tau_rc_ms']:.1f} ms @ 1000 uF)\n"
        f"  Zippy Target: Ø6.35mm x {ch['cap_h_mm']}mm Iron Cap + Ø6.35mm x 2.54mm N52 Magnet\n"
        f"  Masses: Drone Dropper Active = {ch['m_dropper_g']:.2f} g (Core={ch['m_core_g']:.2f}g + Cu={ch['m_cu_g']:.2f}g, NO drone magnet!) | "
        f"Zippy Target = {ch['m_zippy_g']:.2f} g\n"
        f"  Bare-Metal Hold (@ 15um) = {ch['f_hold_gf']:.1f} gf"
    )

    mdl_w, _ = build_design(ch["r_shank_mm"], ch["core_len_mm"], ch["r_tip_mm"], ch["coil_w_mm"], ch["cap_h_mm"])
    sess_w = StaticSolverSession(mdl_w, fe_order=2)
    res_h = sess_w.solve(0.0, compute_diagnostics=True)

    i_18v = 0.84 * (18.0 / (ch["r_coil_ohm"] + 0.152))
    res_r = sess_w.solve(i_18v, compute_diagnostics=True)

    plot_field_map(
        mdl_w,
        res_h,
        str(out_dir / "fig_06_winner_hold_0A.png"),
        r_max_plot_mm=8.0,
        z_min_plot_mm=-6.0,
        z_max_plot_mm=ch["core_len_mm"] + 2.0,
        title=f"Exp 06 Winner (34 AWG, m={ch['m_dropper_g']:.2f}g) — BARE HOLD: {res_h.force_z_n/9.80665*1e3:+.1f} gf",
    )
    plot_field_map(
        mdl_w,
        res_r,
        str(out_dir / "fig_06_winner_repel_18V.png"),
        r_max_plot_mm=8.0,
        z_min_plot_mm=-6.0,
        z_max_plot_mm=ch["core_len_mm"] + 2.0,
        title=f"Exp 06 Winner — 18V ACTIVE REPULSION (I={i_18v:.2f}A): {res_r.force_z_n/9.80665*1e3:+.1f} gf",
    )

    ckt_18 = CircuitParams(v0_cap_v=18.0, c1_f=cfg.CAP_C1_F, pulse_ms=cfg.PULSE_MS)
    tr_18 = solve_transient_pulse(mdl_w, ckt_18, t_end_ms=42.0, dt_ms=0.2, zippy_mass_g=27.0)
    plot_transient_summary(
        tr_18,
        str(out_dir / "fig_06_winner_transient_18V.png"),
        title=f"Exp 06 Winner (34 AWG, m={ch['m_dropper_g']:.2f}g, 1000 uF @ 18V): Bare-Metal Launch!",
    )

print(f"All Exp 06 plots saved to {out_dir}")
