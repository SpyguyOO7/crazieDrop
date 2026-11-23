import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import ngsolve as ngs

from solvers import CircuitParams, StaticSolverSession, solve_transient_pulse
from plotting import plot_field_map, plot_transient_summary
from run_exp08 import build_tip_candidate
import config as cfg

ngs.SetNumThreads(1)

out_dir = Path(__file__).resolve().parent
s_csv = out_dir / "static_tip_curves.csv"
t_csv = out_dir / "transient_tip_comparison.csv"

if not s_csv.exists() or not t_csv.exists():
    raise FileNotFoundError("Run run_exp08.py first!")

s_rows = []
with s_csv.open("r", newline="") as f:
    for r in csv.DictReader(f):
        s_rows.append({k: (v if k in ("name", "family") else float(v)) for k, v in r.items()})

t_rows = []
with t_csv.open("r", newline="") as f:
    for r in csv.DictReader(f):
        t_rows.append({k: (v if k in ("name", "family") else float(v)) for k, v in r.items()})


# --- Plot 1: Contact Area Trade-Off (Center Horn vs Flat vs Outer Rim at Matched Contact Areas!) ---

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.4), dpi=150)

# Extract summary per candidate at I = 0 A (Hold) and I = 2.4 A (~18V peak current in 34 AWG)
names = [c[0] for c in cfg.TIP_CANDIDATES]
fam_summary = {"horn": [], "flat": [], "rim": []}
for nm in names:
    sub = [r for r in s_rows if r["name"] == nm]
    if not sub:
        continue
    r0 = sub[0]
    f_hold = r0["force_gf"]
    f_rep_2p4 = [r["force_gf"] for r in sub if abs(r["current_a"] - 2.4) < 1e-3][0]
    f_rep_max = min(r["force_gf"] for r in sub)
    fam_summary[r0["family"]].append((nm, r0["contact_area_mm2"], f_hold, -f_rep_2p4, -f_rep_max))

# Sort standard taper=1.5mm horns and cup=0.8mm rims by contact area for clean curves
horns_std = [x for x in fam_summary["horn"] if "t1p5" in x[0]]
horns_std.sort(key=lambda x: x[1])
rims_std = [x for x in fam_summary["rim"] if "cup0p8" in x[0]]
rims_std.sort(key=lambda x: x[1])
flats = sorted(fam_summary["flat"], key=lambda x: x[1])

ax1.plot([x[1] for x in horns_std], [x[2] for x in horns_std], "o-", color="#1b9e77", lw=2.1, label="Center Horn (varying tip D)")
ax1.plot([x[1] for x in rims_std],  [x[2] for x in rims_std],  "s-", color="#d95f02", lw=2.1, label="Outer Rim Cup (varying rim width)")
ax1.plot([x[1] for x in flats],     [x[2] for x in flats],     "^-", color="#7570b3", lw=2.1, ms=8, label="No-Horn Flat Cylinder")
ax1.set_xlabel("Pole Contact Area A_contact (mm^2)")
ax1.set_ylabel("Passive Holding Force @ I=0A (gf)")
ax1.set_title("Passive Hold vs Contact Area: Horn vs Rim vs Flat")
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=8.5)

ax2.plot([x[1] for x in horns_std], [x[3] for x in horns_std], "o-", color="#1b9e77", lw=2.1, label="Center Horn @ I=2.4A (~18V)")
ax2.plot([x[1] for x in rims_std],  [x[3] for x in rims_std],  "s-", color="#d95f02", lw=2.1, label="Outer Rim Cup @ I=2.4A (~18V)")
ax2.plot([x[1] for x in flats],     [x[3] for x in flats],     "^-", color="#7570b3", lw=2.1, ms=8, label="No-Horn Flat @ I=2.4A (~18V)")
ax2.axhline(0.0, color="black", ls="--", lw=1.0)
ax2.set_xlabel("Pole Contact Area A_contact (mm^2)")
ax2.set_ylabel("Active Downward Repulsion -F_z @ 2.4A (gf)")
ax2.set_title("Active Repulsive Shove vs Contact Area")
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=8.5)

fig.tight_layout()
fig.savefig(out_dir / "fig_08_area_horn_vs_rim_vs_flat.png")
plt.close(fig)


# --- Plot 2: Static F_z(I) Curves & Moving-Gap Transient Trimpot Launch Window ---

fig2, (ax_cur, ax_v) = plt.subplots(2, 1, figsize=(9.4, 8.2), dpi=150)

arch_styles = {
    "horn_d2p2_t1p5":    ("#66a61e", "--", "Narrow Horn (d=2.2mm)"),
    "horn_d3p4_t1p5":    ("#1b9e77", "-",  "Medium Horn (d=3.4mm)"),
    "horn_d3p8_t1p5":    ("#2b8cbe", "-",  "Wide Horn (d=3.8mm, Exp06 Winner)"),
    "flat_no_horn_d4p4": ("#7570b3", "-.", "No-Horn Flat Cylinder (d=4.4mm)"),
    "rim_w0p35_cup0p8":  ("#e7298a", ":",  "Thin Outer Rim (w=0.35mm)"),
    "rim_w0p60_cup0p8":  ("#d95f02", "-",  "Medium Outer Rim (w=0.60mm)"),
    "rim_w0p90_cup0p8":  ("#a65628", "--", "Wide Outer Rim (w=0.90mm)"),
}

for arch, (col, ls, lbl) in arch_styles.items():
    sub_s = [r for r in s_rows if r["name"] == arch]
    if sub_s:
        ax_cur.plot([r["current_a"] for r in sub_s], [r["force_gf"] for r in sub_s], color=col, ls=ls, lw=1.9, label=lbl)

    sub_t = [r for r in t_rows if r["name"] == arch]
    if sub_t:
        sub_t.sort(key=lambda r: r["v0_v"])
        n_ok = sum(int(r["dropped"]) for r in sub_t)
        hold0 = sub_t[0]["hold_gf"]
        ax_v.plot(
            [r["v0_v"] for r in sub_t],
            [r["sep_42ms_mm"] for r in sub_t],
            marker="o",
            color=col,
            ls=ls,
            lw=2.0,
            label=f"{lbl} (Hold={hold0:.0f}gf) -> {n_ok}/{len(sub_t)} Launch",
        )

ax_cur.axhline(0.0, color="black", ls="--", lw=1.1, label="0 gf (Negative = Repulsion)")
ax_cur.set_xlabel("34 AWG Coil Current I (A)")
ax_cur.set_ylabel("Force on Zippy F_z (gf)")
ax_cur.set_title("Static Push-Off Curves: Center Horn vs Flat vs Outer Rim")
ax_cur.grid(True, alpha=0.3)
ax_cur.legend(fontsize=7.8, loc="upper right")

ax_v.axhline(1.5, color="gray", ls=":", lw=1.1)
ax_v.set_xlabel("Trimpot Capacitor Voltage V0 (V)  [1000 uF, 34 AWG]")
ax_v.set_ylabel("Separation @ 42 ms (mm)")
ax_v.set_title("Moving-Gap Transient Launch Window: Horn vs Flat vs Outer Rim!")
ax_v.grid(True, alpha=0.3)
ax_v.legend(fontsize=7.8, loc="upper left")

fig2.tight_layout()
fig2.savefig(out_dir / "fig_08_transient_launch_comparison.png")
plt.close(fig2)


# --- Plot 3: Render 2D Field Cross-Sections of No-Horn Flat & Best Outer-Rim Contact (@ Hold & 18V Repel) ---

for showcase_name, tag in [("flat_no_horn_d4p4", "flat_no_horn"), ("rim_w0p60_cup0p8", "outer_rim_w0p60")]:
    mdl, _, _ = build_tip_candidate(showcase_name)
    sess = StaticSolverSession(mdl, fe_order=2)
    r_c = mdl.parts[1].coil_res_ohm
    i_18v = 0.84 * (18.0 / (r_c + 0.152))

    res_0 = sess.solve(0.0, compute_diagnostics=True)
    res_18 = sess.solve(i_18v, compute_diagnostics=True)

    plot_field_map(
        mdl,
        res_0,
        str(out_dir / f"fig_08_{tag}_hold_0A.png"),
        r_max_plot_mm=8.0,
        z_min_plot_mm=-5.5,
        z_max_plot_mm=10.0,
        title=f"{showcase_name} — PASSIVE HOLD (I=0A): F_z = {res_0.force_z_n/9.80665*1e3:+.1f} gf",
    )
    plot_field_map(
        mdl,
        res_18,
        str(out_dir / f"fig_08_{tag}_repel_18V.png"),
        r_max_plot_mm=8.0,
        z_min_plot_mm=-5.5,
        z_max_plot_mm=10.0,
        title=f"{showcase_name} — 18V ACTIVE REPULSION (I={i_18v:.2f}A): F_z = {res_18.force_z_n/9.80665*1e3:+.1f} gf",
    )

print(f"All Exp 08 plots saved to {out_dir}")
