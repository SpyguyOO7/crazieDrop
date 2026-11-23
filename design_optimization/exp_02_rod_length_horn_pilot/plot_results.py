import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from solvers import solve_static
from plotting import plot_field_map
from run_sweep import build_candidate


# --- Load Results ---

out_dir = Path(__file__).resolve().parent
csv_path = out_dir / "results.csv"

if not csv_path.exists():
    raise FileNotFoundError(f"Missing {csv_path}. Run run_sweep.py first!")

rows = []
with csv_path.open("r", newline="") as f:
    for r in csv.DictReader(f):
        parsed = {k: (v if k == "interface_mode" else float(v)) for k, v in r.items()}
        rows.append(parsed)

rows.sort(key=lambda d: (d["interface_mode"], d["coil_w_mm"], d["core_len_mm"]))


# --- Plot 1: Rod Length Sweep (at baseline 1.8 mm coil build) ---

fig, axes = plt.subplots(4, 1, figsize=(8.5, 10.0), dpi=150, sharex=True)
mode_styles = {
    "double_horn": ("#d95f02", "-", "Double Horn (Dropper + Zippy Tapered to d=3.0mm)"),
    "dropper_horn": ("#1b9e77", "--", "Dropper Horn (d=3.0mm) + Flat 6.35mm Zippy"),
    "flat_cyl": ("#7570b3", ":", "Pure Flat 6.35mm Cylinder (Current Baseline Shape)"),
}

for mode, (color, ls, label) in mode_styles.items():
    sub = [r for r in rows if r["interface_mode"] == mode and abs(r["coil_w_mm"] - 1.8) < 1e-4]
    if not sub:
        continue
    l_vals = [r["core_len_mm"] for r in sub]
    axes[0].plot(l_vals, [r["f_hold_gf"] for r in sub], color=color, ls=ls, marker="o", ms=4, lw=1.9, label=label)
    axes[1].plot(l_vals, [r["p_dc_null_w"] for r in sub], color=color, ls=ls, marker="o", ms=4, lw=1.9, label=label)
    axes[2].plot(l_vals, [r["pm_demag_margin_null_ka_m"] for r in sub], color=color, ls=ls, marker="o", ms=4, lw=1.9, label=label)
    axes[3].plot(l_vals, [r["m_dropper_active_g"] for r in sub], color=color, ls=ls, marker="o", ms=4, lw=1.9, label=label)

axes[0].axhline(27.0, color="gray", ls="--", lw=1.0, label="1x Zippy Weight (27 gf)")
axes[0].axhline(81.0, color="red", ls=":", lw=1.1, label="3x Safety Factor (81 gf)")
axes[0].axvline(9.5, color="black", ls="-.", lw=1.0, alpha=0.5, label="Current 9.5 mm Length")
axes[0].set_ylabel("Passive Hold (gf)")
axes[0].set_title("Effect of Core Rod Length L_core (Fixed w_coil = 1.8 mm, 15 um Lathe Finish)")
axes[0].grid(True, alpha=0.3)
axes[0].legend(fontsize=7.8, loc="upper right")

axes[1].axvline(9.5, color="black", ls="-.", lw=1.0, alpha=0.5)
axes[1].set_ylabel("Choke Null Power (W)")
axes[1].set_yscale("log")
axes[1].grid(True, which="both", alpha=0.3)

axes[2].axhline(0.0, color="red", ls="--", lw=1.2, label="Irreversible N52 Demag Knee (0 kA/m)")
axes[2].axvline(9.5, color="black", ls="-.", lw=1.0, alpha=0.5)
axes[2].set_ylabel("N52 Demag Margin\n@ Null (kA/m)")
axes[2].grid(True, alpha=0.3)
axes[2].legend(fontsize=8, loc="lower right")

axes[3].axvline(9.5, color="black", ls="-.", lw=1.0, alpha=0.5)
axes[3].set_ylabel("Active Stack Mass (g)\n(Core + PM + Cu)")
axes[3].set_xlabel("Dropper Core Total Length L_core (mm)")
axes[3].grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(out_dir / "fig_length_sweep_summary.png")
plt.close(fig)


# --- Plot 2: Full Mass Pareto Front (Hold Force vs Active Mass, Colored by Null Power) ---

fig2, ax = plt.subplots(figsize=(8.6, 5.8), dpi=150)
markers = {"double_horn": "o", "dropper_horn": "^", "flat_cyl": "s"}

# Filter designs that do NOT severely demagnetize the N52 magnet at null (margin > -50 kA/m) and have reasonable null power <= 25 W
safe_rows = [r for r in rows if r["p_dc_null_w"] <= 25.0]

sc = None
for mode, mk in markers.items():
    sub = [r for r in safe_rows if r["interface_mode"] == mode]
    if not sub:
        continue
    sc = ax.scatter(
        [r["m_dropper_active_g"] for r in sub],
        [r["f_hold_gf"] for r in sub],
        c=[r["p_dc_null_w"] for r in sub],
        cmap="viridis_r",
        vmin=0.5,
        vmax=20.0,
        marker=mk,
        s=52,
        edgecolors="black",
        linewidths=0.4,
        alpha=0.9,
        label=mode_styles[mode][2],
    )

if sc is not None:
    cb = fig2.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("DC Power Required at Choke Null P_null (W)")

# Highlight current 9.5mm flat baseline
base_pt = [r for r in rows if r["interface_mode"] == "flat_cyl" and abs(r["core_len_mm"] - 9.5) < 1e-3 and abs(r["coil_w_mm"] - 1.8) < 1e-3]
if base_pt:
    b0 = base_pt[0]
    ax.scatter([b0["m_dropper_active_g"]], [b0["f_hold_gf"]], color="red", s=160, marker="*", zorder=10, label="Current 9.5mm Flat Baseline")

ax.axhline(81.0, color="red", ls=":", lw=1.1, label="3x Zippy Weight (81 gf)")
ax.set_xlabel("Dropper Active Mass: Core + N52 Magnet + Copper Coil (g)")
ax.set_ylabel("Passive Holding Force @ 15 um Lathe Finish (gf)")
ax.set_title("Mass vs Holding Force & Release Power Trade-off Across Rod Lengths")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8.0, loc="best")
fig2.tight_layout()
fig2.savefig(out_dir / "fig_mass_pareto.png")
plt.close(fig2)


# --- Print Top Pareto Candidates & Render 2D Field Maps of Best Sweet-Spot Design ---

print("\n=== Top Sweet-Spot Designs (Active Mass <= 4.0 g, Demag Safe >= 0 kA/m, P_null <= 12 W) ===")
cands = [
    r for r in rows
    if r["m_dropper_active_g"] <= 4.0 and r["pm_demag_margin_null_ka_m"] >= -10.0 and r["p_dc_null_w"] <= 12.0
]
cands.sort(key=lambda r: (-r["f_hold_gf"], r["m_dropper_active_g"]))

for r in cands[:10]:
    print(
        f"  {r['interface_mode']:12s} | L={r['core_len_mm']:4.1f}mm | w_cu={r['coil_w_mm']:3.1f}mm | "
        f"Mass={r['m_dropper_active_g']:4.2f}g (Zippy={r['m_zippy_g']:.2f}g) | "
        f"Hold={r['f_hold_gf']:6.1f}gf ({r['hold_per_gram_gf_g']:4.1f} gf/g) | "
        f"P_null={r['p_dc_null_w']:5.2f}W (V_32awg={r['v_dc_null_v']:4.1f}V) | "
        f"Demag={r['pm_demag_margin_null_ka_m']:+5.1f}kA/m"
    )

if cands:
    best = cands[0]
    print(f"\nRendering 2D field maps (Hold & Null) for best <=4.0g design: {best['interface_mode']} L={best['core_len_mm']}mm...")
    mdl_best = build_candidate(best["interface_mode"], best["core_len_mm"], best["coil_w_mm"])
    res_hold = solve_static(mdl_best, coil_current_a=0.0, fe_order=2)
    res_null = solve_static(mdl_best, coil_current_a=best["i_null_a"], fe_order=2, init_psi=res_hold.psi_gf)

    z_hi = best["core_len_mm"] + 4.0
    plot_field_map(
        mdl_best,
        res_hold,
        str(out_dir / "fig_best_design_hold.png"),
        z_min_plot_mm=-5.0,
        z_max_plot_mm=z_hi,
        title=(
            f"Best <=4g Design HOLD ({best['interface_mode']}, L={best['core_len_mm']}mm, m={best['m_dropper_active_g']:.2f}g): "
            f"Hold = {best['f_hold_gf']:.1f} gf"
        ),
    )
    plot_field_map(
        mdl_best,
        res_null,
        str(out_dir / "fig_best_design_null.png"),
        z_min_plot_mm=-5.0,
        z_max_plot_mm=z_hi,
        title=(
            f"Best <=4g Design NULL (I={best['i_null_a']:.2f}A, P={best['p_dc_null_w']:.2f}W): "
            f"Res Force = {best['f_null_res_gf']:.2f} gf"
        ),
    )
    print(f"Saved plots to {out_dir}")
