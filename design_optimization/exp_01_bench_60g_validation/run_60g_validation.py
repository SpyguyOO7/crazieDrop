from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from materials import STEEL_1018, MAG_IRON_4000, LATHE_FINISH_GAP_MM, Material, make_ndfeb, make_wire_spec
from geometry import AxiPart, COIL_MAT, rect_rz, build_axi_model
from solvers import CircuitParams, solve_static, solve_transient_pulse
from plotting import plot_field_map, plot_transient_summary


# --- Parameters & Output Setup ---

out_dir = Path(__file__).resolve().parent

r_rod_mm = 6.35 / 2.0               # 1/4" (6.35 mm) McMaster magnetic iron rod -> r = 3.175 mm
core_len_mm = 9.5                   # 9.5 mm dropper core length
pm_h_real_mm = 2.54                 # 1/10" (2.54 mm) McMaster N52 magnet height
lathe_gap_mm = LATHE_FINISH_GAP_MM  # 15 um careful lathe-cut contact finish baseline

# 60g leftover piece of McMaster 1/4" magnetic iron rod (~241 mm physical length; magnetic return asymptotes by ~35 mm)
rod_60g_len_mm = 35.0

# Range of lathe surface finishes: fine polished lathe (6 um) -> careful lathe (15 um) -> rough / slight dust/tilt (50-200 um)
gaps_mm = np.array([0.006, 0.010, 0.015, 0.025, 0.050, 0.100, 0.200, 0.400])


# --- Helper to Build Dropper + Target ---

def build_case(
    pm_h_mm: float,
    pm_grade: str,
    core_mat: Material,
    gap_mm: float,
    target_pts_rz: tuple[tuple[float, float], ...],
    include_coil: bool = False,
):
    parts = [
        rect_rz("core_rod", 0.0, r_rod_mm, 0.0, core_len_mm, core_mat, maxh_mm=0.22, group="dropper"),
        rect_rz(
            "top_pm",
            0.0,
            r_rod_mm,
            core_len_mm,
            core_len_mm + pm_h_mm,
            make_ndfeb(pm_grade, sign=1.0),
            maxh_mm=0.18,
            group="dropper",
        ),
        AxiPart("target_iron", target_pts_rz, core_mat, maxh_mm=0.20, group="target"),
    ]
    if include_coil:
        wire_32 = make_wire_spec(awg=32, fill_factor=0.58)
        parts.append(
            rect_rz(
                "choke_coil",
                r_rod_mm + 0.35,
                r_rod_mm + 0.35 + 1.8,
                0.3,
                9.2,
                COIL_MAT,
                maxh_mm=0.32,
                group="dropper",
                wire=wire_32,
                coil_sign=-1.0,
            )
        )
    return build_axi_model(parts, air_pad_r_mm=30.0, air_pad_z_mm=35.0, air_maxh_mm=2.5)


def rect_target_pts(r_in: float, r_out: float, length_mm: float, gap_mm: float):
    return (
        (r_in, -length_mm - gap_mm),
        (r_out, -length_mm - gap_mm),
        (r_out, -gap_mm),
        (r_in, -gap_mm),
    )


# --- Experiment 1A: Holding Force Across Lathe Surface Finish, Material (1018 vs Mag Iron 4000) & Target Size ---

print("=== Experiment 1A: Passive Holding Force vs Lathe Finish Gap, Material & Target Length ===")

configs = [
    ("Old Assumed (1.5mm N42 + 1018 Steel + 3mm Slug)", 1.50, "N42", STEEL_1018, 3.0),
    ("McMaster N52 (0.1\") + 1018 Steel (mu=850) + 3mm Zippy", 2.54, "N52", STEEL_1018, 3.0),
    ("McMaster N52 (0.1\") + Mag Iron (mu=4000) + 3mm Zippy", 2.54, "N52", MAG_IRON_4000, 3.0),
    ("McMaster N52 (0.1\") + Mag Iron (mu=4000) + 60g Bench Rod", 2.54, "N52", MAG_IRON_4000, rod_60g_len_mm),
]

force_curves = {}
for label, pm_h, grade, c_mat, tgt_len in configs:
    f_vals = []
    for g in gaps_mm:
        mdl = build_case(pm_h, grade, c_mat, float(g), rect_target_pts(0.0, r_rod_mm, tgt_len, float(g)))
        res = solve_static(mdl, coil_current_a=0.0, fe_order=2)
        f_vals.append(res.force_z_n)
    force_curves[label] = np.array(f_vals)
    print(f"\n{label}:")
    for g, fv in zip(gaps_mm, f_vals):
        grams = fv / 9.80665 * 1000.0
        tag = " <-- Nominal 15um Lathe Cut" if abs(g - lathe_gap_mm) < 1e-6 else ""
        print(f"  gap = {g * 1000:5.0f} um ({g:5.3f} mm) -> F_z = {fv:6.3f} N ({grams:5.1f} gf){tag}")

fig, ax = plt.subplots(figsize=(8.4, 5.4), dpi=150)
for label, f_vals in force_curves.items():
    ax.plot(gaps_mm * 1e3, f_vals / 9.80665 * 1e3, marker="o", lw=1.9, label=label)

ax.axvline(lathe_gap_mm * 1e3, color="gray", ls="-.", lw=1.2, label="Careful Lathe Finish (~15 um)")
ax.axhline(60.0, color="red", ls="--", lw=1.3, label="60 g Bench Rod Weight")
ax.axhline(27.0, color="darkorange", ls=":", lw=1.4, label="27 g Zippy Weight")
ax.set_xlabel("Effective Contact Air Gap / Lathe Roughness (um)")
ax.set_ylabel("Passive Holding Force (grams-force)")
ax.set_title("Passive Hold Force: McMaster N52 + Mag Iron (mu_r=4000) vs Lathe Finish")
ax.set_xscale("log")
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=8.2)
fig.tight_layout()
fig.savefig(out_dir / "fig_hold_vs_gap_and_target.png")
plt.close(fig)


# --- Experiment 1B: Rim Flux Concentration & Flush vs Tilted-Edge Contact (@ 15um Lathe Finish) ---

print("\n=== Experiment 1B: Rim Flux Concentration & Edge/Tip Contact (@ 15 um Lathe Finish) ===")

mdl_flush = build_case(
    pm_h_real_mm, "N52", MAG_IRON_4000, lathe_gap_mm, rect_target_pts(0.0, r_rod_mm, rod_60g_len_mm, lathe_gap_mm)
)
res_flush = solve_static(mdl_flush, coil_current_a=0.0, fe_order=2)
plot_field_map(
    mdl_flush,
    res_flush,
    str(out_dir / "fig_60g_flush_contact.png"),
    z_min_plot_mm=-12.0,
    z_max_plot_mm=13.5,
    title=f"McMaster N52 + Mag Iron 60g Rod (Flush @ 15um): F_z = {res_flush.force_z_n:.3f} N ({res_flush.force_z_n/9.80665*1e3:.1f} gf)",
)

# Shallow tilt (~1 deg across face -> 15um at outer rim, ~100um at center)
shallow_tilt_pts = (
    (0.0, -rod_60g_len_mm - lathe_gap_mm),
    (r_rod_mm + 0.8, -rod_60g_len_mm - lathe_gap_mm),
    (r_rod_mm + 0.8, -0.09),
    (r_rod_mm, -lathe_gap_mm),
    (r_rod_mm - 1.1, -lathe_gap_mm),
    (0.0, -0.10),
)
mdl_shallow_tilt = build_case(pm_h_real_mm, "N52", MAG_IRON_4000, lathe_gap_mm, shallow_tilt_pts)
res_shallow_tilt = solve_static(mdl_shallow_tilt, coil_current_a=0.0, fe_order=2)
plot_field_map(
    mdl_shallow_tilt,
    res_shallow_tilt,
    str(out_dir / "fig_60g_rim_tilted_contact.png"),
    z_min_plot_mm=-12.0,
    z_max_plot_mm=13.5,
    title=f"Slight Rim-Edge Contact (~1 deg tilt @ 15um rim): F_z = {res_shallow_tilt.force_z_n:.3f} N ({res_shallow_tilt.force_z_n/9.80665*1e3:.1f} gf)",
)

# Concentrated small tip contact (r <= 1.5 mm at 15 um)
tip_pts = (
    (0.0, -rod_60g_len_mm - lathe_gap_mm),
    (r_rod_mm, -rod_60g_len_mm - lathe_gap_mm),
    (r_rod_mm, -0.50),
    (1.5, -lathe_gap_mm),
    (0.0, -lathe_gap_mm),
)
mdl_tip = build_case(pm_h_real_mm, "N52", MAG_IRON_4000, lathe_gap_mm, tip_pts)
res_tip = solve_static(mdl_tip, coil_current_a=0.0, fe_order=2)
plot_field_map(
    mdl_tip,
    res_tip,
    str(out_dir / "fig_60g_concentrated_tip_contact.png"),
    z_min_plot_mm=-12.0,
    z_max_plot_mm=13.5,
    title=f"Concentrated Tip Contact (r<=1.5mm @ 15um): F_z = {res_tip.force_z_n:.3f} N ({res_tip.force_z_n/9.80665*1e3:.1f} gf)",
)

print(f"  Flush 6.35mm lathe face (15um):      F_z = {res_flush.force_z_n:.3f} N ({res_flush.force_z_n/9.80665*1e3:.1f} gf)")
print(f"  Slight rim-edge contact (~1 deg):    F_z = {res_shallow_tilt.force_z_n:.3f} N ({res_shallow_tilt.force_z_n/9.80665*1e3:.1f} gf)")
print(f"  Concentrated tip funnel (r<=1.5mm):  F_z = {res_tip.force_z_n:.3f} N ({res_tip.force_z_n/9.80665*1e3:.1f} gf)")


# --- Experiment 1C: 20V vs Tuned-V0 Transient Discharge (with mu_r=4000, sigma=1e7 S/m Solid Iron Eddy Currents!) ---

print("\n=== Experiment 1C: Transient Capacitor Discharge on McMaster Mag Iron (mu_r=4000, 15um Lathe Finish) ===")
mdl_zippy = build_case(
    pm_h_real_mm,
    "N52",
    MAG_IRON_4000,
    lathe_gap_mm,
    rect_target_pts(0.0, r_rod_mm, 3.0, lathe_gap_mm),
    include_coil=True,
)

for v0_test in (20.0, 9.0):
    ckt = CircuitParams(v0_cap_v=v0_test, c1_f=1000e-6, pulse_ms=25.0)
    tr = solve_transient_pulse(mdl_zippy, ckt, t_end_ms=45.0, dt_ms=0.2, zippy_mass_g=27.0)
    out_file = out_dir / f"fig_transient_{int(v0_test)}V.png"
    plot_transient_summary(
        tr,
        str(out_file),
        title=f"McMaster N52 + Mag Iron (mu_r=4000) @ 15um Lathe Gap: V0 = {v0_test:.0f} V, C1 = 1000 uF",
    )
    print(
        f"  V0 = {v0_test:4.1f} V -> Peak I = {tr.i_coil_a.max():.2f} A | "
        f"Init Hold = {tr.force_z_n[0]:.3f} N ({tr.force_z_n[0]/9.80665*1e3:.1f} gf) | "
        f"Peak Overdrive F_z = {tr.force_z_n[1:80].max():.3f} N ({tr.force_z_n[1:80].max()/9.80665*1e3:.1f} gf) | "
        f"Min F_z = {tr.force_z_n.min():.4f} N | Sep @ 45ms = {tr.zippy_disp_mm[-1]:.2f} mm | Dropped = {tr.dropped}"
    )
