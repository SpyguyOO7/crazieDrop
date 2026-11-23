import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import ngsolve as ngs

from materials import MAG_IRON_4000, LATHE_FINISH_GAP_MM, make_ndfeb, make_wire_spec
from geometry import AxiPart, COIL_MAT, rect_rz, build_axi_model
from solvers import CircuitParams, StaticSolverSession, solve_transient_pulse
from plotting import plot_field_map, plot_transient_summary

ngs.SetNumThreads(1)


# --- Experiment 03 Parameters: Sweet-Spot Dropper Horn + Real Hand-Windable Wire Gauge Check ---

out_dir = Path(__file__).resolve().parent

r_stock_mm = 6.35 / 2.0               # 1/4" (6.35 mm) McMaster Magnetic Iron (mu_r = 4000)
pm_h_mm = 2.54                        # 1/10" (2.54 mm) McMaster Grade N52 disk magnet
gap_mm = LATHE_FINISH_GAP_MM          # 15 um careful lathe finish

# Candidate sweet-spot core geometry from Exp 02:
#   L_core = 7.0 mm, bottom 1.5 mm conical lathe taper to r_tip = 1.50 mm (3.0 mm contact land),
#   mating against a flat 6.35 mm x 3.0 mm Zippy iron target (allows +/-1.6 mm lateral landing tolerance!).
core_len_mm = 7.0
taper_len_mm = 1.5
r_tip_mm = 1.50

bobbin_wall_mm = 0.35
coil_z_bot_mm = bobbin_wall_mm
coil_z_top_mm = core_len_mm - bobbin_wall_mm   # 6.30 mm axial winding slot
coil_w_mm = 1.80                               # 1.80 mm radial winding build (OD = 10.65 mm, well inside 18x22mm PCB!)
coil_r_in_mm = r_stock_mm + bobbin_wall_mm
coil_r_out_mm = coil_r_in_mm + coil_w_mm

# Realistic hand/lathe layer-winding packing check across standard stocked magnet wire gauges (28 .. 36 AWG):
# Exact physical layer/turn count inside [6.30 mm axial x 1.80 mm radial] window using insulated wire diameter d_ins_mm!
CANDIDATE_AWGS = [28, 30, 32, 34, 36]


# --- 1. Physical Winding Fit & PCB Impedance Matching Table ---

print("=== 1. Real Physical Winding Fit & Driver-Board Matching (Window: 6.30 mm L x 1.80 mm W) ===")
win_len_mm = coil_z_top_mm - coil_z_bot_mm

core_pts = (
    (0.0, 0.0),
    (r_tip_mm, 0.0),
    (r_stock_mm, taper_len_mm),
    (r_stock_mm, core_len_mm),
    (0.0, core_len_mm),
)
tgt_pts = (
    (0.0, -3.0 - gap_mm),
    (r_stock_mm, -3.0 - gap_mm),
    (r_stock_mm, -gap_mm),
    (0.0, -gap_mm),
)


def build_horn_with_awg(awg: int):
    # Compute integer layers & integer turns/layer using realistic orthocyclic/neat-hand layer factor (0.92 * floor)
    w_spec_raw = make_wire_spec(awg=awg, fill_factor=0.58, ins_ratio=1.09)
    n_layers = int(np.floor(coil_w_mm / w_spec_raw.d_ins_mm))
    turns_per_layer = int(np.floor(win_len_mm / w_spec_raw.d_ins_mm))
    # Practical hand-wound neatness factor (~0.88 of ideal square-grid count to account for layer crossovers)
    n_real_turns = int(np.round(0.88 * n_layers * turns_per_layer))
    eff_fill = (n_real_turns * w_spec_raw.area_cu_mm2) / (coil_w_mm * win_len_mm)
    w_spec = make_wire_spec(awg=awg, fill_factor=eff_fill, ins_ratio=1.09)

    parts = [
        AxiPart("core_rod", core_pts, MAG_IRON_4000, maxh_mm=0.20, group="dropper"),
        rect_rz("top_pm", 0.0, r_stock_mm, core_len_mm, core_len_mm + pm_h_mm, make_ndfeb("N52", 1.0), maxh_mm=0.18, group="dropper"),
        rect_rz(
            "choke_coil",
            coil_r_in_mm,
            coil_r_out_mm,
            coil_z_bot_mm,
            coil_z_top_mm,
            COIL_MAT,
            maxh_mm=0.28,
            group="dropper",
            wire=w_spec,
            coil_sign=-1.0,
        ),
        AxiPart("zippy_iron", tgt_pts, MAG_IRON_4000, maxh_mm=0.20, group="target"),
    ]
    mdl = build_axi_model(parts, air_pad_r_mm=28.0, air_pad_z_mm=28.0, air_maxh_mm=2.5)
    return mdl, w_spec, n_layers, turns_per_layer, n_real_turns


# First solve static MMF requirement once using 32 AWG reference
mdl_32, w32, _, _, _ = build_horn_with_awg(32)
sess_32 = StaticSolverSession(mdl_32, fe_order=2)
res_hold = sess_32.solve(0.0, compute_diagnostics=True)

# Find exact static choke null NI_null
n32 = mdl_32.parts[2].n_turns
i_scan = np.linspace(0.5, 4.5, 17)
best_null_res = None
for amp in i_scan:
    r_tmp = sess_32.solve(float(amp), compute_diagnostics=False)
    if best_null_res is None or abs(r_tmp.force_z_n) < abs(best_null_res.force_z_n):
        best_null_res = r_tmp

# Fine polish around best_null_res
assert best_null_res is not None
for amp in np.linspace(best_null_res.coil_current_a - 0.25, best_null_res.coil_current_a + 0.25, 11):
    r_tmp = sess_32.solve(float(amp), compute_diagnostics=True)
    if abs(r_tmp.force_z_n) < abs(best_null_res.force_z_n):
        best_null_res = r_tmp

ni_null_req = best_null_res.coil_current_a * n32
print(f"\nCandidate Core Hold Force (I=0): {res_hold.force_z_n:.3f} N ({res_hold.force_z_n/9.80665*1e3:.1f} gf)")
print(f"Static Choke Null Requirement:   NI_null = {ni_null_req:.1f} Amp-turns (Residual F_z = {best_null_res.force_z_n/9.80665*1e3:.2f} gf)\n")

print(
    f"{'AWG':>3s} | {'d_bare':>7s} | {'Layers x T/L':>12s} | {'N_real':>6s} | {'m_Cu':>6s} | "
    f"{'R_coil':>7s} | {'I_null':>6s} | {'V_dc_null':>9s} | {'tau_RC (1mF)':>12s} | {'Hand-Windability / Verdict'}"
)
print("-" * 116)

gauge_models = {}
for awg in CANDIDATE_AWGS:
    mdl_g, ws, n_lay, tpl, n_real = build_horn_with_awg(awg)
    coil_p = mdl_g.parts[2]
    r_c = coil_p.coil_res_ohm
    i_n = ni_null_req / coil_p.n_turns
    v_n = i_n * r_c + i_n * (0.12 + 0.032)  # Include C1 ESR (0.12) + AO3400A Rdson (0.032)
    tau_rc_ms = (r_c + 0.152) * 1000e-6 * 1e3
    gauge_models[awg] = (mdl_g, v_n, i_n)

    if awg == 28:
        note = "Thick/stiff, fast 1.1ms C1 dump (< eddy time!), high ESR loss"
    elif awg == 30:
        note = "Easy wind, good current, fast 1.7ms C1 droop"
    elif awg == 32:
        note = "SWEET SPOT: Strong wire (0.20mm), ~178 turns, V0~9-11V, 2.7ms RC!"
    elif awg == 34:
        note = "GREAT ALT: Easy wind (0.16mm), ~277 turns, V0~13-16V, 6.2ms flat RC!"
    else:
        note = "Fragile hair-thin wire (0.127mm), needs ~21-25V near 28V boost limit"

    print(
        f"{awg:3d} | {ws.d_bare_mm:5.3f}mm | {n_lay:2d} x {tpl:2d}      | {n_real:6d} | {coil_p.mass_g:5.2f}g | "
        f"{r_c:5.2f} Ω | {i_n:5.2f}A | {v_n:7.2f} V | {tau_rc_ms:9.1f} ms | {note}"
    )


# --- 2. Render 2D Magnetic Field Cross-Sections for Selected 32 AWG Design (Matching Exp 01 Style) ---

print("\n=== 2. Rendering 2D Magnetic Field Maps for 32 AWG Sweet-Spot Design ===")
mdl_sel, v_dc_32, i_null_32 = gauge_models[32]
masses_32 = mdl_sel.mass_breakdown_g("dropper")
print(f"32 AWG Active Dropper Mass Breakdown (g): { {k: round(v, 3) for k, v in masses_32.items()} }")
print(f"Zippy Flat Iron Target Mass (g):          {mdl_sel.parts[3].mass_g:.3f} g")

sess_sel = StaticSolverSession(mdl_sel, fe_order=2)
res_sel_hold = sess_sel.solve(0.0, compute_diagnostics=True)
res_sel_null = sess_sel.solve(i_null_32, compute_diagnostics=True)
res_sel_over = sess_sel.solve(4.2, compute_diagnostics=True)  # ~Peak current if over-driven at 20V

plot_field_map(
    mdl_sel,
    res_sel_hold,
    str(out_dir / "fig_03_field_passive_hold.png"),
    r_max_plot_mm=8.5,
    z_min_plot_mm=-4.5,
    z_max_plot_mm=11.5,
    title=(
        f"Exp 03 Horn (L=7.0mm, 32 AWG, m={masses_32['total']:.2f}g) — PASSIVE HOLD (I=0 A): "
        f"F_z = {res_sel_hold.force_z_n:.3f} N ({res_sel_hold.force_z_n/9.80665*1e3:.1f} gf)"
    ),
)

plot_field_map(
    mdl_sel,
    res_sel_null,
    str(out_dir / "fig_03_field_choke_null.png"),
    r_max_plot_mm=8.5,
    z_min_plot_mm=-4.5,
    z_max_plot_mm=11.5,
    title=(
        f"Exp 03 Horn — ACTIVE CHOKE NULL (I={i_null_32:.2f} A, NI={ni_null_req:.0f} A-t): "
        f"F_z = {res_sel_null.force_z_n:+.4f} N ({res_sel_null.force_z_n/9.80665*1e3:+.2f} gf)"
    ),
)

plot_field_map(
    mdl_sel,
    res_sel_over,
    str(out_dir / "fig_03_field_overdrive_4p2A.png"),
    r_max_plot_mm=8.5,
    z_min_plot_mm=-4.5,
    z_max_plot_mm=11.5,
    title=(
        f"Exp 03 Horn — OVERDRIVEN @ 4.2 A (Untuned 20V Spike): "
        f"Reverse Re-Attraction F_z = {res_sel_over.force_z_n:+.3f} N ({res_sel_over.force_z_n/9.80665*1e3:.1f} gf)"
    ),
)


# --- 3. Full Driver-Board Capacitor + Eddy-Current + Moving-Gap Drop Dynamics ---
# Compare:
#   A) 32 AWG (your current wire!) at Tuned V0 = 11.5 V vs Untuned V0 = 20.0 V (1000 uF C1, 25 ms pulse)
#   B) 34 AWG (slightly finer, still easy to wind: 0.16 mm bare) at Tuned V0 = 16.5 V (slower 6.2 ms RC flat-top!)
#   C) Bonus Mass-Saving Check: What if we swap the bulky 1000 uF THT cap (~2.0g) for a tiny 330 uF cap (~0.6g) with 34 AWG?

print("\n=== 3. Full Driver-Board Capacitor + Solid-Iron Eddy-Current + Drop Kinematics ===")

sim_cases = [
    ("32AWG_Tuned_11p5V_1000uF", 32, 11.5, 1000e-6, 20.0),
    ("32AWG_Untuned_20V_1000uF", 32, 20.0, 1000e-6, 25.0),
    ("34AWG_Tuned_16V_1000uF",   34, 16.0, 1000e-6, 20.0),
    ("34AWG_LightCap_19V_330uF", 34, 19.0, 330e-6,  20.0),
]

results_tr = {}
for tag, awg, v0, cap_f, p_ms in sim_cases:
    mdl_c = gauge_models[awg][0]
    ckt = CircuitParams(v0_cap_v=v0, c1_f=cap_f, pulse_ms=p_ms)
    tr = solve_transient_pulse(mdl_c, ckt, t_end_ms=45.0, dt_ms=0.2, zippy_mass_g=27.0)
    results_tr[tag] = tr

    png_name = f"fig_03_transient_{tag}.png"
    plot_transient_summary(
        tr,
        str(out_dir / png_name),
        title=f"Exp 03 Drop Sim ({tag}): V0={v0}V, C1={int(cap_f*1e6)}uF, Pulse={p_ms:.0f}ms",
    )

    # Measure time when Zippy first clears 0.5 mm and 2.0 mm
    idx_05 = np.where(tr.zippy_disp_mm >= 0.5)[0]
    t_05 = f"{tr.t_s[idx_05[0]]*1e3:.1f} ms" if len(idx_05) > 0 else "never"

    print(
        f"  {tag:26s} -> Peak I = {tr.i_coil_a.max():4.2f} A | "
        f"Init Hold = {tr.force_z_n[0]/9.80665*1e3:5.1f} gf | "
        f"Max Overdrive = {tr.force_z_n[2:80].max()/9.80665*1e3:5.1f} gf | "
        f"0.5mm clear @ {t_05:>7s} | Sep @ 45ms = {tr.zippy_disp_mm[-1]:5.2f} mm | Dropped = {tr.dropped}"
    )

# Side-by-side overlay plot of all 4 circuit/gauge transient runs
fig, axes = plt.subplots(4, 1, figsize=(8.5, 9.2), dpi=150, sharex=True)
colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a"]
for (tag, awg, v0, cap_f, p_ms), col in zip(sim_cases, colors):
    tr = results_tr[tag]
    t_ms = tr.t_s * 1e3
    lbl = f"{awg} AWG, V0={v0}V, C1={int(cap_f*1e6)}uF"
    axes[0].plot(t_ms, tr.i_coil_a, lw=1.8, color=col, label=lbl)
    axes[1].plot(t_ms, tr.v_cap_v, lw=1.8, color=col, label=lbl)
    axes[2].plot(t_ms, tr.force_z_n / 9.80665 * 1e3, lw=1.8, color=col, label=lbl)
    axes[3].plot(t_ms, tr.zippy_disp_mm, lw=1.8, color=col, label=lbl)

axes[0].set_ylabel("Coil Current (A)")
axes[0].set_title("Exp 03 Wire Gauge & Capacitor Tuning Comparison (L_core = 7.0 mm Horn)")
axes[0].grid(True, alpha=0.3)
axes[0].legend(fontsize=8.2, loc="upper right")

axes[1].set_ylabel("Cap Voltage (V)")
axes[1].grid(True, alpha=0.3)

axes[2].axhline(27.0, color="red", ls="--", lw=1.2, label="27g Zippy Release Threshold")
axes[2].set_ylabel("Mag Pull (gf)")
axes[2].set_ylim(-5.0, 330.0)
axes[2].grid(True, alpha=0.3)
axes[2].legend(fontsize=8.2, loc="upper right")

axes[3].set_ylabel("Drop Dist (mm)")
axes[3].set_xlabel("Time (ms)")
axes[3].grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(out_dir / "fig_03_gauge_and_cap_comparison.png")
plt.close(fig)

print(f"\nAll Exp 03 field maps and transient plots saved in: {out_dir}")
