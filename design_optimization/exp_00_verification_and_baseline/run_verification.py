from pathlib import Path
import numpy as np

from materials import MU_0, STEEL_1018, MAG_IRON_4000, LATHE_FINISH_GAP_MM, make_ndfeb, make_wire_spec
from geometry import AxiPart, COIL_MAT, rect_rz, build_axi_model
from solvers import CircuitParams, solve_static, solve_transient_pulse, _eggshell_weight_cfs
from plotting import plot_field_map, plot_transient_summary
import ngsolve as ngs


# --- Output Directory Setup ---

out_dir = Path(__file__).resolve().parent


# --- 1. Analytical Benchmark 1: Standalone Cylindrical Permanent Magnet On-Axis B_z(z) ---

print("--- Benchmark 1: Cylindrical Permanent Magnet On-Axis Field ---")
pm_r_mm = 3.175
pm_h_mm = 3.0
br_test = 1.30

pm_ideal = make_ndfeb("N42", sign=1.0)
pm_ideal = type(pm_ideal)(
    name="pm_ideal",
    kind="pm",
    density_kg_m3=7500.0,
    mu_r=1.0,
    br_z_t=br_test,
    h_knee_a_m=900e3,
)

pm_part = rect_rz("pm", 0.0, pm_r_mm, 0.0, pm_h_mm, pm_ideal, maxh_mm=0.2)
bm1_model = build_axi_model([pm_part], air_pad_r_mm=45.0, air_pad_z_mm=45.0, air_maxh_mm=3.0)
bm1_res = solve_static(bm1_model, coil_current_a=0.0, fe_order=3)

z_probe_mm = np.array([0.5, 1.5, 3.0, 5.0]) + pm_h_mm
d_mm = z_probe_mm - pm_h_mm
bz_exact = 0.5 * br_test * (
    (d_mm + pm_h_mm) / np.sqrt(pm_r_mm**2 + (d_mm + pm_h_mm) ** 2)
    - d_mm / np.sqrt(pm_r_mm**2 + d_mm**2)
)
bz_fem = np.array([float(bm1_res.b_vec_cf[1](bm1_model.mesh(2e-4, z * 1e-3))) for z in z_probe_mm])
rel_err_pm = np.max(np.abs(bz_fem - bz_exact) / bz_exact) * 100.0

print(f"Exact B_z (T): {np.round(bz_exact, 5)}")
print(f"FEM   B_z (T): {np.round(bz_fem, 5)}")
print(f"Max relative error: {rel_err_pm:.3f}%")
assert rel_err_pm < 1.5, f"PM benchmark error too high: {rel_err_pm:.3f}%"


# --- 2. Analytical Benchmark 2: Finite Thick Air-Core Solenoid On-Axis Field ---

print("\n--- Benchmark 2: Thick Solenoid Center Field ---")
sol_r1, sol_r2, sol_len = 3.2, 5.2, 9.0
n_turns_test = 200.0
i_test = 2.0

sol_part = rect_rz(
    "coil", sol_r1, sol_r2, -0.5 * sol_len, 0.5 * sol_len, COIL_MAT, maxh_mm=0.25, turns=n_turns_test, coil_sign=1.0
)
bm2_model = build_axi_model([sol_part], air_pad_r_mm=45.0, air_pad_z_mm=45.0, air_maxh_mm=3.0)
bm2_res = solve_static(bm2_model, coil_current_a=i_test, fe_order=3)

j_dens = (n_turns_test * i_test) / ((sol_r2 - sol_r1) * sol_len * 1e-6)
bz_sol_exact = (
    0.5
    * MU_0
    * j_dens
    * (sol_len * 1e-3)
    * np.log(
        (sol_r2 + np.sqrt(sol_r2**2 + (0.5 * sol_len) ** 2))
        / (sol_r1 + np.sqrt(sol_r1**2 + (0.5 * sol_len) ** 2))
    )
)
bz_sol_fem = float(bm2_res.b_vec_cf[1](bm2_model.mesh(2e-4, 0.0)))
rel_err_sol = abs(bz_sol_fem - bz_sol_exact) / bz_sol_exact * 100.0

print(f"Exact B_z(0) (T): {bz_sol_exact:.6f} | FEM B_z(0) (T): {bz_sol_fem:.6f} | Error: {rel_err_sol:.3f}%")
assert rel_err_sol < 1.0, f"Solenoid benchmark error too high: {rel_err_sol:.3f}%"


# --- 3. Non-Cylindrical Radially Symmetric Profile & Virtual-Work Cross-Check ---

print("\n--- Benchmark 3: Non-Cylindrical Tapered Core Profile & MST vs Virtual-Work Cross-Check ---")
gap0_mm = 0.25
dz_vw_m = 1e-5  # 10 um smooth mesh deformation for zero-remesh-noise Virtual Work check

core_pts = (
    (0.0, 0.0),
    (2.3, 0.0),
    (3.175, 2.5),
    (3.175, 9.5),
    (0.0, 9.5),
)

taper_parts = [
    AxiPart("tapered_core", core_pts, STEEL_1018, maxh_mm=0.22, group="dropper"),
    rect_rz("top_mag", 0.0, 3.175, 9.5, 11.0, make_ndfeb("N42", sign=1.0), maxh_mm=0.2, group="dropper"),
    rect_rz("zippy_slug", 0.0, 3.175, -3.0 - gap0_mm, -gap0_mm, STEEL_1018, maxh_mm=0.2, group="target"),
]
taper_model = build_axi_model(taper_parts, air_pad_r_mm=30.0, air_pad_z_mm=25.0, air_maxh_mm=2.5)
taper_res = solve_static(taper_model, coil_current_a=0.0, fe_order=3)

# Deform the EXACT SAME mesh smoothly by +/- dz * w(r, z) in z so remeshing noise is identically zero!
w_cf, _, _ = _eggshell_weight_cfs(taper_model)
fes_def = ngs.VectorH1(taper_model.mesh, order=1)
gf_def = ngs.GridFunction(fes_def)

gf_def.Set(ngs.CF((0.0, +dz_vw_m * w_cf)))  # Shift target UP by +dz_vw_m (smaller gap)
taper_model.mesh.SetDeformation(gf_def)
pi_up = solve_static(taper_model, coil_current_a=0.0, fe_order=3, init_psi=taper_res.psi_gf).magnetic_energy_j

gf_def.Set(ngs.CF((0.0, -dz_vw_m * w_cf)))  # Shift target DOWN by -dz_vw_m (larger gap)
taper_model.mesh.SetDeformation(gf_def)
pi_down = solve_static(taper_model, coil_current_a=0.0, fe_order=3, init_psi=taper_res.psi_gf).magnetic_energy_j
taper_model.mesh.UnsetDeformation()

# Force in +z is -dPi / dz_target
f_vw = -(pi_up - pi_down) / (2.0 * dz_vw_m)
vw_err_pct = abs(taper_res.force_z_n - f_vw) / abs(f_vw) * 100.0

print(
    f"Tapered core mass: {taper_model.parts[0].mass_g:.3f} g | "
    f"Eggshell MST F_z = {taper_res.force_z_n:.5f} N vs Exact Virtual-Work F_z = {f_vw:.5f} N (diff: {vw_err_pct:.3f}%)"
)
assert vw_err_pct < 2.0, f"MST vs Virtual Work force mismatch too high: {vw_err_pct:.3f}%"
plot_field_map(
    taper_model,
    taper_res,
    str(out_dir / "fig_tapered_core_demo.png"),
    title=f"Tapered Core Demo (Passive Hold @ {gap0_mm}mm gap): F_z = {taper_res.force_z_n:+.3f} N",
)


# --- 3.5. Benchmark 4: McMaster 1/4" x 1/10" N52 Magnet Direct Pull Rating Check (Rated 2.8 lbs = 12.45 N) ---
# Standard Ni-Cu-Ni magnet plating is ~15-20 um thick per side; check direct pull on thick iron plate:

print("\n--- Benchmark 4: McMaster 1/4\" x 1/10\" N52 Magnet Direct Pull vs 2.8 lb (12.45 N) Spec ---")
ni_coat_gap_mm = 0.018
bm4_parts = [
    rect_rz("mcmaster_pm", 0.0, 3.175, 0.0, 2.54, make_ndfeb("N52", sign=1.0), maxh_mm=0.16, group="dropper"),
    rect_rz("thick_plate", 0.0, 12.0, -8.0 - ni_coat_gap_mm, -ni_coat_gap_mm, MAG_IRON_4000, maxh_mm=0.35, group="target"),
]
bm4_model = build_axi_model(bm4_parts, air_pad_r_mm=30.0, air_pad_z_mm=25.0, air_maxh_mm=2.5)
bm4_res = solve_static(bm4_model, coil_current_a=0.0, fe_order=3)
pull_lbs = bm4_res.force_z_n / 4.44822
print(f"Simulated direct pull (@ {ni_coat_gap_mm*1e3:.0f} um Ni-plating gap): {bm4_res.force_z_n:.2f} N ({pull_lbs:.2f} lbs) vs McMaster 2.80 lbs spec")


# --- 4. Current Real Baseline crazieDrop Geometry (McMaster N52 + mu_r=4000 Mag Iron + Lathe Finish, 20V Cap) ---

print("\n--- Baseline crazieDrop Model Evaluation (McMaster 1/4\"x1/10\" N52 + mu_r=4000 Iron + 15um Lathe Finish) ---")
r_rod = 6.35 / 2.0
pm_height_mm = 2.54                 # 0.10 in McMaster N52 magnet
air_gap_mm = LATHE_FINISH_GAP_MM    # 15 um careful lathe-cut contact surface finish
wire_32 = make_wire_spec(awg=32, fill_factor=0.58)

base_parts = [
    rect_rz("core_rod", 0.0, r_rod, 0.0, 9.5, MAG_IRON_4000, maxh_mm=0.25, group="dropper"),
    rect_rz("top_pm", 0.0, r_rod, 9.5, 9.5 + pm_height_mm, make_ndfeb("N52", sign=1.0), maxh_mm=0.20, group="dropper"),
    rect_rz(
        "choke_coil",
        r_rod + 0.35,
        r_rod + 0.35 + 1.8,
        0.3,
        9.2,
        COIL_MAT,
        maxh_mm=0.35,
        group="dropper",
        wire=wire_32,
        coil_sign=-1.0,
    ),
    rect_rz("zippy_iron", 0.0, r_rod, -3.0 - air_gap_mm, -air_gap_mm, MAG_IRON_4000, maxh_mm=0.22, group="target"),
]

base_model = build_axi_model(base_parts, air_pad_r_mm=30.0, air_pad_z_mm=30.0, air_maxh_mm=2.5)
masses = base_model.mass_breakdown_g("dropper")
coil_part = base_parts[2]

print(f"Dropper active mass breakdown (g): { {k: round(v, 3) for k, v in masses.items()} }")
print(f"32 AWG Coil: N = {coil_part.n_turns:.1f} turns, R_dc = {coil_part.coil_res_ohm:.2f} Ohm")

currents = [0.0, 0.5, 1.0, 1.25, 1.5, 2.0, 3.0, 4.5]
prev_psi = None
best_null_res = None

for amp in currents:
    res = solve_static(base_model, coil_current_a=amp, fe_order=2, init_psi=prev_psi)
    prev_psi = res.psi_gf
    ni = amp * coil_part.n_turns
    print(
        f"  I = {amp:4.2f} A (NI = {ni:6.0f} A-t) -> "
        f"F_z = {res.force_z_n:+7.3f} N | B_gap = {res.b_gap_center_t:+6.3f} T | "
        f"B_core_max = {res.b_core_max_t:5.2f} T | PM_demag_margin = {res.pm_demag_margin_a_m/1e3:6.1f} kA/m"
    )
    if amp == 0.0:
        plot_field_map(
            base_model,
            res,
            str(out_dir / "fig_baseline_hold.png"),
            title=f"Baseline Passive Hold (I = 0 A): F_z = {res.force_z_n:+.3f} N, B_gap = {res.b_gap_center_t:+.3f} T",
        )
    if best_null_res is None or abs(res.force_z_n) < abs(best_null_res.force_z_n):
        best_null_res = res
    if amp == 4.5:
        plot_field_map(
            base_model,
            res,
            str(out_dir / "fig_baseline_overdrive.png"),
            title=f"Baseline Overdriven (I = {amp:.1f} A): F_z = {res.force_z_n:+.3f} N (Re-Attraction!), B_gap = {res.b_gap_center_t:+.3f} T",
        )

assert best_null_res is not None
plot_field_map(
    base_model,
    best_null_res,
    str(out_dir / "fig_baseline_choke_null.png"),
    title=(
        f"Baseline Choke Null (I = {best_null_res.coil_current_a:.2f} A): "
        f"F_z = {best_null_res.force_z_n:+.3f} N, B_gap = {best_null_res.b_gap_center_t:+.3f} T"
    ),
)

# --- 5. Tier 2 Moving-Gap Transient Simulation (V0 = 20 V on 1000 uF Cap + Eddy Currents) ---

print("\n--- Tier 2 Moving-Gap Transient Pulse Evaluation (V0 = 20.0 V) ---")
ckt_20v = CircuitParams(v0_cap_v=20.0, c1_f=1000e-6, pulse_ms=25.0)
tr_20v = solve_transient_pulse(base_model, ckt_20v, t_end_ms=45.0, dt_ms=0.2, zippy_mass_g=27.0)
plot_transient_summary(
    tr_20v,
    str(out_dir / "fig_baseline_transient.png"),
    title="Baseline Moving-Gap Transient (V0 = 20.0 V, C1 = 1000 uF, Pulse = 25 ms)",
)

print(
    f"V0 = 20.0 V | Peak I = {tr_20v.i_coil_a.max():.2f} A | "
    f"Initial Hold F_z = {tr_20v.force_z_n[0]:.3f} N | Min F_z = {tr_20v.force_z_n.min():.3f} N | "
    f"End F_z = {tr_20v.force_z_n[-1]:.4f} N | End V_cap = {tr_20v.v_cap_v[-1]:.2f} V | "
    f"Zippy separation @ 45ms = {tr_20v.zippy_disp_mm[-1]:.2f} mm | Dropped = {tr_20v.dropped}"
)
print("\nALL VERIFICATION CHECKS PASSED.")
