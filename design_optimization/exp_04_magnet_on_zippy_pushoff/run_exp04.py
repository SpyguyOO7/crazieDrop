import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse
import csv
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import ngsolve as ngs

from materials import MAG_IRON_4000, Material, make_ndfeb, make_wire_spec
from geometry import AxiPart, COIL_MAT, rect_rz, build_axi_model
from solvers import CircuitParams, StaticSolverSession, solve_transient_pulse

ngs.SetNumThreads(1)

import config as cfg


# --- Helper: Weak Ferrite Magnet Material (Cannot Saturate Iron Core!) ---

def make_ferrite_c8(sign: float = 1.0) -> Material:
    return Material(
        name="ferrite_C8",
        kind="pm",
        density_kg_m3=4900.0,
        mu_r=1.08,
        sigma_s_m=0.0,
        br_z_t=sign * 0.39,
        h_knee_a_m=240e3,
    )


# --- Build the 6 Candidate Architectures ---
# All use the 7.0mm Horn Dropper Core + Existing 32 AWG Coil (~197 turns, 2.95 Ohm) so comparison is apples-to-apples!
#
#   1. "inv_n52_touch_15um":      NO top magnet! Bare N52 on Zippy @ 15 um contact (tests raw mirror-image trap!)
#   2. "inv_n52_spacer_0p55mm":   NO top magnet! N52 on Zippy + 0.55 mm 3D-printed tip spacer (~300 gf hold!)
#   3. "inv_n52_iron_behind":     NO top magnet! N52 on Zippy + 1.0mm iron backing washer behind PM @ 0.65mm spacer
#   4. "inv_iron_cap_on_n52":     NO top magnet! 1.5mm iron cap IN FRONT of Zippy's N52 @ 15um (your saturation-on-Zippy test!)
#   5. "inv_weaker_n35_0p40mm":   NO top magnet! Weaker N35 magnet on Zippy @ 0.40 mm spacer (~300 gf hold)
#   6. "inv_ferrite_touch_15um":  NO top magnet! Weak Ceramic C8 magnet (Br=0.39T, can't saturate core) @ 15um contact
#   7. "dual_n52_spacer_0p85mm":  Top N52 on Dropper + N52 on Zippy @ 0.85 mm spacer (~300 gf hold)

CANDIDATE_CONFIGS = [
    "inv_n52_touch_15um",
    "inv_n52_spacer_0p55mm",
    "inv_n52_iron_behind",
    "inv_iron_cap_on_n52",
    "inv_weaker_n35_0p40mm",
    "inv_ferrite_touch_15um",
    "dual_n52_spacer_0p85mm",
]


def build_architecture(arch_name: str):
    r_s = cfg.R_STOCK_MM
    c_len = cfg.CORE_LEN_MM
    t_len = cfg.TAPER_LEN_MM
    r_tip = cfg.R_TIP_MM

    # Realistic integer-layer 32 AWG packing in 6.30 x 1.80 mm window (197 turns, R = 2.95 Ohm)
    win_len = c_len - 2.0 * cfg.BOBBIN_WALL_MM
    ws_raw = make_wire_spec(cfg.WIRE_AWG, fill_factor=0.58, ins_ratio=1.09)
    n_lay = int(np.floor(cfg.COIL_W_MM / ws_raw.d_ins_mm))
    tpl = int(np.floor(win_len / ws_raw.d_ins_mm))
    n_real = int(np.round(0.88 * n_lay * tpl))
    eff_fill = (n_real * ws_raw.area_cu_mm2) / (cfg.COIL_W_MM * win_len)
    w32 = make_wire_spec(cfg.WIRE_AWG, fill_factor=eff_fill, ins_ratio=1.09)

    core_pts = (
        (0.0, 0.0),
        (r_tip, 0.0),
        (r_s, t_len),
        (r_s, c_len),
        (0.0, c_len),
    )

    parts = [
        AxiPart("core_rod", core_pts, MAG_IRON_4000, maxh_mm=0.20, group="dropper"),
        rect_rz(
            "choke_coil",
            r_s + cfg.BOBBIN_WALL_MM,
            r_s + cfg.BOBBIN_WALL_MM + cfg.COIL_W_MM,
            cfg.BOBBIN_WALL_MM,
            c_len - cfg.BOBBIN_WALL_MM,
            COIL_MAT,
            maxh_mm=0.28,
            group="dropper",
            wire=w32,
            coil_sign=-1.0,  # Pushes -z flux down the core to oppose +z target magnet!
        ),
    ]

    if arch_name == "inv_n52_touch_15um":
        gap = 0.018  # 18 um Ni-plating / lathe contact
        parts.append(rect_rz("zippy_pm", 0.0, r_s, -cfg.PM_HEIGHT_MM - gap, -gap, make_ndfeb("N52", +1.0), maxh_mm=0.18, group="target"))

    elif arch_name == "inv_n52_spacer_0p55mm":
        gap = 0.55
        parts.append(rect_rz("zippy_pm", 0.0, r_s, -cfg.PM_HEIGHT_MM - gap, -gap, make_ndfeb("N52", +1.0), maxh_mm=0.18, group="target"))

    elif arch_name == "inv_n52_iron_behind":
        gap = 0.65
        parts.append(rect_rz("zippy_pm", 0.0, r_s, -cfg.PM_HEIGHT_MM - gap, -gap, make_ndfeb("N52", +1.0), maxh_mm=0.18, group="target"))
        parts.append(rect_rz("zippy_back_iron", 0.0, r_s, -cfg.PM_HEIGHT_MM - 1.0 - gap, -cfg.PM_HEIGHT_MM - gap, MAG_IRON_4000, maxh_mm=0.20, group="target"))

    elif arch_name == "inv_iron_cap_on_n52":
        gap = 0.015  # Iron cap on top of N52 magnet touching dropper horn at 15 um!
        cap_h = 1.50
        parts.append(rect_rz("zippy_iron_cap", 0.0, r_s, -cap_h - gap, -gap, MAG_IRON_4000, maxh_mm=0.18, group="target"))
        parts.append(rect_rz("zippy_pm", 0.0, r_s, -cap_h - cfg.PM_HEIGHT_MM - gap, -cap_h - gap, make_ndfeb("N52", +1.0), maxh_mm=0.18, group="target"))

    elif arch_name == "inv_weaker_n35_0p40mm":
        gap = 0.40
        parts.append(rect_rz("zippy_pm", 0.0, r_s, -cfg.PM_HEIGHT_MM - gap, -gap, make_ndfeb("N35", +1.0), maxh_mm=0.18, group="target"))

    elif arch_name == "inv_ferrite_touch_15um":
        gap = 0.015
        parts.append(rect_rz("zippy_pm", 0.0, r_s, -3.0 - gap, -gap, make_ferrite_c8(+1.0), maxh_mm=0.20, group="target"))

    elif arch_name == "dual_n52_spacer_0p85mm":
        gap = 0.85
        parts.append(rect_rz("top_pm", 0.0, r_s, c_len, c_len + cfg.PM_HEIGHT_MM, make_ndfeb("N52", +1.0), maxh_mm=0.18, group="dropper"))
        parts.append(rect_rz("zippy_pm", 0.0, r_s, -cfg.PM_HEIGHT_MM - gap, -gap, make_ndfeb("N52", +1.0), maxh_mm=0.18, group="target"))

    else:
        raise ValueError(f"Unknown architecture: {arch_name}")

    return build_axi_model(parts, air_pad_r_mm=28.0, air_pad_z_mm=28.0, air_maxh_mm=2.5)


# --- Task 1 Worker: Full Static Current Curve F_z(I) from 0 A to 8 A ---

def eval_static_curve(arch_name: str) -> list[dict]:
    mdl = build_architecture(arch_name)
    sess = StaticSolverSession(mdl, fe_order=2)
    m_drop = mdl.mass_breakdown_g("dropper")["total"]
    m_zip = mdl.mass_breakdown_g("target")["total"]
    coil_p = [p for p in mdl.parts if p.mat.kind == "coil"][0]

    rows = []
    for amp in cfg.CURRENTS_A:
        res = sess.solve(float(amp), compute_diagnostics=True)
        rows.append({
            "arch": arch_name,
            "current_a": float(amp),
            "ni_at": round(float(amp) * coil_p.n_turns, 1),
            "force_z_n": round(res.force_z_n, 5),
            "force_z_gf": round(res.force_z_n / 9.80665 * 1000.0, 2),
            "b_gap_t": round(res.b_gap_center_t, 4),
            "b_core_max_t": round(res.b_core_max_t, 3),
            "pm_demag_margin_ka_m": round(res.pm_demag_margin_a_m / 1e3, 1),
            "m_dropper_g": round(m_drop, 3),
            "m_zippy_g": round(m_zip, 3),
        })
    return rows


# --- Task 2 Worker: Transient 32 AWG + 1000 uF Trimpot Voltage Margin Sweep ---

def eval_transient_job(job: tuple[str, float]) -> dict:
    arch_name, v0 = job
    t0 = time.perf_counter()
    mdl = build_architecture(arch_name)
    ckt = CircuitParams(v0_cap_v=v0, c1_f=cfg.CAP_C1_F, pulse_ms=cfg.PULSE_MS)
    tr = solve_transient_pulse(mdl, ckt, t_end_ms=40.0, dt_ms=0.2, zippy_mass_g=27.0)

    min_f_gf = float(tr.force_z_n.min() / 9.80665 * 1000.0)
    hold_gf = float(tr.force_z_n[0] / 9.80665 * 1000.0)
    idx_1mm = np.where(tr.zippy_disp_mm >= 1.0)[0]
    t_1mm_ms = float(tr.t_s[idx_1mm[0]] * 1e3) if len(idx_1mm) > 0 else -1.0

    return {
        "arch": arch_name,
        "v0_v": v0,
        "hold_gf": round(hold_gf, 1),
        "peak_i_a": round(float(tr.i_coil_a.max()), 2),
        "min_force_gf": round(min_f_gf, 1),  # Negative = ACTIVE MAGNETIC REPULSION PUSH-OFF!
        "t_to_1mm_ms": round(t_1mm_ms, 2),
        "sep_40ms_mm": round(float(tr.zippy_disp_mm[-1]), 2),
        "exit_vel_m_s": round(float(tr.zippy_vel_m_s[-1]), 3),
        "dropped": int(tr.dropped),
        "wall_s": round(time.perf_counter() - t0, 2),
    }


def main():
    parser = argparse.ArgumentParser(description="Exp 04: Magnet-on-Zippy & Active Push-Off Sweep")
    parser.add_argument("--workers", type=int, default=cfg.DEFAULT_WORKERS)
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    t_start = time.perf_counter()

    # Phase 1: Static F_z(I) curves across all 7 architectures in parallel
    print(f"=== Phase 1: Solving Static Push-Off Curves (7 architectures x {len(cfg.CURRENTS_A)} currents) ===")
    all_static_rows = []
    with ProcessPoolExecutor(max_workers=min(args.workers, len(CANDIDATE_CONFIGS))) as pool:
        futs = {pool.submit(eval_static_curve, name): name for name in CANDIDATE_CONFIGS}
        for fut in as_completed(futs):
            rows = fut.result()
            all_static_rows.extend(rows)
            arch = rows[0]["arch"]
            f0 = rows[0]["force_z_gf"]
            f_min = min(r["force_z_gf"] for r in rows)
            min_demag = min(r["pm_demag_margin_ka_m"] for r in rows)
            print(
                f"  [Static Done] {arch:24s} | m_drop={rows[0]['m_dropper_g']:.2f}g, m_zip={rows[0]['m_zippy_g']:.2f}g | "
                f"Hold(0A)={f0:6.1f} gf -> Min F_z(0..8A)={f_min:+7.1f} gf | Worst Demag={min_demag:+6.1f} kA/m"
            )

    all_static_rows.sort(key=lambda r: (CANDIDATE_CONFIGS.index(r["arch"]), r["current_a"]))
    static_csv = out_dir / "static_curves.csv"
    with static_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_static_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_static_rows)

    # Phase 2: Transient 32 AWG + 1000 uF Trimpot Margin Sweep across the most promising architectures
    transient_archs = [
        "inv_n52_spacer_0p55mm",
        "inv_weaker_n35_0p40mm",
        "inv_iron_cap_on_n52",
        "dual_n52_spacer_0p85mm",
    ]
    jobs = [(a, float(v)) for a in transient_archs for v in cfg.TRIMPOT_VOLTAGES_V]
    print(f"\n=== Phase 2: Running {len(jobs)} Moving-Gap Transient Capacitor Drop/Launch Sims on {args.workers} Workers ===")

    tr_rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(eval_transient_job, j): j for j in jobs}
        done = 0
        for fut in as_completed(futs):
            row = fut.result()
            tr_rows.append(row)
            done += 1
            status = "DROPPED/LAUNCHED!" if row["dropped"] else "STUCK"
            print(
                f"  [{done:2d}/{len(jobs):2d}] {row['arch']:24s} @ V0={row['v0_v']:4.1f}V -> "
                f"PeakI={row['peak_i_a']:4.2f}A | MinF={row['min_force_gf']:+6.1f}gf | "
                f"1mm@{row['t_to_1mm_ms']:5.1f}ms | Sep40ms={row['sep_40ms_mm']:5.2f}mm | {status}"
            )

    tr_rows.sort(key=lambda r: (transient_archs.index(r["arch"]), r["v0_v"]))
    tr_csv = out_dir / "transient_margin.csv"
    with tr_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(tr_rows[0].keys()))
        writer.writeheader()
        writer.writerows(tr_rows)

    print(f"\nAll Exp 04 sweeps completed in {time.perf_counter() - t_start:.1f} s! Run plot_exp04.py next.")


if __name__ == "__main__":
    main()
