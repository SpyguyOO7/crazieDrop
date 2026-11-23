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

from materials import MAG_IRON_4000, make_ndfeb, make_wire_spec
from geometry import AxiPart, COIL_MAT, rect_rz, build_axi_model
from solvers import CircuitParams, StaticSolverSession, solve_transient_pulse

ngs.SetNumThreads(1)

import config as cfg


# --- Candidate Builder: Flush Top + Recessed Bottom + Scramble Hand-Wind Fill ---

def build_recessed_design(
    core_len_mm: float,
    coil_w_mm: float,
    z_coil_bot_mm: float,
    fill_factor: float,
):
    gap = cfg.CONTACT_GAP_MM
    r_pm = cfg.R_STOCK_MM

    core_pts = (
        (0.0, 0.0),
        (cfg.R_TIP_MM, 0.0),
        (cfg.R_SHANK_MM, cfg.TAPER_LEN_MM),
        (cfg.R_SHANK_MM, core_len_mm),
        (0.0, core_len_mm),
    )

    z_coil_top_mm = core_len_mm - cfg.TOP_GLUE_GAP_MM
    coil_r_in = cfg.R_SHANK_MM + cfg.BOBBIN_RADIAL_WALL_MM
    coil_r_out = coil_r_in + coil_w_mm

    w_spec = make_wire_spec(awg=cfg.WIRE_AWG, fill_factor=fill_factor, ins_ratio=1.09)

    parts = [
        AxiPart("core_rod", core_pts, MAG_IRON_4000, maxh_mm=0.20, group="dropper"),
        rect_rz(
            "choke_coil",
            coil_r_in,
            coil_r_out,
            z_coil_bot_mm,
            z_coil_top_mm,
            COIL_MAT,
            maxh_mm=0.26,
            group="dropper",
            wire=w_spec,
            coil_sign=-1.0,
        ),
        rect_rz("zippy_iron_cap", 0.0, r_pm, -cfg.ZIPPY_CAP_H_MM - gap, -gap, MAG_IRON_4000, maxh_mm=0.18, group="target"),
        rect_rz(
            "zippy_pm",
            0.0,
            r_pm,
            -cfg.ZIPPY_CAP_H_MM - cfg.PM_HEIGHT_MM - gap,
            -cfg.ZIPPY_CAP_H_MM - gap,
            make_ndfeb("N52", +1.0),
            maxh_mm=0.18,
            group="target",
        ),
    ]

    return build_axi_model(parts, air_pad_r_mm=26.0, air_pad_z_mm=26.0, air_maxh_mm=2.5)


# --- Phase 1 Static Worker (60 Combinations of Recession x Scramble Fill x Window) ---

def eval_static_task(task: tuple[int, float, float, float, float]) -> dict:
    idx, c_len, coil_w, z_bot, fill = task
    t0 = time.perf_counter()

    mdl = build_recessed_design(c_len, coil_w, z_bot, fill)
    core_p, coil_p, cap_p, pm_p = mdl.parts

    m_drop = core_p.mass_g + coil_p.mass_g
    m_zip = cap_p.mass_g + pm_p.mass_g
    r_coil = coil_p.coil_res_ohm
    r_loop = r_coil + 0.152
    tau_rc_ms = r_loop * cfg.CAP_C1_F * 1e3

    i_18v = 0.84 * (18.0 / r_loop)
    i_22v = 0.84 * (22.0 / r_loop)

    sess = StaticSolverSession(mdl, fe_order=2)
    res_0 = sess.solve(0.0, compute_diagnostics=False)
    res_18 = sess.solve(i_18v, compute_diagnostics=False)
    res_22 = sess.solve(i_22v, compute_diagnostics=True)

    f_hold_gf = res_0.force_z_n / 9.80665 * 1000.0
    f_18v_gf = res_18.force_z_n / 9.80665 * 1000.0
    f_22v_gf = res_22.force_z_n / 9.80665 * 1000.0

    return {
        "id": idx,
        "core_len_mm": c_len,
        "coil_w_mm": coil_w,
        "z_coil_bot_mm": z_bot,
        "fill_factor": fill,
        "coil_od_mm": round(2.0 * (cfg.R_SHANK_MM + cfg.BOBBIN_RADIAL_WALL_MM + coil_w), 2),
        "n_turns": int(round(coil_p.n_turns)),
        "r_coil_ohm": round(r_coil, 3),
        "tau_rc_ms": round(tau_rc_ms, 2),
        "m_core_g": round(core_p.mass_g, 3),
        "m_cu_g": round(coil_p.mass_g, 3),
        "m_dropper_g": round(m_drop, 3),
        "m_zippy_g": round(m_zip, 3),
        "f_hold_gf": round(f_hold_gf, 1),
        "i_18v_a": round(i_18v, 2),
        "f_18v_gf": round(f_18v_gf, 1),
        "f_22v_gf": round(f_22v_gf, 1),
        "impulse_18v_gf_ms": round(max(0.0, -f_18v_gf) * tau_rc_ms, 1),
        "demag_22v_ka_m": round(res_22.pm_demag_margin_a_m / 1e3, 1),
        "eval_s": round(time.perf_counter() - t0, 2),
    }


# --- Phase 2 Transient Worker: Direct Stress-Test of Bottom Recession & Messy 0.35 vs 0.42 Fill ---

def eval_transient_task(job: tuple[str, float, float, float, float, float]) -> dict:
    label, c_len, coil_w, z_bot, fill, v0 = job
    t0 = time.perf_counter()

    mdl = build_recessed_design(c_len, coil_w, z_bot, fill)
    coil_p = mdl.parts[1]
    m_drop = mdl.parts[0].mass_g + coil_p.mass_g

    ckt = CircuitParams(v0_cap_v=v0, c1_f=cfg.CAP_C1_F, pulse_ms=cfg.PULSE_MS)
    tr = solve_transient_pulse(mdl, ckt, t_end_ms=42.0, dt_ms=0.2, zippy_mass_g=27.0)

    idx_1mm = np.where(tr.zippy_disp_mm >= 1.0)[0]
    t_1mm = float(tr.t_s[idx_1mm[0]] * 1e3) if len(idx_1mm) > 0 else -1.0

    return {
        "series": label,
        "core_len_mm": c_len,
        "coil_w_mm": coil_w,
        "z_coil_bot_mm": z_bot,
        "fill_factor": fill,
        "n_turns": int(round(coil_p.n_turns)),
        "r_coil_ohm": round(coil_p.coil_res_ohm, 2),
        "m_dropper_g": round(m_drop, 3),
        "v0_v": v0,
        "hold_gf": round(float(tr.force_z_n[0] / 9.80665 * 1000.0), 1),
        "peak_i_a": round(float(tr.i_coil_a.max()), 2),
        "min_f_gf": round(float(tr.force_z_n.min() / 9.80665 * 1000.0), 1),
        "t_1mm_ms": round(t_1mm, 2),
        "sep_42ms_mm": round(float(tr.zippy_disp_mm[-1]), 2),
        "dropped": int(tr.dropped),
        "wall_s": round(time.perf_counter() - t0, 1),
    }


def main():
    parser = argparse.ArgumentParser(description="Exp 07: Scramble Hand-Winding Fill & Bottom Recession Sweep")
    parser.add_argument("--workers", type=int, default=cfg.DEFAULT_WORKERS)
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    t_start = time.perf_counter()

    # Phase 1: Full 60-point Static Grid
    tasks = []
    idx = 0
    for c_len in cfg.CORE_LEN_LIST_MM:
        for c_w in cfg.COIL_W_LIST_MM:
            for z_bot in cfg.Z_COIL_BOT_LIST_MM:
                for fill in cfg.SCRAMBLE_FILL_LIST:
                    tasks.append((idx, c_len, c_w, z_bot, fill))
                    idx += 1

    print(f"=== Phase 1: Running {len(tasks)} Static Scramble-Wind & Recession Evaluations ({args.workers} workers) ===")
    s_rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(eval_static_task, t): t for t in tasks}
        for fut in as_completed(futs):
            row = fut.result()
            s_rows.append(row)

    s_rows.sort(key=lambda r: r["id"])
    with (out_dir / "static_recession_grid.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(s_rows[0].keys()))
        writer.writeheader()
        writer.writerows(s_rows)

    print("Phase 1 Summary (@ fill=0.42 typical scramble wind, V0=18V):")
    for r in s_rows:
        if abs(r["fill_factor"] - 0.42) < 1e-4 and r["z_coil_bot_mm"] in (0.5, 1.3, 2.2):
            print(
                f"  L={r['core_len_mm']}mm, w={r['coil_w_mm']}mm, z_bot={r['z_coil_bot_mm']}mm | "
                f"N={r['n_turns']:3d}t ({r['r_coil_ohm']:4.2f}Ω, m={r['m_dropper_g']:.2f}g) -> "
                f"Hold={r['f_hold_gf']:.0f}gf | F_18V={r['f_18v_gf']:+6.1f}gf | Impulse={r['impulse_18v_gf_ms']:5.0f}gf*ms"
            )

    # Phase 2: Moving-Gap Transient Launch Matrix (30 jobs — runs in ~2 batches on 16 workers!)
    # Directly answers:
    #   Q1: How far can we recess the bottom coil (z_bot = 0.5, 0.9, 1.3, 1.7, 2.2 mm) at typical scramble fill (0.42)?
    #   Q2: What if the hand-wind is SUPER messy (fill = 0.35) with a generous 1.3 mm bottom recession?
    #   Q3: Does widening the bobbin slightly (w = 2.3 mm, OD = 9.7 mm) make even a super-messy 0.35 fill + 1.3 mm recession bulletproof?
    tr_configs = [
        ("Recess_0p5mm_Fill0p42_w1p8", 8.5, 1.8, 0.50, 0.42),
        ("Recess_0p9mm_Fill0p42_w1p8", 8.5, 1.8, 0.90, 0.42),
        ("Recess_1p3mm_Fill0p42_w1p8", 8.5, 1.8, 1.30, 0.42),
        ("Recess_1p7mm_Fill0p42_w1p8", 8.5, 1.8, 1.70, 0.42),
        ("MessyFill_0p35_Recess1p3_w1p8", 8.5, 1.8, 1.30, 0.35),
        ("MessyFill_0p35_Recess1p3_w2p3", 8.5, 2.3, 1.30, 0.35),
    ]
    tr_jobs = [
        (lbl, cl, cw, zb, fl, float(v))
        for (lbl, cl, cw, zb, fl) in tr_configs
        for v in cfg.TRANSIENT_VOLTAGES_V
    ]

    print(f"\n=== Phase 2: Running {len(tr_jobs)} Moving-Gap Transient Launch Sims Across {args.workers} Workers ===")
    tr_rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(eval_transient_task, j): j for j in tr_jobs}
        done = 0
        for fut in as_completed(futs):
            row = fut.result()
            tr_rows.append(row)
            done += 1
            status = "LAUNCHED!" if row["dropped"] else "stuck"
            print(
                f"  [{done:2d}/{len(tr_jobs):2d}] {row['series']:30s} @ {row['v0_v']:4.1f}V | "
                f"N={row['n_turns']}t ({row['m_dropper_g']:.2f}g) -> "
                f"MinF={row['min_f_gf']:+6.1f}gf | Sep42ms={row['sep_42ms_mm']:5.2f}mm | {status}"
            )

    order_map = {c[0]: i for i, c in enumerate(tr_configs)}
    tr_rows.sort(key=lambda r: (order_map[r["series"]], r["v0_v"]))
    with (out_dir / "transient_recession_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(tr_rows[0].keys()))
        writer.writeheader()
        writer.writerows(tr_rows)

    print(f"\nExp 07 complete in {time.perf_counter() - t_start:.1f} s! Run plot_exp07.py next.")


if __name__ == "__main__":
    main()
