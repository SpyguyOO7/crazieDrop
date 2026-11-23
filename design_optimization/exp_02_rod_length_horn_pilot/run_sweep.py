import os

# Keep each worker single-threaded so N parallel processes scale linearly across Zen 5 cores
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse
import csv
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np

from materials import MAG_IRON_4000, make_ndfeb, make_wire_spec
from geometry import AxiPart, COIL_MAT, rect_rz, build_axi_model
from solvers import StaticSolverSession
import ngsolve as ngs

ngs.SetNumThreads(1)

import config as cfg


# --- Geometry Builder for Any (interface_mode, core_len_mm, coil_width_mm) ---

def build_candidate(interface_mode: str, core_len_mm: float, coil_w_mm: float):
    r_stock = cfg.R_STOCK_MM
    r_shank = cfg.R_SHANK_MM
    r_tip = cfg.R_TIP_MM
    gap = cfg.CONTACT_GAP_MM

    if interface_mode == "flat_cyl":
        core_pts = (
            (0.0, 0.0),
            (r_stock, 0.0),
            (r_stock, core_len_mm),
            (0.0, core_len_mm),
        )
        coil_z_bot = cfg.BOBBIN_WALL_MM
    else:
        # Bottom horn taper (0 .. HORN_TAPER_LEN_MM) + straight shank + top collector pad
        z_taper = cfg.HORN_TAPER_LEN_MM
        z_coll_bot = core_len_mm - cfg.TOP_COLLECTOR_LEN_MM
        if abs(r_shank - r_stock) < 1e-6:
            core_pts = (
                (0.0, 0.0),
                (r_tip, 0.0),
                (r_shank, z_taper),
                (r_stock, core_len_mm),
                (0.0, core_len_mm),
            )
        else:
            core_pts = (
                (0.0, 0.0),
                (r_tip, 0.0),
                (r_shank, z_taper),
                (r_shank, z_coll_bot),
                (r_stock, core_len_mm),
                (0.0, core_len_mm),
            )
        # Coil sits on straight shank above the bottom horn taper (or starts just above z=0.35mm if outer OD clears)
        coil_z_bot = cfg.BOBBIN_WALL_MM

    coil_z_top = core_len_mm - cfg.BOBBIN_WALL_MM
    coil_r_in = max(r_shank, r_stock) + cfg.BOBBIN_WALL_MM
    coil_r_out = coil_r_in + coil_w_mm

    if interface_mode == "double_horn":
        tgt_pts = (
            (0.0, -cfg.ZIPPY_LEN_MM - gap),
            (r_stock, -cfg.ZIPPY_LEN_MM - gap),
            (r_stock, -cfg.ZIPPY_TAPER_LEN_MM - gap),
            (r_tip, -gap),
            (0.0, -gap),
        )
    else:
        tgt_pts = (
            (0.0, -cfg.ZIPPY_LEN_MM - gap),
            (r_stock, -cfg.ZIPPY_LEN_MM - gap),
            (r_stock, -gap),
            (0.0, -gap),
        )

    wire = make_wire_spec(awg=cfg.WIRE_AWG, fill_factor=cfg.WIRE_FILL)

    parts = [
        AxiPart("core_rod", core_pts, MAG_IRON_4000, maxh_mm=0.22, group="dropper"),
        rect_rz(
            "top_pm",
            0.0,
            r_stock,
            core_len_mm,
            core_len_mm + cfg.PM_HEIGHT_MM,
            make_ndfeb(cfg.PM_GRADE, sign=1.0),
            maxh_mm=0.18,
            group="dropper",
        ),
        rect_rz(
            "choke_coil",
            coil_r_in,
            coil_r_out,
            coil_z_bot,
            coil_z_top,
            COIL_MAT,
            maxh_mm=0.30,
            group="dropper",
            wire=wire,
            coil_sign=-1.0,
        ),
        AxiPart("zippy_iron", tgt_pts, MAG_IRON_4000, maxh_mm=0.20, group="target"),
    ]

    return build_axi_model(parts, air_pad_r_mm=28.0, air_pad_z_mm=28.0, air_maxh_mm=2.5)


# --- Single Worker Evaluation (Passive Hold + Exact Choke-Null Secant Root-Find) ---

def eval_single_task(task: tuple[int, str, float, float]) -> dict:
    idx, mode, core_len_mm, coil_w_mm = task
    t0 = time.perf_counter()

    model = build_candidate(mode, core_len_mm, coil_w_mm)
    core_p, pm_p, coil_p, zippy_p = model.parts
    sess = StaticSolverSession(model, fe_order=2)

    # 1. Passive Hold State (I = 0 A)
    res0 = sess.solve(coil_current_a=0.0, compute_diagnostics=True)
    f_hold_n = res0.force_z_n
    f_hold_gf = f_hold_n / 9.80665 * 1000.0
    b_gap_hold_t = res0.b_gap_center_t

    # 2. Secant Root-Find on B_gap_center(I) == 0 reusing the exact same FE space & mesh
    n_turns = coil_p.n_turns
    r_coil = coil_p.coil_res_ohm

    i_a = 0.0
    b_a = b_gap_hold_t

    i_b = 380.0 / max(n_turns, 1.0)
    res_b = sess.solve(coil_current_a=i_b, compute_diagnostics=False)
    b_b = res_b.b_gap_center_t

    best_i = i_b
    for _ in range(5):
        if abs(b_b) < 2e-3 or abs(b_b - b_a) < 1e-6:
            break
        i_next = i_b - b_b * (i_b - i_a) / (b_b - b_a)
        i_next = float(np.clip(i_next, 0.02, 4000.0 / max(n_turns, 1.0)))
        res_next = sess.solve(coil_current_a=i_next, compute_diagnostics=False)
        i_a, b_a = i_b, b_b
        i_b, b_b = i_next, res_next.b_gap_center_t
        best_i = i_next

    # Full diagnostics computed on the final null operating point
    best_res = sess.solve(coil_current_a=best_i, compute_diagnostics=True)
    for mult in (0.96, 1.04):
        cand_res = sess.solve(coil_current_a=best_i * mult, compute_diagnostics=False)
        if abs(cand_res.force_z_n) < abs(best_res.force_z_n):
            best_res = sess.solve(coil_current_a=best_i * mult, compute_diagnostics=True)

    i_null_a = best_res.coil_current_a
    ni_null_at = i_null_a * n_turns
    v_dc_null_v = i_null_a * r_coil
    p_dc_null_w = (i_null_a ** 2) * r_coil
    f_null_gf = best_res.force_z_n / 9.80665 * 1000.0

    m_dropper_g = core_p.mass_g + pm_p.mass_g + coil_p.mass_g
    elapsed_s = time.perf_counter() - t0

    return {
        "task_id": idx,
        "interface_mode": mode,
        "core_len_mm": round(core_len_mm, 3),
        "coil_w_mm": round(coil_w_mm, 3),
        "m_core_g": round(core_p.mass_g, 4),
        "m_pm_g": round(pm_p.mass_g, 4),
        "m_coil_g": round(coil_p.mass_g, 4),
        "m_dropper_active_g": round(m_dropper_g, 4),
        "m_zippy_g": round(zippy_p.mass_g, 4),
        "n_turns_32awg": round(n_turns, 1),
        "r_coil_32awg_ohm": round(r_coil, 3),
        "f_hold_n": round(f_hold_n, 5),
        "f_hold_gf": round(f_hold_gf, 2),
        "hold_per_gram_gf_g": round(f_hold_gf / m_dropper_g, 2),
        "b_gap_hold_t": round(b_gap_hold_t, 4),
        "b_core_max_hold_t": round(res0.b_core_max_t, 3),
        "i_null_a": round(i_null_a, 4),
        "ni_null_at": round(ni_null_at, 1),
        "v_dc_null_v": round(v_dc_null_v, 3),
        "p_dc_null_w": round(p_dc_null_w, 3),
        "f_null_res_gf": round(f_null_gf, 3),
        "pm_demag_margin_null_ka_m": round(best_res.pm_demag_margin_a_m / 1e3, 1),
        "eval_time_s": round(elapsed_s, 3),
    }


# --- CLI Multi-Core Runner ---

def main():
    parser = argparse.ArgumentParser(description="Multi-core rod-length & horn-tip optimization sweep.")
    parser.add_argument("--workers", type=int, default=cfg.DEFAULT_WORKERS, help="Number of parallel worker processes")
    parser.add_argument("--resume", action="store_true", help="Skip (mode, core_len, coil_w) combinations already in results.csv")
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    csv_path = out_dir / "results.csv"

    done_keys = set()
    if args.resume and csv_path.exists():
        with csv_path.open("r", newline="") as f:
            for row in csv.DictReader(f):
                done_keys.add((row["interface_mode"], float(row["core_len_mm"]), float(row["coil_w_mm"])))

    tasks = []
    idx = 0
    for mode in cfg.INTERFACE_MODES:
        for l_core in cfg.CORE_LENGTHS_MM:
            for w_coil in cfg.COIL_WIDTHS_MM:
                key = (mode, round(float(l_core), 3), round(float(w_coil), 3))
                if key not in done_keys:
                    tasks.append((idx, mode, float(l_core), float(w_coil)))
                idx += 1

    total = len(tasks)
    if total == 0:
        print("All sweep points already exist in results.csv! Run plot_results.py to generate plots.")
        return

    print(f"Launching {total} nonlinear FEM evaluations across {args.workers} parallel processes...")
    t_start = time.perf_counter()

    fieldnames = [
        "task_id", "interface_mode", "core_len_mm", "coil_w_mm",
        "m_core_g", "m_pm_g", "m_coil_g", "m_dropper_active_g", "m_zippy_g",
        "n_turns_32awg", "r_coil_32awg_ohm",
        "f_hold_n", "f_hold_gf", "hold_per_gram_gf_g", "b_gap_hold_t", "b_core_max_hold_t",
        "i_null_a", "ni_null_at", "v_dc_null_v", "p_dc_null_w", "f_null_res_gf",
        "pm_demag_margin_null_ka_m", "eval_time_s",
    ]

    write_header = not (args.resume and csv_path.exists())
    completed = 0

    with csv_path.open("a" if args.resume else "w", newline="") as f_csv:
        writer = csv.DictWriter(f_csv, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()

        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(eval_single_task, t): t for t in tasks}
            for fut in as_completed(futures):
                row = fut.result()
                writer.writerow(row)
                f_csv.flush()
                completed += 1
                wall = time.perf_counter() - t_start
                rate = completed / max(wall, 1e-6)
                print(
                    f"[{completed:3d}/{total:3d} | {rate:5.1f} eval/s] "
                    f"{row['interface_mode']:12s} | L={row['core_len_mm']:4.1f}mm | w={row['coil_w_mm']:3.1f}mm | "
                    f"Mass={row['m_dropper_active_g']:4.2f}g | Hold={row['f_hold_gf']:6.1f}gf | "
                    f"NI_null={row['ni_null_at']:5.0f}At | P_null={row['p_dc_null_w']:5.2f}W | "
                    f"Demag={row['pm_demag_margin_null_ka_m']:+6.1f}kA/m"
                )

    total_wall = time.perf_counter() - t_start
    print(f"\nSweep finished in {total_wall:.2f} s ({total / max(total_wall, 1e-6):.1f} designs/sec). Saved to {csv_path}")


if __name__ == "__main__":
    main()
