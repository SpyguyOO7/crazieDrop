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


# --- Candidate Geometry Builder ---

def build_design(
    r_shank_mm: float,
    core_len_mm: float,
    r_tip_mm: float,
    coil_w_mm: float,
    cap_h_mm: float,
    gap_mm: float = cfg.CONTACT_GAP_MM,
    awg: int = cfg.WIRE_AWG,
):
    r_pm = cfg.R_STOCK_MM
    taper_len = min(1.5, 0.35 * core_len_mm)

    core_pts = (
        (0.0, 0.0),
        (r_tip_mm, 0.0),
        (r_shank_mm, taper_len),
        (r_shank_mm, core_len_mm),
        (0.0, core_len_mm),
    )

    win_len = core_len_mm - 2.0 * cfg.BOBBIN_WALL_MM
    ws_raw = make_wire_spec(awg, fill_factor=0.58, ins_ratio=1.09)
    n_lay = max(1, int(np.floor(coil_w_mm / ws_raw.d_ins_mm)))
    tpl = max(1, int(np.floor(win_len / ws_raw.d_ins_mm)))
    n_real = max(1, int(np.round(0.88 * n_lay * tpl)))
    eff_fill = (n_real * ws_raw.area_cu_mm2) / (coil_w_mm * win_len)
    w_spec = make_wire_spec(awg, fill_factor=eff_fill, ins_ratio=1.09)

    coil_r_in = r_shank_mm + cfg.BOBBIN_WALL_MM
    coil_r_out = coil_r_in + coil_w_mm

    parts = [
        AxiPart("core_rod", core_pts, MAG_IRON_4000, maxh_mm=0.20, group="dropper"),
        rect_rz(
            "choke_coil",
            coil_r_in,
            coil_r_out,
            cfg.BOBBIN_WALL_MM,
            core_len_mm - cfg.BOBBIN_WALL_MM,
            COIL_MAT,
            maxh_mm=0.26,
            group="dropper",
            wire=w_spec,
            coil_sign=-1.0,
        ),
        rect_rz("zippy_iron_cap", 0.0, r_pm, -cap_h_mm - gap_mm, -gap_mm, MAG_IRON_4000, maxh_mm=0.18, group="target"),
        rect_rz(
            "zippy_pm",
            0.0,
            r_pm,
            -cap_h_mm - cfg.PM_HEIGHT_MM - gap_mm,
            -cap_h_mm - gap_mm,
            make_ndfeb("N52", +1.0),
            maxh_mm=0.18,
            group="target",
        ),
    ]

    mdl = build_axi_model(parts, air_pad_r_mm=26.0, air_pad_z_mm=26.0, air_maxh_mm=2.5)
    return mdl, n_real


# --- Stage 1 Static Evaluator (Checks Contact Hold, Repulsion, AND 1.0mm Separation Pullback!) ---

def eval_static_candidate(task: tuple[int, float, float, float, float, float]) -> dict:
    idx, r_shank, c_len, r_tip, coil_w, cap_h = task
    t0 = time.perf_counter()

    mdl, n_real = build_design(r_shank, c_len, r_tip, coil_w, cap_h, gap_mm=cfg.CONTACT_GAP_MM)
    core_p, coil_p, cap_p, pm_p = mdl.parts

    m_dropper = core_p.mass_g + coil_p.mass_g
    m_zippy = cap_p.mass_g + pm_p.mass_g
    r_coil = coil_p.coil_res_ohm
    r_loop = r_coil + 0.12 + 0.032
    tau_rc_ms = r_loop * cfg.CAP_C1_F * 1e3

    # Peak coil current reached during 1000 uF discharge into 34 AWG coil (~0.84 factor due to longer RC vs L/R):
    i_16v = 0.84 * (16.0 / r_loop)
    i_20v = 0.84 * (20.0 / r_loop)
    i_24v = 0.84 * (24.0 / r_loop)

    sess = StaticSolverSession(mdl, fe_order=2)
    res_0 = sess.solve(0.0, compute_diagnostics=False)
    res_16 = sess.solve(i_16v, compute_diagnostics=False)
    res_20 = sess.solve(i_20v, compute_diagnostics=False)
    res_24 = sess.solve(i_24v, compute_diagnostics=True)

    f_hold_gf = res_0.force_z_n / 9.80665 * 1000.0
    f_16v_gf = res_16.force_z_n / 9.80665 * 1000.0
    f_20v_gf = res_20.force_z_n / 9.80665 * 1000.0
    f_24v_gf = res_24.force_z_n / 9.80665 * 1000.0

    # Also measure passive pullback force at gap = 1.0 mm so we penalize designs that lasso Zippy back!
    mdl_1mm, _ = build_design(r_shank, c_len, r_tip, coil_w, cap_h, gap_mm=1.00)
    res_1mm = StaticSolverSession(mdl_1mm, fe_order=2).solve(0.0, compute_diagnostics=False)
    f_pullback_1mm_gf = res_1mm.force_z_n / 9.80665 * 1000.0

    # True Escape Score: Repulsive impulse (-F_16V & -F_20V * tau_RC) minus pullback penalty!
    net_launch_score = (max(0.0, -f_16v_gf) * 1.5 + max(0.0, -f_20v_gf)) * tau_rc_ms - 8.0 * f_pullback_1mm_gf

    return {
        "id": idx,
        "r_shank_mm": r_shank,
        "core_len_mm": c_len,
        "r_tip_mm": r_tip,
        "coil_w_mm": coil_w,
        "cap_h_mm": cap_h,
        "coil_od_mm": round(2.0 * (r_shank + cfg.BOBBIN_WALL_MM + coil_w), 2),
        "n_turns_34awg": n_real,
        "r_coil_ohm": round(r_coil, 3),
        "tau_rc_ms": round(tau_rc_ms, 2),
        "m_core_g": round(core_p.mass_g, 3),
        "m_cu_g": round(coil_p.mass_g, 3),
        "m_dropper_g": round(m_dropper, 3),
        "m_zippy_g": round(m_zippy, 3),
        "f_hold_gf": round(f_hold_gf, 1),
        "f_pullback_1mm_gf": round(f_pullback_1mm_gf, 1),
        "f_16v_gf": round(f_16v_gf, 1),
        "f_20v_gf": round(f_20v_gf, 1),
        "f_24v_gf": round(f_24v_gf, 1),
        "net_launch_score": round(net_launch_score, 1),
        "demag_24v_ka_m": round(res_24.pm_demag_margin_a_m / 1e3, 1),
        "eval_s": round(time.perf_counter() - t0, 2),
    }


# --- Stage 2 Transient Verification Worker ---

def eval_transient_job(job: tuple[str, dict, float, float, float]) -> dict:
    label, cand, v0, cap_f, gap_mm = job
    t0 = time.perf_counter()
    mdl, _ = build_design(
        cand["r_shank_mm"],
        cand["core_len_mm"],
        cand["r_tip_mm"],
        cand["coil_w_mm"],
        cand["cap_h_mm"],
        gap_mm=gap_mm,
        awg=cfg.WIRE_AWG,
    )
    ckt = CircuitParams(v0_cap_v=v0, c1_f=cap_f, pulse_ms=cfg.PULSE_MS)
    tr = solve_transient_pulse(mdl, ckt, t_end_ms=42.0, dt_ms=0.2, zippy_mass_g=27.0)

    idx_1mm = np.where(tr.zippy_disp_mm >= 1.0)[0]
    t_1mm = float(tr.t_s[idx_1mm[0]] * 1e3) if len(idx_1mm) > 0 else -1.0

    return {
        "design_label": label,
        "id": int(cand["id"]),
        "m_dropper_g": cand["m_dropper_g"],
        "m_zippy_g": cand["m_zippy_g"],
        "cap_uf": int(round(cap_f * 1e6)),
        "gap_um": int(round(gap_mm * 1e3)),
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
    parser = argparse.ArgumentParser(description="Exp 06: 34 AWG Bare-Contact Iron-Cap + N52 Optimization")
    parser.add_argument("--workers", type=int, default=cfg.DEFAULT_WORKERS)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    csv_static = out_dir / "sweep_static.csv"

    done_ids = set()
    existing_rows = []
    if args.resume and csv_static.exists():
        with csv_static.open("r", newline="") as f:
            for r in csv.DictReader(f):
                done_ids.add(int(r["id"]))
                existing_rows.append({k: float(v) for k, v in r.items()})

    tasks = []
    idx = 0
    for r_sh in cfg.R_SHANK_LIST_MM:
        for c_len in cfg.CORE_LEN_LIST_MM:
            for r_tip in cfg.R_TIP_LIST_MM:
                if r_tip > r_sh:
                    continue
                for c_w in cfg.COIL_W_LIST_MM:
                    for cap_h in cfg.ZIPPY_CAP_H_LIST_MM:
                        if idx not in done_ids:
                            tasks.append((idx, r_sh, c_len, r_tip, c_w, cap_h))
                        idx += 1

    t_start = time.perf_counter()
    all_rows = list(existing_rows)

    if tasks:
        print(f"=== Phase 1: Evaluating {len(tasks)} 34-AWG Bare-Contact Designs Across {args.workers} Workers ===")
        fieldnames = [
            "id", "r_shank_mm", "core_len_mm", "r_tip_mm", "coil_w_mm", "cap_h_mm",
            "coil_od_mm", "n_turns_34awg", "r_coil_ohm", "tau_rc_ms",
            "m_core_g", "m_cu_g", "m_dropper_g", "m_zippy_g",
            "f_hold_gf", "f_pullback_1mm_gf", "f_16v_gf", "f_20v_gf", "f_24v_gf",
            "net_launch_score", "demag_24v_ka_m", "eval_s",
        ]
        write_hdr = not (args.resume and csv_static.exists())
        with csv_static.open("a" if args.resume else "w", newline="") as f_out:
            writer = csv.DictWriter(f_out, fieldnames=fieldnames)
            if write_hdr:
                writer.writeheader()
            done = 0
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                futs = {pool.submit(eval_static_candidate, t): t for t in tasks}
                for fut in as_completed(futs):
                    row = fut.result()
                    writer.writerow(row)
                    f_out.flush()
                    all_rows.append(row)
                    done += 1
                    if done % 16 == 0 or done == len(tasks):
                        rate = done / max(time.perf_counter() - t_start, 1e-6)
                        print(
                            f"  [{done:3d}/{len(tasks):3d} | {rate:4.1f} eval/s] "
                            f"m_drop={row['m_dropper_g']:.2f}g | Hold={row['f_hold_gf']:5.0f}gf | "
                            f"Pull@1mm={row['f_pullback_1mm_gf']:4.0f}gf | F_16V={row['f_16v_gf']:+6.1f}gf | "
                            f"tau={row['tau_rc_ms']:.1f}ms"
                        )

    valid = [
        r for r in all_rows
        if r["f_hold_gf"] >= cfg.MIN_HOLD_GF
        and r["m_zippy_g"] <= cfg.MAX_ZIPPY_TARGET_G
        and r["f_20v_gf"] < -40.0
    ]

    champions = []
    buckets = [
        ("Featherweight_Sub2p2g", 0.0, 2.25),
        ("Balanced_Sub2p8g",      2.25, 2.85),
        ("HeavyDuty_Sub3p6g",     2.85, 3.65),
    ]
    for b_name, m_lo, m_hi in buckets:
        sub = [r for r in valid if m_lo < r["m_dropper_g"] <= m_hi]
        if sub:
            sub.sort(key=lambda r: r["net_launch_score"], reverse=True)
            champions.append((b_name, sub[0]))

    print("\n=== Top 3 Stage-1 34-AWG Champions Selected for Moving-Gap Transient Verification ===")
    for b_name, ch in champions:
        print(
            f"  {b_name:22s} (ID={int(ch['id'])}): r_shank={ch['r_shank_mm']}mm, L={ch['core_len_mm']}mm, "
            f"r_tip={ch['r_tip_mm']}mm, w_cu={ch['coil_w_mm']}mm, cap_h={ch['cap_h_mm']}mm | "
            f"m_drop={ch['m_dropper_g']:.2f}g, R={ch['r_coil_ohm']:.2f}Ω (tau={ch['tau_rc_ms']:.1f}ms) | "
            f"Hold={ch['f_hold_gf']:.0f}gf, Pull@1mm={ch['f_pullback_1mm_gf']:.0f}gf, F_16V={ch['f_16v_gf']:+.0f}gf"
        )

    jobs = []
    # 1. Full 12V..25V voltage sweep at bare 15um contact with 1000 uF across all 3 mass-tier champions
    for b_name, ch in champions:
        for v in cfg.VERIFY_VOLTAGES_V:
            jobs.append((b_name, ch, float(v), cfg.CAP_C1_F, cfg.CONTACT_GAP_MM))

    # 2. On Balanced_Sub2p8g champion, ALSO test:
    #    - Smaller 470 uF capacitor across [16, 20, 24] V (saves ~1.2g on PCB!)
    #    - Optional Kapton bench shims (65 um & 140 um) at 16 V & 20 V to verify your shimming reasoning!
    bal_champ = champions[1][1] if len(champions) > 1 else champions[0][1]
    for v in (16.0, 20.0, 24.0):
        jobs.append(("Balanced_LightCap_470uF", bal_champ, float(v), 470e-6, cfg.CONTACT_GAP_MM))
    for g_shim in (0.065, 0.140):
        for v in (14.0, 18.0, 22.0):
            jobs.append((f"Balanced_Shim_{int(g_shim*1e3)}um", bal_champ, float(v), cfg.CAP_C1_F, g_shim))

    print(f"\n=== Phase 2: Running {len(jobs)} Moving-Gap Transient Launch Sims Across {args.workers} Workers ===")
    tr_rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(eval_transient_job, j): j for j in jobs}
        done = 0
        for fut in as_completed(futs):
            row = fut.result()
            tr_rows.append(row)
            done += 1
            status = "LAUNCHED!" if row["dropped"] else "stuck"
            print(
                f"  [{done:2d}/{len(jobs):2d}] {row['design_label']:24s} | V0={row['v0_v']:4.1f}V | "
                f"Hold={row['hold_gf']:5.0f}gf | PeakI={row['peak_i_a']:4.2f}A | "
                f"MinF={row['min_f_gf']:+6.1f}gf | Sep42ms={row['sep_42ms_mm']:5.2f}mm | {status}"
            )

    tr_rows.sort(key=lambda r: (r["design_label"], r["v0_v"]))
    csv_tr = out_dir / "transient_champions.csv"
    with csv_tr.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(tr_rows[0].keys()))
        writer.writeheader()
        writer.writerows(tr_rows)

    print(f"\nExp 06 complete in {time.perf_counter() - t_start:.1f} s! Run plot_exp06.py next.")


if __name__ == "__main__":
    main()
