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


# --- Tip Geometry Builder (Exact Exp 06 Winding & Target Settings) ---

def build_tip_candidate(cand_name: str, gap_mm: float = cfg.CONTACT_GAP_MM):
    match = [c for c in cfg.TIP_CANDIDATES if c[0] == cand_name]
    if not match:
        raise ValueError(f"Unknown candidate '{cand_name}'")
    _, family, prm = match[0]

    r_sh = cfg.R_SHANK_MM
    c_len = cfg.CORE_LEN_MM
    r_pm = cfg.R_STOCK_MM

    if family == "horn":
        r_tip = prm["r_tip"]
        t_len = prm["taper_len"]
        core_pts = (
            (0.0, 0.0),
            (r_tip, 0.0),
            (r_sh, t_len),
            (r_sh, c_len),
            (0.0, c_len),
        )
        contact_area_mm2 = np.pi * (r_tip ** 2)

    elif family == "flat":
        r_foot = prm["r_foot"]
        f_h = prm["foot_h"]
        if f_h <= 0.0 or abs(r_foot - r_sh) < 1e-6:
            core_pts = (
                (0.0, 0.0),
                (r_sh, 0.0),
                (r_sh, c_len),
                (0.0, c_len),
            )
        else:
            core_pts = (
                (0.0, 0.0),
                (r_foot, 0.0),
                (r_foot, f_h),
                (r_sh, f_h + 0.20),
                (r_sh, c_len),
                (0.0, c_len),
            )
        contact_area_mm2 = np.pi * (r_foot ** 2)

    elif family == "rim":
        # Annular rim contact from r_in = (r_sh - w_rim) out to r_sh!
        # Center cup is recessed upward to z = recess_h (with slight 118-deg drill-tip cone at center)
        w_rim = prm["w_rim"]
        rec_h = prm["recess_h"]
        r_in = r_sh - w_rim
        z_apex = rec_h + 0.30 * r_in  # Standard tailstock drill-point cone inside cup!
        core_pts = (
            (0.0, z_apex),
            (r_in, rec_h),
            (r_in, 0.0),
            (r_sh, 0.0),
            (r_sh, c_len),
            (0.0, c_len),
        )
        contact_area_mm2 = np.pi * (r_sh ** 2 - r_in ** 2)

    else:
        raise ValueError(f"Unknown family '{family}'")

    # Identical Exp 06 34-AWG coil on all candidates!
    win_len = c_len - 2.0 * cfg.BOBBIN_WALL_MM
    ws_raw = make_wire_spec(cfg.WIRE_AWG, fill_factor=0.58, ins_ratio=1.09)
    n_lay = int(np.floor(cfg.COIL_W_MM / ws_raw.d_ins_mm))
    tpl = int(np.floor(win_len / ws_raw.d_ins_mm))
    n_real = int(np.round(0.88 * n_lay * tpl))
    eff_fill = (n_real * ws_raw.area_cu_mm2) / (cfg.COIL_W_MM * win_len)
    w34 = make_wire_spec(cfg.WIRE_AWG, fill_factor=eff_fill, ins_ratio=1.09)

    coil_r_in = r_sh + cfg.BOBBIN_WALL_MM
    coil_z_bot = max(cfg.BOBBIN_WALL_MM, 0.55 if (family == "flat" and prm.get("foot_h", 0.0) > 0) else cfg.BOBBIN_WALL_MM)

    parts = [
        AxiPart("core_rod", core_pts, MAG_IRON_4000, maxh_mm=0.18, group="dropper"),
        rect_rz(
            "choke_coil",
            coil_r_in,
            coil_r_in + cfg.COIL_W_MM,
            coil_z_bot,
            c_len - cfg.BOBBIN_WALL_MM,
            COIL_MAT,
            maxh_mm=0.26,
            group="dropper",
            wire=w34,
            coil_sign=-1.0,
        ),
        rect_rz("zippy_iron_cap", 0.0, r_pm, -cfg.ZIPPY_CAP_H_MM - gap_mm, -gap_mm, MAG_IRON_4000, maxh_mm=0.18, group="target"),
        rect_rz(
            "zippy_pm",
            0.0,
            r_pm,
            -cfg.ZIPPY_CAP_H_MM - cfg.PM_HEIGHT_MM - gap_mm,
            -cfg.ZIPPY_CAP_H_MM - gap_mm,
            make_ndfeb("N52", +1.0),
            maxh_mm=0.18,
            group="target",
        ),
    ]

    mdl = build_axi_model(parts, air_pad_r_mm=26.0, air_pad_z_mm=26.0, air_maxh_mm=2.5)
    return mdl, family, contact_area_mm2


# --- Phase 1 Worker: Full Static F_z(I) Curve (0 .. 3.6 A) for Each Tip Geometry ---

def eval_static_curve(cand_name: str) -> list[dict]:
    mdl, family, area_mm2 = build_tip_candidate(cand_name)
    sess = StaticSolverSession(mdl, fe_order=2)
    m_drop = mdl.parts[0].mass_g + mdl.parts[1].mass_g
    coil_p = mdl.parts[1]

    rows = []
    for amp in cfg.CURRENTS_A:
        res = sess.solve(float(amp), compute_diagnostics=False)
        rows.append({
            "name": cand_name,
            "family": family,
            "contact_area_mm2": round(area_mm2, 2),
            "m_core_g": round(mdl.parts[0].mass_g, 3),
            "m_dropper_g": round(m_drop, 3),
            "r_coil_ohm": round(coil_p.coil_res_ohm, 3),
            "current_a": float(amp),
            "ni_at": round(float(amp) * coil_p.n_turns, 1),
            "force_gf": round(res.force_z_n / 9.80665 * 1000.0, 2),
            "b_gap_center_t": round(res.b_gap_center_t, 4),
        })
    return rows


# --- Phase 2 Worker: Moving-Gap Transient Voltage Sweep Across Key Tip Archetypes ---

def eval_transient_job(job: tuple[str, float]) -> dict:
    cand_name, v0 = job
    t0 = time.perf_counter()
    mdl, family, area_mm2 = build_tip_candidate(cand_name)
    m_drop = mdl.parts[0].mass_g + mdl.parts[1].mass_g

    ckt = CircuitParams(v0_cap_v=v0, c1_f=cfg.CAP_C1_F, pulse_ms=cfg.PULSE_MS)
    tr = solve_transient_pulse(mdl, ckt, t_end_ms=42.0, dt_ms=0.2, zippy_mass_g=27.0)

    idx_1mm = np.where(tr.zippy_disp_mm >= 1.0)[0]
    t_1mm = float(tr.t_s[idx_1mm[0]] * 1e3) if len(idx_1mm) > 0 else -1.0

    return {
        "name": cand_name,
        "family": family,
        "contact_area_mm2": round(area_mm2, 2),
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
    parser = argparse.ArgumentParser(description="Exp 08: Horn vs Flat vs Outer-Rim Contact Geometries (Exp 06 Settings)")
    parser.add_argument("--workers", type=int, default=cfg.DEFAULT_WORKERS)
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    t_start = time.perf_counter()

    names = [c[0] for c in cfg.TIP_CANDIDATES]
    print(f"=== Phase 1: Evaluating Static Curves for {len(names)} Tip/Horn/Rim Geometries ({args.workers} workers) ===")

    s_rows = []
    with ProcessPoolExecutor(max_workers=min(args.workers, len(names))) as pool:
        futs = {pool.submit(eval_static_curve, n): n for n in names}
        for fut in as_completed(futs):
            rows = fut.result()
            s_rows.extend(rows)
            r0 = rows[0]
            f_min = min(r["force_gf"] for r in rows)
            f_2p4a = [r["force_gf"] for r in rows if abs(r["current_a"] - 2.4) < 1e-3][0]
            print(
                f"  [Static] {r0['name']:22s} ({r0['family']:4s}, A={r0['contact_area_mm2']:5.2f}mm^2, m={r0['m_dropper_g']:.2f}g) -> "
                f"Hold(0A)={r0['force_gf']:6.1f}gf | F(2.4A)={f_2p4a:+6.1f}gf | MaxRepel={f_min:+6.1f}gf"
            )

    order_idx = {n: i for i, n in enumerate(names)}
    s_rows.sort(key=lambda r: (order_idx[r["name"]], r["current_a"]))
    with (out_dir / "static_tip_curves.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(s_rows[0].keys()))
        writer.writeheader()
        writer.writerows(s_rows)

    jobs = [(arch, float(v)) for arch in cfg.TRANSIENT_ARCHETYPES for v in cfg.TRANSIENT_VOLTAGES_V]
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
                f"  [{done:2d}/{len(jobs):2d}] {row['name']:20s} @ {row['v0_v']:4.1f}V | "
                f"Hold={row['hold_gf']:5.0f}gf | MinF={row['min_f_gf']:+6.1f}gf | "
                f"1mm@{row['t_1mm_ms']:5.1f}ms | Sep42ms={row['sep_42ms_mm']:5.2f}mm | {status}"
            )

    arch_order = {a: i for i, a in enumerate(cfg.TRANSIENT_ARCHETYPES)}
    tr_rows.sort(key=lambda r: (arch_order[r["name"]], r["v0_v"]))
    with (out_dir / "transient_tip_comparison.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(tr_rows[0].keys()))
        writer.writeheader()
        writer.writerows(tr_rows)

    print(f"\nExp 08 complete in {time.perf_counter() - t_start:.1f} s! Run plot_exp08.py next.")


if __name__ == "__main__":
    main()
