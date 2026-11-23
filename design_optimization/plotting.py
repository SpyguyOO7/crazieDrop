import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

from geometry import AxiModel
from solvers import StaticResult, TransientResult


def _overlay_parts(ax, model: AxiModel) -> None:
    color_map = {
        "pm": "#ff4d4d",
        "steel": "#d9d9d9",
        "coil": "#ff9933",
        "linear": "#66cc66",
    }
    for p in model.parts:
        ec = color_map.get(p.mat.kind, "white")
        pts_pos = np.asarray(p.pts_rz_mm)
        pts_neg = pts_pos.copy()
        pts_neg[:, 0] *= -1.0
        ax.add_patch(MplPolygon(pts_pos, closed=True, fill=False, edgecolor=ec, linewidth=1.4))
        ax.add_patch(MplPolygon(pts_neg, closed=True, fill=False, edgecolor=ec, linewidth=1.4))


def plot_field_map(
    model: AxiModel,
    res: StaticResult,
    out_path: str,
    r_max_plot_mm: float = 9.0,
    z_min_plot_mm: float = -5.0,
    z_max_plot_mm: float = 13.0,
    title: str = "",
) -> None:
    nr, nz = 120, 200
    r_mm = np.linspace(0.03, r_max_plot_mm, nr)
    z_mm = np.linspace(z_min_plot_mm, z_max_plot_mm, nz)
    rr_mm, zz_mm = np.meshgrid(r_mm, z_mm)

    b_r_grid = np.zeros_like(rr_mm)
    b_z_grid = np.zeros_like(rr_mm)
    b_mag = np.zeros_like(rr_mm)
    psi_grid = np.zeros_like(rr_mm)

    mesh = model.mesh
    for i in range(nz):
        for j in range(nr):
            mip = mesh(rr_mm[i, j] * 1e-3, zz_mm[i, j] * 1e-3)
            br, bz = res.b_vec_cf(mip)
            b_r_grid[i, j] = br
            b_z_grid[i, j] = bz
            b_mag[i, j] = np.hypot(br, bz)
            psi_grid[i, j] = res.psi_gf(mip)

    r_full = np.hstack([-r_mm[::-1], r_mm])
    b_full = np.hstack([b_mag[:, ::-1], b_mag])
    bz_full = np.hstack([b_z_grid[:, ::-1], b_z_grid])
    br_full = np.hstack([-b_r_grid[:, ::-1], b_r_grid])
    psi_full = np.hstack([psi_grid[:, ::-1], psi_grid])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 6.8), dpi=150, sharey=True)

    # Left panel: |B| magnitude + exact flux lines (contours of psi = r * A_phi)
    vmax_mag = max(1.6, float(np.percentile(b_mag, 99.5)))
    im1 = ax1.pcolormesh(r_full, z_mm, b_full, shading="auto", cmap="inferno", vmin=0.0, vmax=vmax_mag)
    cb1 = fig.colorbar(im1, ax=ax1, pad=0.02, fraction=0.046)
    cb1.set_label("|B| Magnitude (T)")

    p_max = np.max(np.abs(psi_grid))
    if p_max > 1e-12:
        levels = np.linspace(-0.92 * p_max, 0.92 * p_max, 24)
        levels = levels[np.abs(levels) > 0.025 * p_max]
        if len(levels) > 0:
            ax1.contour(r_full, z_mm, psi_full, levels=levels, colors="cyan", linewidths=0.65, alpha=0.8)

    _overlay_parts(ax1, model)
    ax1.axvline(0.0, color="white", linestyle="--", linewidth=0.5, alpha=0.4)
    ax1.set_aspect("equal")
    ax1.set_xlim(-r_max_plot_mm, r_max_plot_mm)
    ax1.set_ylim(z_min_plot_mm, z_max_plot_mm)
    ax1.set_xlabel("Radius r (mm)")
    ax1.set_ylabel("Axial Position z (mm)")
    ax1.set_title("|B| & Flux Contours")

    # Right panel: Signed axial B_z (red = +z, blue = -z) + normalized B direction arrows
    vmax_bz = max(1.2, float(np.percentile(np.abs(b_z_grid), 99.5)))
    im2 = ax2.pcolormesh(r_full, z_mm, bz_full, shading="auto", cmap="RdBu_r", vmin=-vmax_bz, vmax=vmax_bz)
    cb2 = fig.colorbar(im2, ax=ax2, pad=0.02, fraction=0.046)
    cb2.set_label("Signed Axial B_z (T)")

    step_z, step_r = 12, 12
    rr_full = np.hstack([-rr_mm[:, ::-1], rr_mm])
    zz_full = np.hstack([zz_mm[:, ::-1], zz_mm])
    sub_r = rr_full[::step_z, ::step_r]
    sub_z = zz_full[::step_z, ::step_r]
    sub_br = br_full[::step_z, ::step_r]
    sub_bz = bz_full[::step_z, ::step_r]
    sub_norm = np.hypot(sub_br, sub_bz)
    mask = sub_norm > 0.015 * vmax_mag
    ax2.quiver(
        sub_r[mask],
        sub_z[mask],
        (sub_br / np.maximum(sub_norm, 1e-9))[mask],
        (sub_bz / np.maximum(sub_norm, 1e-9))[mask],
        color="black",
        alpha=0.65,
        scale=26,
        width=0.0035,
    )

    _overlay_parts(ax2, model)
    ax2.axvline(0.0, color="black", linestyle="--", linewidth=0.5, alpha=0.35)
    ax2.set_aspect("equal")
    ax2.set_xlim(-r_max_plot_mm, r_max_plot_mm)
    ax2.set_ylim(z_min_plot_mm, z_max_plot_mm)
    ax2.set_xlabel("Radius r (mm)")
    ax2.set_title("Signed B_z & Field Direction")

    fig.suptitle(
        title or f"I = {res.coil_current_a:.2f} A | F_z = {res.force_z_n:+.3f} N | B_gap = {res.b_gap_center_t:+.3f} T",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_transient_summary(tr: TransientResult, out_path: str, title: str = "") -> None:
    t_ms = tr.t_s * 1e3
    fig, axes = plt.subplots(4, 1, figsize=(7.5, 8.5), dpi=150, sharex=True)

    axes[0].plot(t_ms, tr.i_coil_a, color="#d95f02", lw=1.8, label="I_coil (A)")
    axes[0].set_ylabel("Coil Current (A)")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc="upper right")

    axes[1].plot(t_ms, tr.v_cap_v, color="#1b9e77", lw=1.8, label="V_cap (V)")
    axes[1].set_ylabel("C1 Voltage (V)")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(loc="upper right")

    axes[2].plot(t_ms, tr.force_z_n, color="#7570b3", lw=1.8, label="F_mag_z (N)")
    axes[2].axhline(0.027 * 9.80665, color="red", ls="--", lw=1.1, label="Zippy Weight (0.265 N)")
    axes[2].axhline(0.0, color="black", ls=":", lw=0.9)
    axes[2].set_ylabel("Force on Zippy (N)")
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(loc="upper right")

    axes[3].plot(t_ms, tr.zippy_disp_mm, color="#e7298a", lw=1.8, label="Drop Distance (mm)")
    axes[3].set_ylabel("Separation (mm)")
    axes[3].set_xlabel("Time (ms)")
    axes[3].grid(True, alpha=0.3)
    axes[3].legend(loc="upper left")

    if title:
        fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
