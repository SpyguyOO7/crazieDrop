from dataclasses import dataclass
import numpy as np
from shapely.geometry import Point as ShpPoint
import ngsolve as ngs

from materials import MU_0, steel_reluctivity_params
from geometry import AxiModel, build_axi_model


# --- Axisymmetric Formulation Notes (Mistake-Proofing) ---
# Coordinates in NGSolve 2D mesh: x = r (meters, r >= 0), y = z (meters).
# Unknown field: u = psi(r, z) = r * A_phi(r, z)  [Wb / rad], with Dirichlet u = 0 on "axis" (r=0) and "outer".
# Magnetic flux density B = (B_r, B_z):
#   B_r = -(1/r) * du/dz
#   B_z = +(1/r) * du/dr
# Notice:
#   curl(H)_phi * (2*pi*r dr dz) integrates by parts cleanly to:
#   int_Omega [ nu(|B|) * (1/r) * (du/dr * dv/dr + du/dz * dv/dz) ] * 2*pi dr dz
# Permanent magnet (axial M_0 = Br_z / (mu_0 * mu_rec)):
#   H = nu_rec * B - H_c_vec,  where H_c_z = Br_z / (mu_0 * mu_rec).
#   Weak form term: - int_PM [ H_c_z * B_z(v) ] * 2*pi*r dr dz = - int_PM [ (Br_z / (mu_0 * mu_rec)) * (dv/dr) ] * 2*pi dr dz
# Azimuthal coil current density J_phi = coil_sign * N * I / Area_window:
#   Weak form RHS: int_Coil [ J_phi * A_phi(v) ] * 2*pi*r dr dz = int_Coil [ J_phi * (v / r) ] * 2*pi*r dr dz = int_Coil [ 2*pi * J_phi * v ] dr dz
# Eddy currents in conductive regions: J_eddy = -sigma * dA_phi/dt = -sigma * (1/r) * du/dt
#   Weak form LHS: + int_Cond [ 2*pi * (sigma / r) * (du/dt) * v ] dr dz
# Notice how ALL 1/r factors in source and magnet terms cancel analytically at r=0, and u=0 on r=0 makes (1/r)*grad(u) regular!


@dataclass
class StaticResult:
    coil_current_a: float
    force_z_n: float
    b_gap_center_t: float
    b_core_max_t: float
    min_pm_b_along_mag_t: float
    pm_demag_margin_a_m: float
    coil_flux_linkage_wb: float
    magnetic_energy_j: float
    psi_gf: ngs.GridFunction
    b_vec_cf: ngs.CoefficientFunction


def _build_region_cfs(model: AxiModel, coil_current_a: float):
    mesh = model.mesh
    r_safe = ngs.IfPos(ngs.x - 1e-9, ngs.x, 1e-9)

    # Default (air) values on every mesh material label
    mat_names = mesh.GetMaterials()
    nu_lin_dict = {m: 1.0 / MU_0 for m in mat_names}
    is_steel_dict = {m: 0.0 for m in mat_names}
    j_sat_dict = {m: 1.85 for m in mat_names}
    h_c_steel_dict = {m: 1000.0 for m in mat_names}
    pm_hc_z_dict = {m: 0.0 for m in mat_names}
    j_phi_unit_dict = {m: 0.0 for m in mat_names}
    sigma_dict = {m: 0.0 for m in mat_names}

    for p in model.parts:
        if p.name not in mat_names:
            continue
        m = p.mat
        sigma_dict[p.name] = m.sigma_s_m

        if m.kind in ("air", "linear", "coil"):
            nu_lin_dict[p.name] = 1.0 / (MU_0 * m.mu_r)
        elif m.kind == "pm":
            nu_lin_dict[p.name] = 1.0 / (MU_0 * m.mu_r)
            pm_hc_z_dict[p.name] = m.br_z_t / (MU_0 * m.mu_r)
        elif m.kind == "steel":
            is_steel_dict[p.name] = 1.0
            js, hc = steel_reluctivity_params(m)
            j_sat_dict[p.name] = js
            h_c_steel_dict[p.name] = hc
            nu_lin_dict[p.name] = 1.0 / (MU_0 * m.mu_r_init)

        if m.kind == "coil":
            area_m2 = p.area_mm2 * 1e-6
            j_phi_unit_dict[p.name] = p.coil_sign * p.n_turns / area_m2

    nu_lin_cf = mesh.MaterialCF(nu_lin_dict)
    is_steel_cf = mesh.MaterialCF(is_steel_dict)
    j_sat_cf = mesh.MaterialCF(j_sat_dict)
    hc_steel_cf = mesh.MaterialCF(h_c_steel_dict)
    pm_hcz_cf = mesh.MaterialCF(pm_hc_z_dict)
    j_phi_unit_cf = mesh.MaterialCF(j_phi_unit_dict)
    sigma_cf = mesh.MaterialCF(sigma_dict)

    return r_safe, nu_lin_cf, is_steel_cf, j_sat_cf, hc_steel_cf, pm_hcz_cf, j_phi_unit_cf, sigma_cf


def _nu_steel_cf(b_sq: ngs.CoefficientFunction, j_sat_cf, hc_steel_cf) -> ngs.CoefficientFunction:
    # Smooth C-infinity reluctivity nu(|B|) from analytical inverse of Frohlich-Kennelly + mu_0*H
    b_mag = ngs.sqrt(b_sq + 1e-16)
    term = b_mag - j_sat_cf - MU_0 * hc_steel_cf
    h_mag = (term + ngs.sqrt(term * term + 4.0 * MU_0 * hc_steel_cf * b_mag)) / (2.0 * MU_0)
    return h_mag / b_mag


def _eggshell_weight_cfs(model: AxiModel):
    assert model.mst_shell_mm is not None
    r_in, r_out, z_t_min, z_bot, z_t_max, z_top = [v * 1e-3 for v in model.mst_shell_mm]

    def _clamp01(expr):
        c_hi = ngs.IfPos(expr - 1.0, 1.0, expr)
        return ngs.IfPos(c_hi, c_hi, 0.0)

    s_r = (r_out - ngs.x) / (r_out - r_in)
    s_top = (z_top - ngs.y) / (z_top - z_t_max)
    s_bot = (ngs.y - z_bot) / (z_t_min - z_bot)

    w_r = _clamp01(s_r)
    w_top = _clamp01(s_top)
    w_bot = _clamp01(s_bot)

    dw_r_dr = ngs.IfPos(s_r, ngs.IfPos(1.0 - s_r, -1.0 / (r_out - r_in), 0.0), 0.0)
    dw_top_dz = ngs.IfPos(s_top, ngs.IfPos(1.0 - s_top, -1.0 / (z_top - z_t_max), 0.0), 0.0)
    dw_bot_dz = ngs.IfPos(s_bot, ngs.IfPos(1.0 - s_bot, +1.0 / (z_t_min - z_bot), 0.0), 0.0)

    # Product w(r, z) = w_r(r) * w_top(z) * w_bot(z):
    #   = 1 everywhere inside target [0..r_in] x [z_t_min..z_t_max]
    #   = 0 everywhere on & outside mst_box [0..r_out] x [z_bot..z_top]
    #   smoothly bilinear across the "air_mst" shell!
    w_cf = w_r * w_top * w_bot
    dw_dr_cf = dw_r_dr * w_top * w_bot
    dw_dz_cf = w_r * (dw_top_dz * w_bot + w_top * dw_bot_dz)
    return w_cf, dw_dr_cf, dw_dz_cf


def compute_mst_force_z(model: AxiModel, b_vec_cf: ngs.CoefficientFunction) -> float:
    # Volume-weighted Eggshell / Arkkio Maxwell Stress Tensor over the 2D "air_mst" shell:
    # By divergence theorem with w = 1 on target boundary and w = 0 on outer mst_box boundary:
    #   F_z = int_target T_iz n_i dS = - int_{air_mst} (T_rz * dw/dr + T_zz * dw/dz) * 2*pi*r dr dz
    if model.mst_shell_mm is None:
        return 0.0

    mesh = model.mesh
    _, dw_dr_cf, dw_dz_cf = _eggshell_weight_cfs(model)

    b_r = b_vec_cf[0]
    b_z = b_vec_cf[1]
    t_zz = (b_z * b_z - b_r * b_r) / (2.0 * MU_0)
    t_rz = (b_r * b_z) / MU_0

    integrand = -2.0 * np.pi * ngs.x * (t_rz * dw_dr_cf + t_zz * dw_dz_cf)
    return float(ngs.Integrate(integrand, mesh, definedon=mesh.Materials("air_mst"), order=6))


class StaticSolverSession:
    """Reusable FE space, bilinear form, and vectors on a fixed mesh for fast multi-current sweeps."""

    def __init__(self, model: AxiModel, fe_order: int = 2):
        self.model = model
        self.mesh = model.mesh
        self.fes = ngs.H1(self.mesh, order=fe_order, dirichlet="axis|outer")
        u, v = self.fes.TnT()

        (
            self.r_safe,
            self.nu_lin_cf,
            self.is_steel_cf,
            self.j_sat_cf,
            self.hc_steel_cf,
            self.pm_hcz_cf,
            self.j_phi_unit_cf,
            _,
        ) = _build_region_cfs(model, 0.0)

        grad_u = ngs.grad(u)
        grad_v = ngs.grad(v)
        b_sq = (grad_u * grad_u) / (self.r_safe * self.r_safe)
        nu_nl = _nu_steel_cf(b_sq, self.j_sat_cf, self.hc_steel_cf)
        nu_eff = (1.0 - self.is_steel_cf) * self.nu_lin_cf + self.is_steel_cf * nu_nl

        self.i_param = ngs.Parameter(0.0)
        self.a = ngs.BilinearForm(self.fes)
        self.a += (2.0 * np.pi / self.r_safe) * nu_eff * (grad_u * grad_v) * ngs.dx
        self.a += -(2.0 * np.pi) * self.pm_hcz_cf * grad_v[0] * ngs.dx
        self.a += -(2.0 * np.pi) * (self.i_param * self.j_phi_unit_cf) * v * ngs.dx

        self.psi = ngs.GridFunction(self.fes)
        self.res = self.psi.vec.CreateVector()
        self.res_free = self.psi.vec.CreateVector()
        self.du = self.psi.vec.CreateVector()
        self.old_vec = self.psi.vec.CreateVector()
        self.free_dofs = self.fes.FreeDofs()
        self._proj = ngs.Projector(self.free_dofs, True)

        # Pre-compute interior diagnostic probe points ONCE per geometry
        self._steel_mips = []
        self._pm_probes = []
        for p in model.parts:
            if p.mat.kind not in ("steel", "pm"):
                continue
            pts_arr = np.asarray(p.pts_rz_mm)
            r_lo, z_lo = pts_arr.min(axis=0) * 1e-3
            r_hi, z_hi = pts_arr.max(axis=0) * 1e-3
            rs = np.linspace(max(2e-5, r_lo + 2e-5), r_hi - 2e-5, 9)
            zs = np.linspace(z_lo + 2e-5, z_hi - 2e-5, 11)
            mips = [
                self.mesh(rr, zz)
                for rr in rs
                for zz in zs
                if p.poly_mm.contains_properly(ShpPoint(rr * 1e3, zz * 1e3))
            ]
            if p.mat.kind == "steel":
                self._steel_mips.extend(mips)
            else:
                self._pm_probes.append((p.mat, mips))

        self._gap_mip = None
        if model.mst_shell_mm is not None:
            z_mid_m = 0.5 * (model.mst_shell_mm[4] + model.mst_shell_mm[5]) * 1e-3
            self._gap_mip = self.mesh(2e-5, z_mid_m)

    def solve(
        self,
        coil_current_a: float = 0.0,
        newton_tol: float = 1e-8,
        newton_maxit: int = 30,
        compute_diagnostics: bool = True,
    ) -> StaticResult:
        self.i_param.Set(coil_current_a)

        for _it in range(newton_maxit):
            self.a.Apply(self.psi.vec, self.res)
            self.res_free.data = self.res
            self._proj.Project(self.res_free)
            norm_res = ngs.Norm(self.res_free)
            if norm_res < newton_tol:
                break

            self.a.AssembleLinearization(self.psi.vec)
            inv = self.a.mat.Inverse(self.free_dofs, inverse="sparsecholesky")
            self.du.data = inv * self.res

            alpha = 1.0
            self.old_vec.data = self.psi.vec
            while alpha > 1e-3:
                self.psi.vec.data = self.old_vec - alpha * self.du
                self.a.Apply(self.psi.vec, self.res)
                self.res_free.data = self.res
                self._proj.Project(self.res_free)
                if ngs.Norm(self.res_free) < (1.0 - 0.1 * alpha) * norm_res or norm_res < 1e-7:
                    break
                alpha *= 0.5

        psi_copy = ngs.GridFunction(self.fes)
        psi_copy.vec.data = self.psi.vec
        g_psi = ngs.grad(psi_copy)
        b_vec_cf = ngs.CF((-g_psi[1] / self.r_safe, g_psi[0] / self.r_safe))

        force_z = compute_mst_force_z(self.model, b_vec_cf)
        b_gap_center = float(b_vec_cf[1](self._gap_mip)) if self._gap_mip is not None else 0.0

        b_core_max = 0.0
        min_pm_b = 0.0
        min_demag_margin = 0.0
        coil_lambda = 0.0
        mag_energy = 0.0

        if compute_diagnostics:
            b_mag_cf = ngs.Norm(b_vec_cf)
            if self._steel_mips:
                b_core_max = max(float(b_mag_cf(m)) for m in self._steel_mips)

            if self._pm_probes:
                min_b_val = np.inf
                min_m_val = np.inf
                for mat, mips in self._pm_probes:
                    sgn = np.sign(mat.br_z_t)
                    br_abs = abs(mat.br_z_t)
                    denom = MU_0 * mat.mu_r
                    for m in mips:
                        b_aligned = sgn * float(b_vec_cf[1](m))
                        h_par = (b_aligned - br_abs) / denom
                        margin = mat.h_knee_a_m + h_par
                        if b_aligned < min_b_val:
                            min_b_val = b_aligned
                        if margin < min_m_val:
                            min_m_val = margin
                if min_b_val < np.inf:
                    min_pm_b = min_b_val
                    min_demag_margin = min_m_val

            coil_lambda = float(ngs.Integrate(2.0 * np.pi * self.j_phi_unit_cf * psi_copy, self.mesh, order=4))

            b_mag_gf = ngs.sqrt(b_vec_cf * b_vec_cf + 1e-16)
            term_gf = b_mag_gf - self.j_sat_cf - MU_0 * self.hc_steel_cf
            h_steel_gf = (term_gf + ngs.sqrt(term_gf * term_gf + 4.0 * MU_0 * self.hc_steel_cf * b_mag_gf)) / (2.0 * MU_0)
            w_co_steel_gf = (
                0.5 * MU_0 * h_steel_gf * h_steel_gf
                + self.j_sat_cf * (h_steel_gf - self.hc_steel_cf * ngs.log(1.0 + h_steel_gf / self.hc_steel_cf))
            )
            w_steel_gf = b_mag_gf * h_steel_gf - w_co_steel_gf
            w_lin_gf = 0.5 * self.nu_lin_cf * (b_vec_cf * b_vec_cf)
            w_dens = (1.0 - self.is_steel_cf) * w_lin_gf + self.is_steel_cf * w_steel_gf
            pot_integrand = (
                2.0 * np.pi * self.r_safe * (w_dens - self.pm_hcz_cf * b_vec_cf[1])
                - 2.0 * np.pi * (coil_current_a * self.j_phi_unit_cf) * psi_copy
            )
            mag_energy = float(ngs.Integrate(pot_integrand, self.mesh, order=5))

        return StaticResult(
            coil_current_a=coil_current_a,
            force_z_n=force_z,
            b_gap_center_t=b_gap_center,
            b_core_max_t=b_core_max,
            min_pm_b_along_mag_t=min_pm_b,
            pm_demag_margin_a_m=min_demag_margin,
            coil_flux_linkage_wb=coil_lambda,
            magnetic_energy_j=mag_energy,
            psi_gf=psi_copy,
            b_vec_cf=b_vec_cf,
        )


def solve_static(
    model: AxiModel,
    coil_current_a: float = 0.0,
    fe_order: int = 2,
    newton_tol: float = 1e-9,
    newton_maxit: int = 40,
    init_psi: ngs.GridFunction | None = None,
    compute_diagnostics: bool = True,
) -> StaticResult:
    sess = StaticSolverSession(model, fe_order=fe_order)
    if init_psi is not None:
        sess.psi.Set(init_psi)
    return sess.solve(
        coil_current_a=coil_current_a,
        newton_tol=newton_tol,
        newton_maxit=newton_maxit,
        compute_diagnostics=compute_diagnostics,
    )


# --- Tier 2: Coupled Transient Eddy-Current + PCB Capacitor/MOSFET/Flyback Circuit Solver ---

@dataclass
class CircuitParams:
    v0_cap_v: float = 15.0        # Initial voltage on C1 (set via R5 trim pot, 3.0V .. 27.2V)
    c1_f: float = 1000e-6         # Bulk capacitor C1 (1000 uF on crazieDrop PCB)
    c1_esr_ohm: float = 0.12      # Typical radial electrolytic ESR
    r_dson_ohm: float = 0.032     # AO3400A NMOS on-resistance at Vgs ~ 3.3V
    v_diode_drop_v: float = 0.38  # SS14 flyback diode forward drop
    pulse_ms: float = 25.0        # Digital GPIO pulse width (also bounded by R1*C2 = 220ms hardware one-shot)
    v_bat_min_v: float = 3.0      # D3 bypass diode clamps C1 from discharging below ~V_bat - Vf


@dataclass
class TransientResult:
    t_s: np.ndarray
    i_coil_a: np.ndarray
    v_cap_v: np.ndarray
    force_z_n: np.ndarray
    b_gap_center_t: np.ndarray
    zippy_disp_mm: np.ndarray
    zippy_vel_m_s: np.ndarray
    dropped: bool


def solve_transient_pulse(
    model: AxiModel,
    circuit: CircuitParams,
    t_end_ms: float = 40.0,
    dt_ms: float = 0.2,
    zippy_mass_g: float = 27.0,
    fe_order: int = 2,
    remesh_step_mm: float = 0.08,
) -> TransientResult:
    # Step 0: Steady-state passive initial condition (I = 0, no eddy currents)
    init_res = solve_static(model, coil_current_a=0.0, fe_order=fe_order)

    r_coil = sum(p.coil_res_ohm for p in model.parts if p.mat.kind == "coil")
    if r_coil <= 0.0:
        raise ValueError("Transient circuit solver requires a coil with a WireSpec (non-zero resistance).")

    dt = dt_ms * 1e-3
    n_steps = int(np.round((t_end_ms * 1e-3) / dt))
    base_parts = list(model.parts)

    def _setup_fe_state(curr_model: AxiModel, old_psi_gf: ngs.GridFunction | None):
        mesh = curr_model.mesh
        fes = ngs.H1(mesh, order=fe_order, dirichlet="axis|outer")
        u, v = fes.TnT()
        cfs = _build_region_cfs(curr_model, 0.0)
        r_safe, nu_lin_cf, is_steel_cf, j_sat_cf, hc_steel_cf, pm_hcz_cf, j_phi_unit_cf, sigma_cf = cfs

        psi = ngs.GridFunction(fes)
        psi_prev = ngs.GridFunction(fes)
        if old_psi_gf is not None:
            # Project magnetic flux potential smoothly onto the new mesh (dropper steel/coil haven't moved!)
            psi.Set(old_psi_gf)
            psi_prev.vec.data = psi.vec

        grad_u = ngs.grad(u)
        grad_v = ngs.grad(v)
        b_sq = (grad_u * grad_u) / (r_safe * r_safe)
        nu_nl = _nu_steel_cf(b_sq, j_sat_cf, hc_steel_cf)
        nu_eff = (1.0 - is_steel_cf) * nu_lin_cf + is_steel_cf * nu_nl

        i_param = ngs.Parameter(0.0)
        a_step = ngs.BilinearForm(fes)
        a_step += (2.0 * np.pi / r_safe) * nu_eff * (grad_u * grad_v) * ngs.dx
        a_step += (2.0 * np.pi * sigma_cf / (r_safe * dt)) * (u - psi_prev) * v * ngs.dx
        a_step += -(2.0 * np.pi) * pm_hcz_cf * grad_v[0] * ngs.dx
        a_step += -(2.0 * np.pi) * (i_param * j_phi_unit_cf) * v * ngs.dx

        return mesh, fes, r_safe, j_phi_unit_cf, psi, psi_prev, i_param, a_step

    curr_model = model
    mesh, fes, r_safe, j_phi_unit_cf, psi, psi_prev, i_param, a_step = _setup_fe_state(curr_model, init_res.psi_gf)
    meshed_disp_mm = 0.0

    t_hist = np.zeros(n_steps + 1)
    i_hist = np.zeros(n_steps + 1)
    v_hist = np.zeros(n_steps + 1)
    f_hist = np.zeros(n_steps + 1)
    b_hist = np.zeros(n_steps + 1)
    z_disp_mm = np.zeros(n_steps + 1)
    z_vel = np.zeros(n_steps + 1)

    v_hist[0] = circuit.v0_cap_v
    f_hist[0] = init_res.force_z_n
    b_hist[0] = init_res.b_gap_center_t

    lambda_prev = init_res.coil_flux_linkage_wb
    v_cap = circuit.v0_cap_v
    i_curr = 0.0

    r_on_tot = r_coil + circuit.c1_esr_ohm + circuit.r_dson_ohm
    m_zippy_kg = zippy_mass_g * 1e-3
    g_acc = 9.80665

    res_vec = psi.vec.CreateVector()
    res_free = psi.vec.CreateVector()
    du_vec = psi.vec.CreateVector()
    free_dofs = fes.FreeDofs()

    for k in range(1, n_steps + 1):
        t_now = k * dt
        t_hist[k] = t_now
        q1_on = (t_now <= circuit.pulse_ms * 1e-3 + 1e-12)

        # Rebuild mesh at Zippy's true instantaneous separation when target moves significantly
        curr_disp = z_disp_mm[k - 1]
        thresh = min(remesh_step_mm, max(0.02, 0.20 * curr_disp))
        if abs(curr_disp - meshed_disp_mm) >= thresh:
            shifted_parts = [
                p.with_z_shift(-curr_disp) if p.group == "target" else p
                for p in base_parts
            ]
            curr_model = build_axi_model(
                shifted_parts,
                air_pad_r_mm=25.0,
                air_pad_z_mm=25.0,
                air_maxh_mm=2.5,
            )
            mesh, fes, r_safe, j_phi_unit_cf, psi, psi_prev, i_param, a_step = _setup_fe_state(curr_model, psi_prev)
            lambda_prev = float(ngs.Integrate(2.0 * np.pi * j_phi_unit_cf * psi_prev, mesh, order=4))
            res_vec = psi.vec.CreateVector()
            res_free = psi.vec.CreateVector()
            du_vec = psi.vec.CreateVector()
            free_dofs = fes.FreeDofs()
            meshed_disp_mm = curr_disp

        # Coupled circuit-field iteration (outer updates with inner Newton steps on free DOFs)
        for _fp in range(6):
            i_param.Set(i_curr)
            for _nt in range(12):
                a_step.Apply(psi.vec, res_vec)
                res_free.data = res_vec
                for idx in range(fes.ndof):
                    if not free_dofs[idx]:
                        res_free[idx] = 0.0
                if ngs.Norm(res_free) < 1e-8:
                    break
                a_step.AssembleLinearization(psi.vec)
                inv = a_step.mat.Inverse(free_dofs, inverse="sparsecholesky")
                du_vec.data = inv * res_vec
                psi.vec.data -= du_vec

            lambda_k = float(ngs.Integrate(2.0 * np.pi * j_phi_unit_cf * psi, mesh, order=4))
            dlambda_dt = (lambda_k - lambda_prev) / dt

            if q1_on:
                i_target = (v_cap - dlambda_dt) / (r_on_tot + dt / circuit.c1_f)
                i_target = max(0.0, i_target)
            else:
                if i_curr > 0.0 or (-dlambda_dt > circuit.v_diode_drop_v):
                    i_target = max(0.0, (-circuit.v_diode_drop_v - dlambda_dt) / r_coil)
                else:
                    i_target = 0.0

            if abs(i_target - i_curr) < 1e-4:
                i_curr = i_target
                break
            i_curr = 0.6 * i_target + 0.4 * i_curr

        if q1_on:
            v_cap = max(circuit.v_bat_min_v, v_cap - (dt / circuit.c1_f) * i_curr)

        lambda_prev = float(ngs.Integrate(2.0 * np.pi * j_phi_unit_cf * psi, mesh, order=4))
        psi_prev.vec.data = psi.vec

        g_psi = ngs.grad(psi)
        b_vec_cf = ngs.CF((-g_psi[1] / r_safe, g_psi[0] / r_safe))
        fz = compute_mst_force_z(curr_model, b_vec_cf)

        b_gap = 0.0
        if curr_model.mst_box_mm is not None:
            b_gap = float(b_vec_cf[1](mesh(2e-5, curr_model.mst_box_mm[3] * 1e-3)))

        i_hist[k] = i_curr
        v_hist[k] = v_cap
        f_hist[k] = fz
        b_hist[k] = b_gap

        # True 1D rigid-body separation kinematics driven by actual FEM force fz at true gap!
        # (+z is UP toward dropper; disp_mm >= 0 is downward separation from initial contact stop.)
        acc_down = g_acc - (fz / m_zippy_kg)
        if z_disp_mm[k - 1] <= 0.0 and acc_down <= 0.0:
            z_vel[k] = 0.0
            z_disp_mm[k] = 0.0
        else:
            v_new = z_vel[k - 1] + acc_down * dt
            x_new_mm = z_disp_mm[k - 1] + v_new * dt * 1e3
            if x_new_mm <= 0.0:
                # Re-captured / slammed back into contact stop!
                z_disp_mm[k] = 0.0
                z_vel[k] = 0.0
            else:
                z_disp_mm[k] = x_new_mm
                z_vel[k] = v_new

    dropped = bool(z_disp_mm[-1] > 1.5 and z_vel[-1] > 0.0)

    return TransientResult(
        t_s=t_hist,
        i_coil_a=i_hist,
        v_cap_v=v_hist,
        force_z_n=f_hist,
        b_gap_center_t=b_hist,
        zippy_disp_mm=z_disp_mm,
        zippy_vel_m_s=z_vel,
        dropped=dropped,
    )
