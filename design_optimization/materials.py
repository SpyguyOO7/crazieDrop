from dataclasses import dataclass
import numpy as np

MU_0 = 4.0 * np.pi * 1e-7
CU_RHO_OHM_M = 1.7241e-8
CU_DENSITY_KG_M3 = 8960.0


# --- Material Definitions ---

@dataclass(frozen=True)
class Material:
    name: str
    kind: str  # "air", "linear", "steel", "pm", "coil"
    density_kg_m3: float = 0.0
    mu_r: float = 1.0
    sigma_s_m: float = 0.0
    br_z_t: float = 0.0
    h_knee_a_m: float = 0.0
    b_sat_t: float = 0.0
    mu_r_init: float = 1.0


AIR = Material(
    name="air",
    kind="air",
    density_kg_m3=0.0,
    mu_r=1.0,
)

FR4 = Material(
    name="fr4",
    kind="linear",
    density_kg_m3=1850.0,
    mu_r=1.0,
)

STEEL_1018 = Material(
    name="steel_1018",
    kind="steel",
    density_kg_m3=7870.0,
    mu_r=1.0,
    sigma_s_m=6.2e6,
    b_sat_t=1.85,
    mu_r_init=850.0,
)

# Exact McMaster-Carr "Multipurpose Low-Carbon Magnetic Iron Rod" (1/4" D, Rockwell B45, annealed, mu_r = 4000)
# Consumable-electrode-melted pure magnetic iron (ASTM A848-type):
#   mu_r = 4000, pure-iron saturation J_sat = 2.15 T, electrical conductivity sigma ~ 1.0e7 S/m, density = 7870 kg/m^3.
MAG_IRON_4000 = Material(
    name="mcmaster_mag_iron_4000",
    kind="steel",
    density_kg_m3=7870.0,
    mu_r=1.0,
    sigma_s_m=1.0e7,
    b_sat_t=2.15,
    mu_r_init=4000.0,
)

# Standard careful lathe-faced contact gap baseline:
#   Careful lathe turning / facing typically achieves Ra ~ 0.8 - 3.2 um (peak-to-valley Rz ~ 4 - 15 um),
#   plus minor flatness/squareness tolerance across a 6.35 mm face -> ~10 - 20 um effective magnetic air gap
#   (15 um = 0.015 mm nominal baseline).
LATHE_FINISH_GAP_MM = 0.015


def make_ndfeb(grade: str = "N52", sign: float = 1.0) -> Material:
    if sign not in (-1.0, 1.0):
        raise ValueError("Magnet sign must be +1.0 (+z) or -1.0 (-z).")

    table = {
        "N35": (1.20, 1.05, 860e3),
        "N42": (1.31, 1.05, 900e3),
        "N48": (1.39, 1.05, 860e3),
        "N52": (1.46, 1.05, 820e3),
    }
    if grade not in table:
        raise ValueError(f"Unknown NdFeB grade '{grade}'. Valid grades: {list(table.keys())}")

    br_mag, mu_rec, h_knee = table[grade]
    return Material(
        name=f"ndfeb_{grade}",
        kind="pm",
        density_kg_m3=7500.0,
        mu_r=mu_rec,
        sigma_s_m=0.7e6,
        br_z_t=sign * br_mag,
        h_knee_a_m=h_knee,
    )


# --- Smooth Analytical Steel Reluctivity Law ---
# Uses smooth Fröhlich-Kennelly saturation + vacuum asymptote:
#   B(H) = mu_0 * H + J_sat * H / (H_c + |H|),  J_sat = B_sat,  H_c = J_sat / (mu_0 * (mu_r_init - 1))
# Inverting for H(B) gives an exact, strictly monotonic, C-infinity reluctivity nu(B) = H(B) / B:
#   nu(B) = [ (B - J_sat - mu_0 * H_c) + sqrt((B - J_sat - mu_0 * H_c)^2 + 4 * mu_0 * H_c * B) ] / (2 * mu_0 * B)
# Notice at B -> 0, nu(0) = 1 / (mu_0 * mu_r_init), and as B -> inf, dB/dH -> mu_0.

def steel_reluctivity_params(mat: Material) -> tuple[float, float]:
    if mat.kind != "steel":
        raise ValueError("Material must have kind='steel'.")
    j_sat = mat.b_sat_t
    h_c = j_sat / (MU_0 * (mat.mu_r_init - 1.0))
    return j_sat, h_c


def eval_steel_bh_numpy(b_vals: np.ndarray, mat: Material = STEEL_1018) -> np.ndarray:
    j_sat, h_c = steel_reluctivity_params(mat)
    b_safe = np.maximum(np.abs(b_vals), 1e-12)
    term = b_safe - j_sat - MU_0 * h_c
    h_mag = (term + np.sqrt(term * term + 4.0 * MU_0 * h_c * b_safe)) / (2.0 * MU_0)
    return np.sign(b_vals) * h_mag


# --- Magnet Wire (AWG) Calculator ---

@dataclass(frozen=True)
class WireSpec:
    awg: int
    d_bare_mm: float
    d_ins_mm: float
    area_cu_mm2: float
    ohm_per_m: float
    fill_factor: float


def make_wire_spec(awg: int, fill_factor: float = 0.58, ins_ratio: float = 1.09) -> WireSpec:
    if not (10 <= awg <= 50):
        raise ValueError(f"AWG {awg} out of valid range [10, 50].")
    if not (0.1 <= fill_factor <= 0.85):
        raise ValueError(f"fill_factor {fill_factor} out of physical range [0.1, 0.85].")

    d_bare_mm = 0.127 * (92.0 ** ((36.0 - awg) / 39.0))
    d_ins_mm = d_bare_mm * ins_ratio
    area_cu_mm2 = 0.25 * np.pi * (d_bare_mm ** 2)
    ohm_per_m = CU_RHO_OHM_M / (area_cu_mm2 * 1e-6)

    return WireSpec(
        awg=awg,
        d_bare_mm=d_bare_mm,
        d_ins_mm=d_ins_mm,
        area_cu_mm2=area_cu_mm2,
        ohm_per_m=ohm_per_m,
        fill_factor=fill_factor,
    )
