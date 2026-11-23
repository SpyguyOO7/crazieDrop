import numpy as np

from materials import LATHE_FINISH_GAP_MM


# --- Exact Exp 06 Ideal Settings Baseline ---

R_STOCK_MM = 6.35 / 2.0               # 1/4" (6.35 mm) McMaster Magnetic Iron (mu_r = 4000)
PM_HEIGHT_MM = 2.54                   # 1/10" (2.54 mm) McMaster Grade N52 on Zippy (0.603 g)
CONTACT_GAP_MM = LATHE_FINISH_GAP_MM  # 15 um bare lathe-cut contact

R_SHANK_MM = 2.20                     # Ø4.40 mm turned iron shank (matches Exp 06 ID=94)
CORE_LEN_MM = 8.00                    # 8.00 mm core length (matches Exp 06 ID=94)
ZIPPY_CAP_H_MM = 1.60                 # Ø6.35 x 1.60 mm iron cap on Zippy N52 (matches Exp 06 ID=94)

BOBBIN_WALL_MM = 0.35
COIL_W_MM = 1.80
WIRE_AWG = 34                         # 34 AWG ideal layer packing (same as Exp 06: 361 turns, R=6.70 Ohm)

CAP_C1_F = 1000e-6
PULSE_MS = 22.0


# --- Tip Geometry Families to Compare ---
# Each entry is (name, family, params_dict):
#   - "horn":     Center conical horn tapering from R_SHANK_MM (2.20mm) to r_tip over taper_len
#   - "flat":     Pure flat cylinder (no horn, r_tip == R_SHANK_MM) or flared flat foot
#   - "rim":      Outer annular rim contact! Center drilled/bored out to depth recess_h with inner radius r_in = R_SHANK - w_rim
#                 (Includes 45-deg drill-point bevel inside cup so flux transitions smoothly into outer rim!)

TIP_CANDIDATES = [
    # Family A: Center Horns across final diameters (d_tip = 1.6, 2.2, 2.8, 3.4, 3.8, 4.1 mm) & taper steepness
    ("horn_d1p6_t1p5", "horn", {"r_tip": 0.80, "taper_len": 1.50}),
    ("horn_d2p2_t1p5", "horn", {"r_tip": 1.10, "taper_len": 1.50}),
    ("horn_d2p8_t1p5", "horn", {"r_tip": 1.40, "taper_len": 1.50}),
    ("horn_d3p4_t1p5", "horn", {"r_tip": 1.70, "taper_len": 1.50}),
    ("horn_d3p8_t1p5", "horn", {"r_tip": 1.90, "taper_len": 1.50}),  # Exact Exp 06 Winner (ID=94)!
    ("horn_d3p4_steep0p7",   "horn", {"r_tip": 1.70, "taper_len": 0.70}),
    ("horn_d3p4_shallow2p4", "horn", {"r_tip": 1.70, "taper_len": 2.40}),

    # Family B: No-Horn Flat Cylinder & Wide Flat Collector Foot
    ("flat_no_horn_d4p4",    "flat", {"r_foot": 2.20, "foot_h": 0.0}),
    ("flat_flared_d5p4",     "flat", {"r_foot": 2.70, "foot_h": 0.30}),  # Fits under coil bobbin ID

    # Family C: Outer Rim / Cupped Annular Contact (r_out = 2.20mm, varying rim width w_rim & cup depth)
    ("rim_w0p35_cup0p8", "rim", {"w_rim": 0.35, "recess_h": 0.80}),  # Thin outer knife-rim (A = 4.45 mm^2 ~ d=2.4mm horn)
    ("rim_w0p60_cup0p8", "rim", {"w_rim": 0.60, "recess_h": 0.80}),  # Medium outer rim     (A = 7.16 mm^2 ~ d=3.0mm horn)
    ("rim_w0p90_cup0p8", "rim", {"w_rim": 0.90, "recess_h": 0.80}),  # Wide outer rim       (A = 9.90 mm^2 ~ d=3.55mm horn)
    ("rim_w1p25_cup0p8", "rim", {"w_rim": 1.25, "recess_h": 0.80}),  # Thick outer rim      (A = 12.4 mm^2 ~ d=4.0mm horn)
    ("rim_w0p60_deep1p6", "rim", {"w_rim": 0.60, "recess_h": 1.60}), # Deep cup check
]

CURRENTS_A = np.round(np.linspace(0.0, 3.6, 19), 3)

# Subset of 6 key archetypes for full moving-gap transient voltage-margin comparison (36 jobs total -> fast on 16 cores!)
TRANSIENT_ARCHETYPES = [
    "horn_d2p2_t1p5",
    "horn_d3p4_t1p5",
    "horn_d3p8_t1p5",
    "flat_no_horn_d4p4",
    "rim_w0p35_cup0p8",
    "rim_w0p60_cup0p8",
    "rim_w0p90_cup0p8",
]
TRANSIENT_VOLTAGES_V = [12.0, 14.0, 16.0, 18.0, 20.0, 23.0]

DEFAULT_WORKERS = 16
