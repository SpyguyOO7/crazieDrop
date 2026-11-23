import numpy as np

from materials import LATHE_FINISH_GAP_MM


# --- Fixed Stock, Magnet, and Existing 32 AWG Driver Constraints ---

R_STOCK_MM = 6.35 / 2.0        # 1/4" (6.35 mm) McMaster Magnetic Iron (mu_r = 4000)
PM_HEIGHT_MM = 2.54            # 1/10" (2.54 mm) McMaster Grade N52 on Zippy (0.603 g)
CONTACT_GAP_MM = LATHE_FINISH_GAP_MM  # 15 um careful lathe finish

BOBBIN_WALL_MM = 0.35
WIRE_AWG = 32                  # Keep your existing 32 AWG magnet wire!

CAP_C1_F = 1000e-6             # Existing 1000 uF capacitor
PULSE_MS = 18.0


# --- Stage 1 Multi-Variable Geometry Sweep Grid (360 Lathe-Turnable Candidates) ---
# Architecture:
#   - Dropper (Drone Side, NO Permanent Magnet!):
#       Turned from 1/4" Mag Iron rod -> Shank radius R_SHANK_MM, length CORE_LEN_MM,
#       conical bottom horn tapering over TAPER_LEN_MM to contact radius R_TIP_MM,
#       wound with 32 AWG coil of radial thickness COIL_W_MM.
#   - Target (Zippy Side):
#       Thin Mag Iron Cap of thickness CAP_H_MM (6.35 mm base, top tapered to R_CAP_TOP = max(R_TIP + 0.5mm, R_STOCK) or flat 6.35mm)
#       glued directly on top of the 1/4" x 1/10" N52 Magnet!

R_SHANK_LIST_MM = [2.20, 2.70, 3.175]     # Turned-down 4.4mm / 5.4mm shank vs untouched 6.35mm stock!
CORE_LEN_LIST_MM = [4.5, 5.5, 6.5, 7.5]   # Dropper core total length
R_TIP_LIST_MM = [1.20, 1.55, 1.95]        # Concentrated horn tip radius (d_tip = 2.4, 3.1, 3.9 mm)
COIL_W_LIST_MM = [1.20, 1.70]             # 32 AWG radial build thickness
ZIPPY_CAP_H_LIST_MM = [0.60, 0.90, 1.20, 1.50, 2.00]  # Zippy front iron-cap thickness!

# Mass & Performance Targets:
MAX_ZIPPY_TARGET_G = 1.10     # Keep Zippy's [Iron Cap + N52] light (~0.75 - 1.05 g)
MIN_HOLD_GF = 250.0           # User requirement: keep strong >= ~280-300 gf flight hold!

# Phase 2 Transient Verification Voltages on Top Pareto Winners:
VERIFY_VOLTAGES_V = [12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 25.0]

DEFAULT_WORKERS = 16
