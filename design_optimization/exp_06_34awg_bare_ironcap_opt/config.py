import numpy as np

from materials import LATHE_FINISH_GAP_MM


# --- Fixed Stock, Bare-Metal Lathe Contact, and 34 AWG Winding Constraints ---

R_STOCK_MM = 6.35 / 2.0               # 1/4" (6.35 mm) McMaster Magnetic Iron (mu_r = 4000)
PM_HEIGHT_MM = 2.54                   # 1/10" (2.54 mm) McMaster Grade N52 on Zippy (0.603 g)
CONTACT_GAP_MM = LATHE_FINISH_GAP_MM  # 15 um bare lathe-cut metal-to-metal contact (NO plastic spacer!)

BOBBIN_WALL_MM = 0.35
WIRE_AWG = 34                         # 34 AWG magnet wire (0.160 mm bare, already on hand!)

CAP_C1_F = 1000e-6
PULSE_MS = 22.0


# --- Stage 1 Multi-Variable Bare-Metal Sweep Grid (288 Lathe-Turnable Candidates) ---

R_SHANK_LIST_MM = [2.20, 2.70, 3.175]     # Turned-down 4.4mm / 5.4mm shank vs untouched 6.35mm stock
CORE_LEN_LIST_MM = [5.0, 6.0, 7.0, 8.0]   # Dropper core total length
R_TIP_LIST_MM = [1.25, 1.55, 1.90]        # Concentrated horn tip radius
COIL_W_LIST_MM = [1.30, 1.80]             # 34 AWG radial build thickness
ZIPPY_CAP_H_LIST_MM = [1.00, 1.30, 1.60, 1.90]  # Thicker iron-cap range to tame 1mm pullback!

MAX_ZIPPY_TARGET_G = 1.15
MIN_HOLD_GF = 300.0

# Stage 2 Full Moving-Gap Transient Voltage-Window Verification:
VERIFY_VOLTAGES_V = [12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 25.0]

# Bonus Checks on the Overall Winner:
#   1. Lighter Capacitor Check (470 uF & 680 uF vs 1000 uF) to see if PCB C1 mass can drop by >1g!
#   2. Optional Bench Kapton Shim Check (15 um bare vs 50 um [2 mil] vs 125 um [5 mil] Kapton)
SHIM_CHECK_GAPS_MM = [0.015, 0.065, 0.140]

DEFAULT_WORKERS = 16
