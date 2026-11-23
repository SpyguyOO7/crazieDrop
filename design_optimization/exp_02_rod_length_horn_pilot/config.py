import numpy as np

from materials import LATHE_FINISH_GAP_MM


# --- Fixed Hardware Stock & Magnet Parameters ---

R_STOCK_MM = 6.35 / 2.0        # 1/4" (6.35 mm) McMaster Magnetic Iron (mu_r = 4000)
PM_HEIGHT_MM = 2.54            # 1/10" (2.54 mm) McMaster Grade N52 disk magnet
PM_GRADE = "N52"
CONTACT_GAP_MM = LATHE_FINISH_GAP_MM  # 15 um careful lathe-cut finish


# --- Lathe Collector + Bottom Horn / Concentrated Funnel Geometry ---

TOP_COLLECTOR_LEN_MM = 1.0     # Top collector pad right below the 1/4" N52 magnet
R_SHANK_MM = R_STOCK_MM        # Straight shank radius (kept at 1/4" stock in Pilot Exp 02)
HORN_TAPER_LEN_MM = 1.5        # Conical bottom taper axial length
R_TIP_MM = 1.50                # Concentrated contact land radius (3.0 mm tip diameter, matching Exp 01)

BOBBIN_WALL_MM = 0.35          # 3D-printed bobbin radial wall + axial end-cheek thickness
WIRE_AWG = 32
WIRE_FILL = 0.58


# --- Zippy Target Parameters ---

ZIPPY_LEN_MM = 3.0
ZIPPY_TAPER_LEN_MM = 0.60      # Horn bevel height on Zippy when target_mode == "tapered_horn"


# --- Sweep Grids ---

# Core total axial length from 3.5 mm (minimum leaving ~0.3 mm coil slot) up to 14.0 mm (9.5 mm is baseline)
CORE_LENGTHS_MM = np.round(np.arange(3.5, 14.01, 0.5), 2)

# Coil radial build thickness (sweeps copper mass across the full 1.8g -> 6.5g total active mass range)
COIL_WIDTHS_MM = np.array([0.8, 1.2, 1.8, 2.4])

# Compare three pole-interface styles across length & mass:
#   1. "double_horn":  Dropper tapered to r=1.5mm + Zippy tapered to r=1.5mm (concentrated funnel throat)
#   2. "dropper_horn": Dropper tapered to r=1.5mm + flat 6.35mm Zippy cylinder (misalignment-tolerant!)
#   3. "flat_cyl":     Pure 6.35mm flat cylinder core + flat 6.35mm Zippy cylinder (original baseline shape)
INTERFACE_MODES = ["double_horn", "dropper_horn", "flat_cyl"]

DEFAULT_WORKERS = 16
