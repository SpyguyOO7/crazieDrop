import numpy as np


# --- Fixed Stock & Existing 32 AWG Winding Parameters ---

R_STOCK_MM = 6.35 / 2.0        # 1/4" (6.35 mm) McMaster Magnetic Iron (mu_r = 4000)
PM_HEIGHT_MM = 2.54            # 1/10" (2.54 mm) McMaster disk magnet height

CORE_LEN_MM = 7.0              # 7.0 mm dropper iron core
TAPER_LEN_MM = 1.5             # Bottom conical horn taper height
R_TIP_MM = 1.50                # 3.0 mm diameter concentrated pole tip

BOBBIN_WALL_MM = 0.35
COIL_W_MM = 1.80               # 1.80 mm radial build
WIRE_AWG = 32                  # Existing 32 AWG magnet wire!

# Current scan up to deep reverse saturation (0 .. 8.0 A -> 0 .. ~1,570 A-turns)
CURRENTS_A = np.round(np.linspace(0.0, 8.0, 33), 3)

# Trimpot voltage sweep to test voltage-window margin width with 1000 uF capacitor + 32 AWG
TRIMPOT_VOLTAGES_V = [10.0, 13.0, 16.0, 19.0, 22.0, 25.0]
CAP_C1_F = 1000e-6
PULSE_MS = 20.0

DEFAULT_WORKERS = 16
