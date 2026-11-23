from materials import LATHE_FINISH_GAP_MM


# --- Fixed Hardware & Target Specs (Based on Exp 06 Winner ID=94) ---

R_STOCK_MM = 6.35 / 2.0               # 1/4" (6.35 mm) McMaster Magnetic Iron (mu_r = 4000)
PM_HEIGHT_MM = 2.54                   # 1/10" (2.54 mm) McMaster Grade N52 on Zippy (0.603 g)
CONTACT_GAP_MM = LATHE_FINISH_GAP_MM  # 15 um bare lathe-cut contact

R_SHANK_MM = 2.20                     # Ø4.40 mm turned iron shank
R_TIP_MM = 1.75                       # Ø3.50 mm concentrated horn tip
TAPER_LEN_MM = 1.50                   # Bottom conical horn height
ZIPPY_CAP_H_MM = 1.50                 # Ø6.35 x 1.50 mm iron cap on top of Zippy's N52 (~0.98 g total target)

BOBBIN_RADIAL_WALL_MM = 0.35          # Inner bobbin tube wall thickness
TOP_GLUE_GAP_MM = 0.05                # Coil goes flush to top (glued straight to flat PCB underside!)
WIRE_AWG = 34                         # 34 AWG magnet wire

CAP_C1_F = 1000e-6
PULSE_MS = 22.0


# --- Sweep Parameters: Scramble Hand-Winding Density & Bottom Recession ---

# Bottom recession z_coil_bot = proud iron tip clearance + bottom bobbin cheek + lead exit slot
Z_COIL_BOT_LIST_MM = [0.50, 0.90, 1.30, 1.70, 2.20]

# Realistic hand scramble-winding bare-copper packing fraction:
#   0.35 = very loose/messy blitz hand wind
#   0.42 = typical fast hand wind with light finger tension
#   0.50 = semi-neat hand wind
#   (Note: insulated wire diameter is ~1.09x bare, so 0.42 bare Cu fill = ~63% gross wire cylinder packing!)
SCRAMBLE_FILL_LIST = [0.35, 0.42, 0.50]

# Radial winding window width & core length options to compensate for looser fill + bottom recession
COIL_W_LIST_MM = [1.80, 2.30]         # Outer coil OD = 8.70 mm or 9.70 mm (both tiny vs 18x22mm PCB!)
CORE_LEN_LIST_MM = [7.50, 8.50]

# Phase 2 Full Moving-Gap Transient Verification Matrix:
# Test how z_coil_bot & worst-case 0.35 scramble fill behave across Trimpot V0 = [16, 18, 20, 22, 25] V!
TRANSIENT_VOLTAGES_V = [16.0, 18.0, 20.0, 22.0, 25.0]

DEFAULT_WORKERS = 16
