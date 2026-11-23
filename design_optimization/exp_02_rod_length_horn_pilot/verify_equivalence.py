from pathlib import Path
import csv

from run_sweep import eval_single_task

# Verify that the fast StaticSolverSession produces the EXACT same physics numbers
# as the already-saved slow run on row 1 of results.csv (task_id=16, double_horn, L=5.5, w=0.8)

csv_path = Path(__file__).resolve().parent / "results.csv"
with csv_path.open("r", newline="") as f:
    ref_row = next(csv.DictReader(f))

new_row = eval_single_task((
    int(ref_row["task_id"]),
    ref_row["interface_mode"],
    float(ref_row["core_len_mm"]),
    float(ref_row["coil_w_mm"]),
))

print("=== Bit-for-Bit Physics Equivalence Check vs Saved Slow Solver ===")
print(f"Old slow eval time: {ref_row['eval_time_s']} s  -->  New fast eval time: {new_row['eval_time_s']} s")
for key in ["f_hold_n", "f_hold_gf", "b_gap_hold_t", "ni_null_at", "p_dc_null_w", "pm_demag_margin_null_ka_m"]:
    old_v = float(ref_row[key])
    new_v = float(new_row[key])
    print(f"  {key:26s}: old = {old_v:10.5f} | new = {new_v:10.5f} | diff = {abs(new_v - old_v):.6f}")
