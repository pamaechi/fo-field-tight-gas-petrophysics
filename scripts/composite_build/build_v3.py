"""Composite v3 = v2 + BS, CAL_BS, BADHOLE, despiked sonic, shifted core data (CPOR, CKL, CRHOG), F-S1 RHOB_RAW; F-O1 RHOB_ALT removed."""
import sys; sys.path.insert(0, "work")
import pandas as pd, numpy as np, lasio, glob, os
from postproc import postprocess
OUT = "work/out/composites_v3"; os.makedirs(OUT, exist_ok=True)
UNITS = {"BS": "in", "CAL_BS": "in", "BADHOLE": "", "CPOR": "V/V", "CKL": "mD", "CRHOG": "g/cm3", "RHOB_RAW": "g/cm3"}
def core_curves(well, grid):
    c = pd.read_csv("work/core_all_shifted.csv"); c = c[c.well == well]
    out = pd.DataFrame(index=grid)
    if c.empty: return out, 0
    idx = np.clip(np.searchsorted(grid, c.depth_log_m.values), 1, len(grid) - 1)
    near = np.where(np.abs(grid[idx - 1] - c.depth_log_m.values) < np.abs(grid[idx] - c.depth_log_m.values), idx - 1, idx)
    c = c.assign(g=grid[near])
    agg = c.groupby("g").agg(CPOR=("phi_pct", "mean"), CKL=("kkl_md", "mean"), CRHOG=("rhog", "mean"))
    agg["CPOR"] /= 100.0
    return out.join(agg), len(c)
def build(well):
    v2 = pd.read_csv(glob.glob(f"work/out/composites/{well}_composite_v2.csv")[0], index_col=0)
    units = {c.mnemonic: c.unit for c in lasio.read(f"work/out/composites/{well}_composite_v2.las").curves}
    d, notes = postprocess(well, v2)
    if "RHOB_ALT" in d: d = d.drop(columns="RHOB_ALT"); notes["RHOB_ALT"] = "removed: core porosity shows +5 pu bias (RHOB bias -0.1 pu)"
    if well == "F-S1" and "DRHO" in d:
        d["RHOB_RAW"] = d["RHOB"] - d["DRHO"]; notes["RHOB_RAW"] = "RHOB - DRHO; matches core porosity (bias -0.4 pu) whereas RHOB reads 6 pu low (DRHO over-correction, barite mud)"
    cc, n = core_curves(well, d.index.values); d = d.join(cc)
    if n: notes["CORE"] = f"{n} core plugs placed at log depth (core depth + per-core shift, see core_shift_table.csv)"
    d.to_csv(f"{OUT}/{well}_composite_v3.csv", float_format="%.5f")
    las = lasio.LASFile(); las.well["WELL"].value = well; las.well["FLD"].value = "F-O Gas Field"; las.well["CTRY"].value = "South Africa"
    las.other = (f"{well} composite v3. Logger depth, m below KB, 0.1524 m step. Built from v2 (spliced + rule-cleaned). " +
                 " | ".join(f"{k}: {v}" for k, v in notes.items()))
    las.append_curve("DEPT", d.index.values, unit="M")
    for c in d.columns: las.append_curve(c, d[c].values, unit=UNITS.get(c, units.get(c, "")))
    las.write(f"{OUT}/{well}_composite_v3.las", version=2.0, fmt="%.5f")
    return notes
if __name__ == "__main__":
    wells = sys.argv[1:] or ["F-O1", "F-O2", "F-O3", "F-O4", "F-O6", "F-O8", "F-R1", "F-S1"]
    rows = []
    for w in wells:
        n = build(w); rows += [dict(well=w, item=k, note=v) for k, v in n.items()]; print(w, "ok")
    pd.DataFrame(rows).to_csv(f"{OUT}/v3_processing_notes.csv", index=False)
