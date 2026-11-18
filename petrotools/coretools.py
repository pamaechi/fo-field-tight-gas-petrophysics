"""Core-log integration tools for the F-O field (Notebook 02).

Conventions
-----------
* Depths in metres below KB. log depth = core depth + shift.
* Core porosity in %, permeability in mD, grain density in g/cc.
* Composites are the v2 CSV files (index = DEPT) from F-O_composites_v2.zip.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------------------------
# 1. Core-to-log depth matching
# ----------------------------------------------------------------------------------------------
def shift_scan(core: pd.DataFrame, log: pd.DataFrame, curve: str = "RHOB", smin: float = -8, smax: float = 8,
               step: float = 0.0762, sign: int = -1, minn: int = 6) -> pd.DataFrame:
    """Correlation between core porosity and a log curve for every trial shift.
    sign=-1 for density/sonic-slowness behaving as -porosity (RHOB), +1 for NPHI / DT (porosity-like)."""
    z = core.depth_m.to_numpy(); p = core.phi_pct.to_numpy(); ok = np.isfinite(p); z, p = z[ok], p[ok]
    L = log[curve].dropna()
    out = []
    for s in np.arange(smin, smax + 1e-9, step):
        v = np.interp(z + s, L.index.to_numpy(), L.to_numpy(), left=np.nan, right=np.nan); m = np.isfinite(v)
        r = np.corrcoef(p[m], sign * v[m])[0, 1] if m.sum() >= minn else np.nan
        out.append((round(s, 4), r))
    return pd.DataFrame(out, columns=["shift_m", "r"])

def best_shift(scan: pd.DataFrame, sep: float = 1.5) -> dict:
    """Best shift, its r, r at zero shift, and the best r more than `sep` m away (ambiguity check)."""
    s = scan.dropna()
    if s.empty: return dict(best_shift=np.nan, r_best=np.nan, r_zero=np.nan, second_r=np.nan)
    b = s.loc[s.r.idxmax()]; far = s[(s.shift_m - b.shift_m).abs() > sep]
    return dict(best_shift=float(b.shift_m), r_best=float(b.r), r_zero=float(s.loc[s.shift_m.abs().idxmin(), "r"]),
                second_r=float(far.r.max()) if len(far) else np.nan)

def choose_shifts(results: pd.DataFrame, overrides: dict | None = None, r_min: float = 0.75, margin: float = 0.10,
                  n_min: int = 6) -> pd.DataFrame:
    """Per-core decision. results needs: well, core, best_shift, r_best, second_r, n.
    Robust = r_best >= r_min, clear of the 2nd peak by `margin`, n >= n_min. Non-robust cores take the well median
    of robust cores; a well with no robust core keeps its own best shift if r >= 0.6 and margin >= 0.15 (MODERATE).
    `overrides` = {(well, core): (shift, reason)} for multi-log consensus decisions."""
    overrides = overrides or {}
    R = results.copy()
    R["robust"] = (R.r_best >= r_min) & ((R.r_best - R.second_r).fillna(1) >= margin) & (R.n >= n_min)
    med = R[R.robust].groupby("well").best_shift.median()
    rows = []
    for _, r in R.iterrows():
        key = (r.well, int(r.core))
        if key in overrides: sh, how = overrides[key]
        elif r.robust: sh, how = r.best_shift, f"core-log correlation (r={r.r_best:.2f}, n={int(r.n)})"
        elif r.well in med.index: sh, how = med[r.well], f"well median of robust cores (unconstrained: r={r.r_best:.2f}, 2nd {r.second_r:.2f}, n={int(r.n)})"
        elif r.r_best >= 0.6 and (r.r_best - r.second_r) >= 0.15: sh, how = r.best_shift, f"core-log correlation, MODERATE (r={r.r_best:.2f}, 2nd {r.second_r:.2f})"
        else: sh, how = np.nan, "no reliable match"
        rows.append(dict(well=r.well, core=int(r.core), shift_m=round(float(sh), 2), method=how))
    return pd.DataFrame(rows)

def apply_shifts(core: pd.DataFrame, shifts: pd.DataFrame) -> pd.DataFrame:
    c = core.merge(shifts[["well", "core", "shift_m"]], on=["well", "core"], how="left")
    c["depth_log_m"] = c.depth_m + c.shift_m
    return c

# ----------------------------------------------------------------------------------------------
# 2. Density calibration check against core
# ----------------------------------------------------------------------------------------------
def density_bias(core: pd.DataFrame, log: pd.DataFrame, curve: str = "RHOB", rhof: float = 1.0, rhoma: float | None = None) -> dict:
    """Mean (phiD - phi_core) in porosity units at shifted core depths. rhoma=None uses plug grain density."""
    g = core[core.phi_pct.notna() & core.depth_log_m.notna()]
    s = log[curve] if isinstance(curve, str) else curve
    L = s.dropna(); rb = np.interp(g.depth_log_m, L.index, L.to_numpy(), left=np.nan, right=np.nan)
    rg = np.full(len(g), rhoma) if rhoma else g.rhog.fillna(g.rhog.median() if g.rhog.notna().any() else 2.67).to_numpy()
    phid = 100 * (rg - rb) / (rg - rhof); m = np.isfinite(phid); d = phid[m] - g.phi_pct.to_numpy()[m]
    return dict(n=int(m.sum()), bias_mean_pu=round(float(d.mean()), 2), bias_median_pu=round(float(np.median(d)), 2),
                rmse_pu=round(float(np.sqrt(np.mean(d ** 2))), 2), r=round(float(np.corrcoef(phid[m], g.phi_pct.to_numpy()[m])[0, 1]), 2))

# ----------------------------------------------------------------------------------------------
# 3. Bit size, bad-hole flag, sonic despike
# ----------------------------------------------------------------------------------------------
SHOE = {  # 9-5/8 in casing shoe (m MD) and evidence
    "F-O1": (3089.0, "run-3 top; caliper 12.5 -> 9-13 in"),
    "F-O2": (3446.0, "casing signature above: LLS 0.15 ohm.m, DT 59.6 us/ft"),
    "F-O3": (3398.0, "well completion report Table 1"),
    "F-O4": (3117.0, "drilling completion report"),
    "F-S1": (3500.0, "casing signature above: DT ~56 us/ft, CALI 8.6 in, DRHO -0.86"),
}
def caliper(d: pd.DataFrame) -> str | None:
    for c in ("CALI", "LCAL", "CALS"):
        if c in d and d[c].notna().sum() > 100: return c

def _break_depth(cal: pd.Series, thr: float, win_m: float = 20.0, step: float = 0.1524):
    n = max(3, int(win_m / step)); p10 = cal.rolling(n, center=True, min_periods=n // 2).quantile(0.1)
    below = (p10 < thr).astype(float).where(p10.notna()).dropna()
    if below.empty: return None
    rev = below[::-1].expanding().mean()[::-1]; ok = rev[rev >= 0.9]
    return float(ok.index.min()) if len(ok) else None

def bit_size(well: str, d: pd.DataFrame):
    if "LCAL" in d and "DCAL" in d:
        bs = (d.LCAL - d.DCAL).round(2); return bs.where(bs.between(5, 18)).ffill().bfill(), "LCAL - DCAL"
    c = caliper(d); bs = pd.Series(12.25, index=d.index)
    z, src = SHOE.get(well, (_break_depth(d[c].loc[3600:4100], 12.0), "caliper break (ASSUMED)"))
    bs[d.index > z] = 8.5; note = f"12.25 in above / 8.5 in below {z:.1f} m ({src})"
    if well == "F-R1":
        z2 = _break_depth(d[c].loc[:3000], 15.5)
        if z2: bs[d.index < z2] = 17.5; note += f"; 17.5 in above {z2:.1f} m (assumed)"
    return bs, note

def despike(s: pd.Series, win: int = 11, thr: float = 15.0):
    med = s.rolling(win, center=True, min_periods=5).median(); bad = (s - med).abs() > thr
    return s.mask(bad), int(bad.sum())

def postprocess(well: str, d: pd.DataFrame, cal_excess: float = 3.0, drho_max: float = 0.10):
    d = d.copy(); notes = {}
    d["BS"], notes["BS"] = bit_size(well, d)
    c = caliper(d); d["CAL_BS"] = d[c] - d.BS; notes["CAL_BS"] = f"{c} - BS (in)"
    bad = d.CAL_BS > cal_excess
    if "DRHO" in d: bad |= d.DRHO.abs() > drho_max
    d["BADHOLE"] = np.where(d[c].isna(), np.nan, bad.astype(float))
    notes["BADHOLE"] = f"1 where {c} - BS > {cal_excess} in" + (f" or |DRHO| > {drho_max} g/cc" if "DRHO" in d else "")
    for s in ("DT", "DTL", "DTCO", "DT4P", "DTSM", "DT4S"):
        if s in d: d[s], n = despike(d[s], thr=30.0 if s in ("DTSM", "DT4S") else 15.0); notes[s] = f"despiked {n} samples"
    return d, notes

def core_curves(core: pd.DataFrame, grid: np.ndarray) -> pd.DataFrame:
    """Place shifted plugs on the log grid as sparse curves CPOR (v/v), CKL (mD), CRHOG (g/cc)."""
    c = core[core.depth_log_m.notna()]
    if c.empty: return pd.DataFrame(index=grid)
    idx = np.clip(np.searchsorted(grid, c.depth_log_m.to_numpy()), 1, len(grid) - 1)
    near = np.where(np.abs(grid[idx - 1] - c.depth_log_m.to_numpy()) < np.abs(grid[idx] - c.depth_log_m.to_numpy()), idx - 1, idx)
    agg = c.assign(g=grid[near]).groupby("g").agg(CPOR=("phi_pct", "mean"), CKL=("kkl_md", "mean"), CRHOG=("rhog", "mean"))
    agg["CPOR"] /= 100.0
    return pd.DataFrame(index=grid).join(agg)
