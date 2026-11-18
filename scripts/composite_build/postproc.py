"""Composite v3 post-processing: bit size (BS), bad-hole flag (BADHOLE), sonic despike.
BS: F-O6/F-O8 exact (LCAL - DCAL); other wells 12.25 in above / 8.5 in below the depth where the caliper 10th percentile
(20 m window) drops below 10.4 in and stays there; F-R1 also 17.5 in above the first break below 15.5 in."""
import numpy as np, pandas as pd
CALPREF = ("CALI", "LCAL", "CALS")
def caliper(d):
    for c in CALPREF:
        if c in d and d[c].notna().sum() > 100: return c
def break_depth(cal, thr, win_m=20.0, step=0.1524):
    n = max(3, int(win_m / step)); p10 = cal.rolling(n, center=True, min_periods=n // 2).quantile(0.1)
    below = (p10 < thr).astype(float).where(p10.notna())
    # first depth after which >= 90 % of valid samples stay below thr
    idx = below.dropna()
    if idx.empty: return None
    rev = idx[::-1].expanding().mean()[::-1]
    ok = rev[rev >= 0.9]
    return float(ok.index.min()) if len(ok) else None
SHOE = {  # 9-5/8 in casing shoe (m MD): documented, or from casing-ID caliper (8.6 in) / steel sonic at run tops
    "F-O1": (3089.0, "run-3 top; caliper 12.5 -> 9-13 in"), "F-O2": (3446.0, "casing signature above: LLS 0.15, DT 59.6"),
    "F-O3": (3398.0, "well completion report Table 1"), "F-O4": (3117.0, "drilling completion report"),
    "F-S1": (3500.0, "casing signature above: DT ~56, CALI 8.6, DRHO -0.86")}
def bit_size(well, d):
    if "LCAL" in d and "DCAL" in d:
        bs = (d["LCAL"] - d["DCAL"]).round(2)
        return bs.where(bs.between(5, 18)).ffill().bfill(), "LCAL - DCAL"
    c = caliper(d); bs = pd.Series(12.25, index=d.index)
    if well in SHOE:
        z85, src = SHOE[well]
    else:                                   # F-R1: not documented in the text reports -> caliper break (assumed)
        z85, src = break_depth(d[c].loc[3600:4100], 12.0), "caliper break (ASSUMED, not documented)"
    bs[d.index > z85] = 8.5
    note = f"12.25 in above / 8.5 in below 9-5/8 in shoe at {z85:.1f} m ({src})"
    if well == "F-R1":
        z175 = break_depth(d[c].loc[:3000], 15.5)
        if z175: bs[d.index < z175] = 17.5; note += f"; 17.5 in above {z175:.1f} m (caliper break, assumed)"
    return bs, note
def despike(s, win=11, thr=15.0):
    med = s.rolling(win, center=True, min_periods=5).median(); bad = (s - med).abs() > thr
    return s.mask(bad), int(bad.sum())
def postprocess(well, d):
    d = d.copy(); notes = {}
    d["BS"], notes["BS"] = bit_size(well, d)
    c = caliper(d)
    d["CAL_BS"] = d[c] - d["BS"]; notes["CAL_BS"] = f"{c} - BS (in)"
    bad = d["CAL_BS"] > 3.0
    if "DRHO" in d: bad |= d["DRHO"].abs() > 0.10
    d["BADHOLE"] = np.where(d[c].isna(), np.nan, bad.astype(float))
    notes["BADHOLE"] = f"1 where {c} - BS > 3.0 in" + (" or |DRHO| > 0.10 g/cc" if "DRHO" in d else "")
    for s in ("DT", "DTL", "DTCO", "DT4P", "DTSM", "DT4S"):
        if s in d:
            d[s], n = despike(d[s], thr=30.0 if s in ("DTSM", "DT4S") else 15.0); notes[s] = f"despiked {n} samples"
    return d, notes
