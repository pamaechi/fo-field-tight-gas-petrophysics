"""Resistivity exactly at tool limit (1950) is KEPT in open hole (means >=1950); only held runs, in-casing data, >20000 spikes and SFLU 9696 fill are removed.
Rule-based log cleaning applied after splicing (v2). Every nulled interval is logged to <well>_cleaning_log.csv."""
import re, numpy as np, pandas as pd
RES = re.compile(r"^(LLD|LLS|ILD|ILM|SFLU|SFL|MSFL|RXO|AT\d0|AF\d0|AO\d0|RLA\d|HDRS|HMRS)$")
CAL = {"CALI", "CALS", "LCAL", "HCAL", "C1", "C2"}
ZERO_OK = {"KTIM", "KSDR", "CMFF", "BFV"}
CASED = {  # (top, base, curves kept as-is besides GR->GR_CH)
    "F-O4": (160.0, 3117.0, set()),
    "F-O6": (2700.0, 3095.0, {"DTSM", "DTCO", "DT4P"}),
}
def _runs(s, minlen):
    s = s.dropna(); v = s.values; idx = s.index.values; out = []; i0 = 0
    for i in range(1, len(v) + 1):
        if i == len(v) or abs(v[i] - v[i0]) > 1e-4 * max(1.0, abs(v[i0])):
            if idx[i - 1] - idx[i0] >= minlen: out.append((idx[i0], idx[i - 1], v[i0]))
            i0 = i
    return out
def clean(well, comp, minlen=3.0):
    log = []; comp = comp.copy()
    def kill(c, mask, why):
        n = int((mask & comp[c].notna()).sum())
        if n:
            dd = comp.index[mask & comp[c].notna()]
            log.append(dict(well=well, curve=c, top_m=round(dd.min(), 2), base_m=round(dd.max(), 2), samples=n, reason=why))
            comp.loc[mask, c] = np.nan
    if well in CASED:
        t, b, keep = CASED[well]; m = (comp.index >= t) & (comp.index <= b)
        if "GR" in comp:
            comp["GR_CH"] = comp.get("GR_CH", pd.Series(np.nan, index=comp.index))
            fill = m & comp["GR_CH"].isna(); comp.loc[fill, "GR_CH"] = comp.loc[fill, "GR"]
        for c in comp.columns:
            if c not in keep and c != "GR_CH": kill(c, m, "logged inside 9-5/8 in casing")
    for c in list(comp.columns):
        s = comp[c]
        if RES.match(c):
            kill(c, s > 20000, "resistivity spike > 20000 ohm.m")
            if c in ("SFLU", "SFL"): kill(c, s >= 9000, "SFLU fill value (9696)")
        if c == "PEF": kill(c, s > 15, "PEF > 15 (junk)")
        if c in CAL: kill(c, s < 6, "caliper < 6 in (arms closed / in casing)")
        if c in ("RHOB", "RHOB_ALT"): kill(c, s < 1.7, "RHOB < 1.7 (tool in casing/junk)")
        for t0, b0, v in _runs(comp[c], minlen):
            if c in ZERO_OK and abs(v) < 1e-3: continue          # NMR perm / FFV floor in tight rock is real
            if c in CAL and v >= 6: continue
            kill(c, (comp.index >= t0) & (comp.index <= b0), f"held constant value {v:.4g} over {b0 - t0:.1f} m")
    comp = comp.dropna(axis=1, how="all").dropna(how="all")
    return comp, pd.DataFrame(log)

PASS_LOG = []
def clean_series(name, s, minlen=3.0, src=""):
    """Same rules as clean(), applied to one source pass before splicing (casing rule is applied later, per well)."""
    s = s.copy()
    def k(mask, why):
        mask = mask & s.notna()
        if mask.any():
            dd = s.index[mask]; PASS_LOG.append(dict(curve=name, source=src, top_m=round(float(dd.min()), 2), base_m=round(float(dd.max()), 2), samples=int(mask.sum()), reason=why))
            s[mask] = np.nan
    if RES.match(name):
        k(s > 20000, "resistivity spike > 20000 ohm.m")
        if name in ("SFLU", "SFL"): k(s >= 9000, "SFLU fill value (9696)")
    if name == "PEF": k(s > 15, "PEF > 15 (junk)")
    if name in CAL: k(s < 6, "caliper < 6 in (arms closed / in casing)")
    if name in ("RHOB", "RHOB_ALT"): k(s < 1.7, "RHOB < 1.7 (junk)")
    ss = s[~s.index.duplicated()].sort_index()
    for t0, b0, v in _runs(ss, minlen):
        if name in ZERO_OK and abs(v) < 1e-3: continue
        if name in CAL and v >= 6: continue
        k((s.index >= t0) & (s.index <= b0), f"held constant value {v:.4g} over {b0 - t0:.1f} m")
    return s
