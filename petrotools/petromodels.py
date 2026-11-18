"""Shale volume, core-calibrated porosity and permeability models for the F-O field (Notebook 03)."""
from __future__ import annotations
import numpy as np
import pandas as pd

# reservoir windows: RPS (2008) M10 (BLVS) to M1; F-O6 / F-O8 have no RPS picks -> 3700 m to TD
WINDOWS = {"F-O1": (3695.0, 3815.5), "F-O2": (3595.0, 3699.5), "F-O3": (3678.0, 3794.0), "F-O4": (3700.0, 3876.0),
           "F-O6": (3700.0, 3930.0), "F-O8": (3700.0, 3916.0), "F-R1": (3836.0, 4035.0), "F-S1": (3596.0, 3712.5)}

def harmonise(well: str, d: pd.DataFrame) -> pd.DataFrame:
    """Unit fixes and the density curve to use. F-O8 POTA is stored as a fraction -> %; F-O6 spectral GR renamed to
    POTA/THOR/URAN (HFK fraction -> %); RHOB_USE = RHOB_RAW for F-S1 (DRHO over-correction), else RHOB."""
    d = d.copy()
    if well == "F-O8" and "POTA" in d and d.POTA.median() < 0.2: d["POTA"] = d.POTA * 100
    if well == "F-O6" and "HFK" in d:
        d["POTA"] = d.HFK * 100; d["THOR"] = d.HTHO; d["URAN"] = d.HURA
    for c, (lo, hi) in {"POTA": (0, 10), "THOR": (0, 40), "URAN": (0, 30)}.items():   # spectral-GR spikes (NGT/HNGS)
        if c in d: d[c] = despike(d[c].where(d[c].between(lo, hi)), thr=3 * float(d[c].std() or 0) or np.inf)
    for c in ("KTIM", "KSDR"):                                                            # NMR permeability spikes, in log space
        if c in d: d[c] = despike_log(d[c])
    d["RHOB_USE"] = d.RHOB_RAW if "RHOB_RAW" in d else d.RHOB
    if well == "F-S1":   # bad hole from caliper only: RHOB_RAW already handles the DRHO problem
        d["BADHOLE"] = np.where(d.CAL_BS.isna(), np.nan, (d.CAL_BS > 3.0).astype(float))
    return d

def despike(s: pd.Series, win: int = 11, thr: float = np.inf) -> pd.Series:
    med = s.rolling(win, center=True, min_periods=5).median(); return s.mask((s - med).abs() > thr)

def despike_log(s: pd.Series, win: int = 11, thr_dex: float = 1.5, kmax: float = 2000.0) -> pd.Series:
    """Remove samples > thr_dex decades from the running median of log10(k), and anything above kmax mD."""
    l = np.log10(s.clip(lower=1e-5)); med = l.rolling(win, center=True, min_periods=5).median()
    return s.mask(((l - med).abs() > thr_dex) | (s > kmax))

def gr_endpoints(d: pd.DataFrame, window: tuple, curve: str = "GR", pad: float = 60.0, lo: float = 0.05, hi: float = 0.95):
    s = d[curve].loc[window[0] - pad: window[1] + pad].dropna()
    return float(s.quantile(lo)), float(s.quantile(hi))

def vsh_linear(x: pd.Series, clean: float, shale: float) -> pd.Series:
    return ((x - clean) / (shale - clean)).clip(0, 1)

def vsh_larionov_old(igr: pd.Series) -> pd.Series:
    return (0.33 * (2 ** (2 * igr) - 1)).clip(0, 1)

def phit_density(rhob: pd.Series, rho_ma: float, rho_f: float = 1.0, calib_pu: float = 0.0) -> pd.Series:
    """Total porosity (v/v) from density; calib_pu = per-well additive calibration to core (porosity units)."""
    return ((rho_ma - rhob) / (rho_ma - rho_f) + calib_pu / 100.0).clip(0, 0.35)

def fit_logk(phi_pct: pd.Series, k_md: pd.Series, kmin: float = 0.001):
    m = phi_pct.notna() & (k_md > kmin)
    b = np.polyfit(phi_pct[m], np.log10(k_md[m]), 1); res = np.log10(k_md[m]) - np.polyval(b, phi_pct[m])
    return dict(a=float(b[1]), b=float(b[0]), sd_dex=float(res.std()), r=float(np.corrcoef(phi_pct[m], np.log10(k_md[m]))[0, 1]), n=int(m.sum()))

def perm_from_phi(phit: pd.Series, a: float, b: float) -> pd.Series:
    return 10 ** (a + b * 100 * phit)

def match_logs(core: pd.DataFrame, logs: dict, curves: list, win: int = 1) -> pd.DataFrame:
    """Log values at each plug's log depth (mean of 2*win+1 samples)."""
    rows = []
    for w, g in core.groupby("well"):
        d = logs[w]; z = d.index.to_numpy()
        for _, r in g.iterrows():
            i = int(np.clip(np.searchsorted(z, r.depth_log_m), win, len(z) - win - 1)); seg = d.iloc[i - win: i + win + 1]
            rows.append({**r.to_dict(), **{c: seg[c].mean() if c in d else np.nan for c in curves}})
    return pd.DataFrame(rows)

def error_stats(pred: pd.Series, obs: pd.Series) -> dict:
    m = pred.notna() & obs.notna(); e = pred[m] - obs[m]
    return dict(n=int(m.sum()), bias=round(float(e.mean()), 2), rmse=round(float(np.sqrt((e ** 2).mean())), 2),
                r=round(float(np.corrcoef(pred[m], obs[m])[0, 1]), 2))
