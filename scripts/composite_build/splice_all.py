"""Splice an IP-ready composite for every well from the raw tapes.
Each segment: (file glob, frame-key prefix, {source: output}, top, base, role)  role: 'ref' = depth reference, else matched to it.
Bulk depth shift applied to a non-reference segment only if GR cross-correlation improves by >=0.03 and |lag| <= 1.0 m."""
import os, sys, json, numpy as np, pandas as pd, lasio, warnings, logging
warnings.filterwarnings("ignore"); logging.getLogger("dlisio").setLevel(logging.CRITICAL)
sys.path.insert(0, "work")
from lagcheck import frame
try:
    from clean import clean_series
except Exception:
    clean_series = None

STEP = 0.1524
GRN = ("GR", "GR_EDTC", "ECGR", "SGR")

def grcol(df):
    for c in GRN:
        if c in df and df[c].notna().sum() > 200: return c

def resample(s, grid, src_step):
    s = s[~s.index.duplicated()].sort_index().astype(float)
    if src_step < STEP * 0.8:                     # finer than grid: smooth to grid resolution first
        w = max(1, int(round(STEP / src_step)))
        s = s.rolling(w, center=True, min_periods=1).mean()
    v = np.interp(grid, s.index.values, s.values, left=np.nan, right=np.nan)
    ok = np.interp(grid, s.index.values, s.notna().values.astype(float), left=0, right=0) > 0.99
    v[~ok] = np.nan
    return v

def best_lag(A, B, maxlag=1.0):
    ca, cb = grcol(A), grcol(B)
    if not ca or not cb: return None
    lo, hi = max(A.index.min(), B.index.min()), min(A.index.max(), B.index.max())
    if hi - lo < 40: return None
    g = np.arange(lo, hi, STEP)
    x = np.interp(g, A.index, A[ca].values); y = np.interp(g, B.index, B[cb].values)
    def r(L):
        xx, yy = (x[L:], y[:len(y) - L]) if L >= 0 else (x[:L], y[-L:])
        m = np.isfinite(xx) & np.isfinite(yy)
        return np.corrcoef(xx[m], yy[m])[0, 1] if m.sum() > 200 else np.nan
    n = int(maxlag / STEP)
    rs = {L: r(L) for L in range(-n, n + 1)}
    L = max(rs, key=lambda k: -1 if np.isnan(rs[k]) else rs[k])
    return dict(lag_m=round(L * STEP, 3), r0=round(rs[0], 3), rbest=round(rs[L], 3))

def splice(well, segs, top, base, notes, units_override=None):
    grid = np.round(np.arange(top, base + STEP / 2, STEP), 4)
    comp = pd.DataFrame(index=pd.Index(grid, name="DEPT"))
    units, prov = {}, []
    ref = None
    for pat, key, cmap, t, b, role in segs:
        df, u = frame(well, pat, key)
        src_step = float(np.median(np.abs(np.diff(df.index.values[np.isfinite(df.index.values)]))))
        shift, lagi = 0.0, None
        if "@shift=" in role:                      # forced bulk shift from windowed GR lag analysis
            shift = float(role.split("@shift=")[1]); lagi = best_lag(ref, df) if ref is not None else None
        elif role == "ref":
            ref = df
        elif ref is not None:
            lagi = best_lag(ref, df)
            if lagi and lagi["rbest"] - lagi["r0"] >= 0.03 and abs(lagi["lag_m"]) > 0:
                shift = lagi["lag_m"]
        d2 = df.copy(); d2.index = d2.index + shift
        d2 = d2.loc[t:b]
        for src, dst in cmap.items():
            if src not in d2 or d2[src].notna().sum() < 10: continue
            ser = d2[src].astype(float).copy()
            if clean_series is not None: ser = clean_series(dst, ser, src=f"{pat} {key}")   # per-pass cleaning so deeper/other passes can fill gaps
            if ser.notna().sum() < 10: continue
            v = resample(ser, grid, src_step)
            v[(grid < t) | (grid > b)] = np.nan
            if dst not in comp: comp[dst] = np.nan
            fill = comp[dst].isna().values & ~np.isnan(v)
            comp.loc[fill, dst] = v[fill]
            units.setdefault(dst, (units_override or {}).get(dst, u.get(src, "")))
            prov.append(dict(well=well, curve=dst, source_file=os.path.basename(__import__('glob').glob(f"data/raw/{well}/Raw_logs/{pat}")[0]),
                             frame=key, source_mnem=src, top_m=t, base_m=b, role=role, depth_shift_m=shift,
                             gr_lag_check=json.dumps(lagi) if lagi else ""))
    comp = comp.dropna(how="all")
    for c in ("NPHI", "TNPH", "NPOR"):          # neutron porosity to v/v
        if c in comp and np.nanmedian(comp[c]) > 1.0:
            comp[c] = comp[c] / 100.0; units[c] = "V/V"
        elif c in comp: units[c] = "V/V"
    return comp, units, pd.DataFrame(prov)

RANGES = {"GR": (5, 250), "RHOB": (1.8, 3.0), "NPHI": (-0.05, 0.6), "TNPH": (-0.05, 0.6), "NPOR": (-0.05, 0.6), "PEF": (1, 10),
          "DT": (40, 160), "DTL": (40, 160), "DTCO": (40, 160), "DT4P": (40, 160), "DTSM": (60, 400), "DT4S": (60, 400),
          "CALI": (5, 20), "CALS": (5, 20), "DCAL": (-2, 8), "LCAL": (5, 20), "DRHO": (-0.2, 0.3)}
def sanity(comp):
    out = []
    for c, (lo, hi) in RANGES.items():
        if c in comp:
            v = comp[c].dropna()
            if len(v): out.append((c, round(float(v.median()), 3), round(float(((v < lo) | (v > hi)).mean() * 100), 1)))
    return out

def write(well, comp, units, prov, notes, outdir):
    os.makedirs(outdir, exist_ok=True)
    comp.to_csv(f"{outdir}/{well}_composite_v2.csv", float_format="%.4f")
    prov.to_csv(f"{outdir}/{well}_provenance.csv", index=False)
    las = lasio.LASFile()
    las.well["WELL"].value = well; las.well["FLD"].value = "F-O Gas Field"; las.well["CTRY"].value = "South Africa"
    las.well["LOC"].value = "Block 9, Bredasdorp Basin, offshore South Africa"
    las.other = (f"{well} composite v2 spliced from raw contractor tapes (LIS/TIF-LIS/DLIS/LAS) with dlisio/lasio. "
                 "Depth = logger measured depth, m below KB, 0.1524 m step. No core depth shift applied. Rule-based cleaning applied (held values, saturated/fill values, in-casing data) - see cleaning log CSV. "
                 "Curve sources, bulk depth shifts and GR lag checks are in the provenance CSV. " + notes)
    las.append_curve("DEPT", comp.index.values, unit="M", descr="Measured depth")
    for c in comp.columns:
        las.append_curve(c, comp[c].values, unit=units.get(c, ""))
    las.write(f"{outdir}/{well}_composite_v2.las", version=2.0, fmt="%.4f")
