"""Small reader layer: any LIS/TIF-LIS/DLIS/LAS frame -> clean pandas DataFrame indexed by depth in metres."""
import os, numpy as np, pandas as pd
from dlisio import lis, dlis
from dlisio.common import ErrorHandler, Actions
import lasio
from .tif_repair import tif_repair

H = ErrorHandler(critical=Actions.LOG_ERROR)
U2M = {"M": 1, "METR": 1, "F": .3048, "FT": .3048, ".1IN": .00254, "0.1IN": .00254, "IN": .0254, ".5MM": .0005, "MM": .001, "CM": .01}
REP = os.path.join(os.path.expanduser("~"), ".fo_repaired_lis")

def _c(s): return (s.decode(errors="ignore") if isinstance(s, bytes) else str(s or "")).strip()

def lis_frames(path, well="x"):
    """Yield (key, DataFrame, units dict) for every DFSR in a LIS/TIF-LIS file (sample-rate-1 curves)."""
    src = path
    if open(path, "rb").read(4) == b"\0\0\0\0":
        os.makedirs(REP, exist_ok=True)
        src = f"{REP}/{well}__{os.path.basename(path)}.lis"
        if not os.path.exists(src): tif_repair(path, src)
    with lis.load(src, error_handler=H) as lf:
        for li, f in enumerate(lf):
            for di, fs in enumerate(f.data_format_specs()):
                rates = sorted({s.samples for i, s in enumerate(fs.specs) if not (fs.depth_mode == 0 and i == 0)})
                rate = 1 if 1 in rates else (rates[0] if rates else 1)
                try: a = lis.curves(f, fs, sample_rate=rate, strict=False)
                except Exception: continue
                if a.size == 0: continue
                idx = a.dtype.names[0]
                u = _c(fs.specs[0].units) if fs.depth_mode == 0 else _c(fs.index_units or fs.depth_units)
                fac = U2M.get(u.upper(), np.nan)
                cols = {}
                units = {}
                for s in fs.specs:
                    m = _c(s.mnemonic)
                    raw = s.mnemonic if isinstance(s.mnemonic, str) else _c(s.mnemonic)
                    for cand in (raw, m):
                        if cand in a.dtype.names and a.dtype[cand].shape == () and a.dtype[cand].kind in "iuf":
                            v = a[cand].astype(float)
                            v[(v == fs.absent_value) | (v <= -999) ] = np.nan
                            cols[m] = v; units[m] = _c(s.units); break
                df = pd.DataFrame(cols)
                df.index = a[idx].astype(float) * fac
                df.index.name = "DEPTH_M"
                df = df.loc[:, [c for c in df.columns if c.upper() not in ("DEPT", "DEPTH")]]
                yield f"L{li}D{di}", df.sort_index(), units

def dlis_frames(path, max_mb=2000):
    with dlis.load(path, error_handler=H) as lf:
        for li, f in enumerate(lf):
            for fr in f.frames:
                chans = [ch for ch in fr.channels if ch.dimension in ([1], [])]
                if len(chans) < 2: continue
                try: a = fr.curves()
                except Exception: continue
                idx = fr.channels[0].name
                u = _c(fr.channels[0].units)
                fac = U2M.get(u.upper(), {"0.1 in": .00254}.get(u, np.nan))
                cols = {ch.name: a[ch.name].astype(float) for ch in chans[1:] if ch.name in a.dtype.names and a.dtype[ch.name].kind in "iuf"}
                df = pd.DataFrame(cols); df.index = a[idx].astype(float) * fac; df.index.name = "DEPTH_M"
                df = df.replace([-999.25, -9999, -999], np.nan)
                yield f"L{li}:{fr.name}", df.sort_index(), {ch.name: _c(ch.units) for ch in chans}

def las_frame(path):
    l = lasio.read(path, ignore_header_errors=True)
    df = l.df(); u = _c(l.curves[0].unit).upper() or "M"
    df.index = df.index.astype(float) * U2M.get(u, 1); df.index.name = "DEPTH_M"
    return df.sort_index(), {c.mnemonic: c.unit for c in l.curves}
