"""Inventory every raw log file in the data room (LIS / TIF-LIS / DLIS / LAS).
Writes work/inventory/all_frames.csv (one row per frame/DFSR) and all_curves.csv."""
import os, glob, sys, logging, warnings
import numpy as np, pandas as pd
from dlisio import lis, dlis
from dlisio.common import ErrorHandler, Actions
import lasio

logging.getLogger("dlisio").setLevel(logging.CRITICAL)
warnings.filterwarnings("ignore")
H = ErrorHandler(critical=Actions.LOG_ERROR)
ROOT = "data/raw"
OUT = "work/inventory"; os.makedirs(OUT, exist_ok=True)
U2M = {"M": 1, "METR": 1, "F": .3048, "FT": .3048, "FEET": .3048, ".1IN": .00254, "0.1IN": .00254, "0.1 IN": .00254, ".5MM": .0005, "MM": .001, "CM": .01, "IN": .0254}

def c(s): return (s.decode(errors="ignore") if isinstance(s, bytes) else str(s or "")).strip()

def sniff(p):
    b = open(p, "rb").read(120)
    if not b: return "EMPTY"
    if b"RECORD" in b[:90]: return "DLIS"
    if b.lstrip()[:2] in (b"~V", b"~v"): return "LAS"
    if b[:4] == b"\0\0\0\0" or b[0] in (0x84, 0x82, 0x80): return "LIS"
    if b[:4] == b"\x01\x19\xf1\xf8": return "PDS"
    if b.startswith(b"%PDF"): return "PDF"
    if b"PrimaryIndex" in b: return "ASC"
    return "OTHER"

frames, curves = [], []

def add_frame(well, path, fmt, key, idx, units, top, base, step, n, names, note=""):
    f = top * U2M.get(units.upper(), np.nan) if top is not None else np.nan
    g = base * U2M.get(units.upper(), np.nan) if base is not None else np.nan
    frames.append(dict(well=well, file=os.path.basename(path), fmt=fmt, frame=key, index=idx, index_units=units,
                       top_m=min(f, g) if np.isfinite(f) else np.nan, base_m=max(f, g) if np.isfinite(g) else np.nan,
                       step=step, n=n, n_curves=len(names), curves=" ".join(names), note=note))

from tif_repair import tif_repair
REP = "work/repaired"; os.makedirs(REP, exist_ok=True)
def do_lis(well, p):
    src = p
    if open(p, "rb").read(4) == b"\0\0\0\0":
        src = f"{REP}/{well}__{os.path.basename(p)}.lis"
        if not os.path.exists(src): tif_repair(p, src)
    with lis.load(src, error_handler=H) as lf:
        for li, f in enumerate(lf):
            for di, fs in enumerate(f.data_format_specs()):
                key = f"L{li}D{di}"
                rates = sorted({s.samples for i, s in enumerate(fs.specs) if not (fs.depth_mode == 0 and i == 0)})
                rate = 1 if 1 in rates else (rates[0] if rates else 1)
                try:
                    a = lis.curves(f, fs, sample_rate=rate, strict=False)
                except Exception as e:
                    add_frame(well, p, "LIS", key, "", "", None, None, None, 0, [], f"ERR {e}"); continue
                if a.size == 0: continue
                idx = a.dtype.names[0]
                u = c(fs.specs[0].units) if fs.depth_mode == 0 else c(fs.index_units or fs.depth_units)
                d = a[idx].astype(float)
                names = [n for n in a.dtype.names[1:] if a.dtype[n].shape == ()]
                add_frame(well, p, "LIS", key, idx, u, np.nanmin(d), np.nanmax(d),
                          float(np.nanmedian(np.abs(np.diff(d)))) if len(d) > 1 else np.nan, len(d), [c(x) for x in names])
                absent = fs.absent_value
                for s in fs.specs:
                    m = c(s.mnemonic)
                    if m in a.dtype.names and a.dtype[m].shape == () and a.dtype[m].kind in "if":
                        v = a[m].astype(float); v[(v == absent) | (v == -999.25)] = np.nan
                        curves.append(dict(well=well, file=os.path.basename(p), frame=key, mnem=m, units=c(s.units),
                                           service=c(s.service_id), pct_valid=round(100 * np.isfinite(v).mean(), 1),
                                           p05=np.nanpercentile(v, 5) if np.isfinite(v).any() else np.nan,
                                           p50=np.nanmedian(v) if np.isfinite(v).any() else np.nan,
                                           p95=np.nanpercentile(v, 95) if np.isfinite(v).any() else np.nan))

def do_dlis(well, p):
    with dlis.load(p, error_handler=H) as lf:
        for li, f in enumerate(lf):
            for fr in f.frames:
                key = f"L{li}:{fr.name}"
                chans = fr.channels
                if not chans: continue
                idx = chans[0]
                names = [ch.name for ch in chans[1:] if (ch.dimension in ([1], [])) ]
                arr_names = [ch.name for ch in chans[1:] if ch.dimension not in ([1], [])]
                u = c(idx.units)
                top, base = fr.index_min, fr.index_max
                try: top, base = float(top), float(base)
                except Exception: top = base = None
                add_frame(well, p, "DLIS", key, idx.name, u, top, base, fr.spacing, None, names,
                          note=("arrays: " + " ".join(arr_names[:15])) if arr_names else "")
                for ch in chans[1:]:
                    curves.append(dict(well=well, file=os.path.basename(p), frame=key, mnem=ch.name, units=c(ch.units),
                                       service="", pct_valid=np.nan, p05=np.nan, p50=np.nan, p95=np.nan,
                                       long_name=c(ch.long_name), dim=str(ch.dimension)))

def do_las(well, p):
    l = lasio.read(p, ignore_header_errors=True, engine="normal")
    d = l.index.astype(float)
    u = c(l.curves[0].unit) or "M"
    names = [cv.mnemonic for cv in l.curves[1:]]
    add_frame(well, p, "LAS", "las", l.curves[0].mnemonic, u, np.nanmin(d), np.nanmax(d),
              float(np.nanmedian(np.abs(np.diff(d)))) if len(d) > 1 else np.nan, len(d), names)
    for cv in l.curves[1:]:
        v = np.asarray(cv.data, dtype=float) if cv.data.dtype.kind in "if" else np.array([np.nan])
        curves.append(dict(well=well, file=os.path.basename(p), frame="las", mnem=cv.mnemonic, units=c(cv.unit),
                           service="", pct_valid=round(100 * np.isfinite(v).mean(), 1),
                           p05=np.nanpercentile(v, 5) if np.isfinite(v).any() else np.nan,
                           p50=np.nanmedian(v) if np.isfinite(v).any() else np.nan,
                           p95=np.nanpercentile(v, 95) if np.isfinite(v).any() else np.nan, long_name=cv.descr))

files = sorted(p for p in glob.glob(f"{ROOT}/*/Raw_logs/**/*", recursive=True) if os.path.isfile(p))
skipped = []
for p in files:
    well = p.split("/")[4]
    fmt = sniff(p)
    try:
        if fmt == "LIS": do_lis(well, p)
        elif fmt == "DLIS": do_dlis(well, p)
        elif fmt == "LAS": do_las(well, p)
        else: skipped.append((well, os.path.basename(p), fmt))
    except Exception as e:
        skipped.append((well, os.path.basename(p), f"{fmt} ERROR {e}"))
    print(f"{well:6s} {fmt:5s} {os.path.basename(p)}", flush=True)

pd.DataFrame(frames).to_csv(f"{OUT}/all_frames.csv", index=False)
pd.DataFrame(curves).to_csv(f"{OUT}/all_curves.csv", index=False)
pd.DataFrame(skipped, columns=["well", "file", "reason"]).to_csv(f"{OUT}/skipped.csv", index=False)
print("frames", len(frames), "curves", len(curves), "skipped", len(skipped))
