"""Build a spliced, depth-regular composite from a list of (file-glob, frame, {curve: new_name}, top, base) segments."""
import glob, os, numpy as np, pandas as pd
from .logio import lis_frames, dlis_frames, las_frame

def load_frame(raw_dir, well, pattern, frame):
    hits = sorted(glob.glob(os.path.join(str(raw_dir), pattern)))
    if not hits: raise FileNotFoundError(f"no file matching {pattern} in {raw_dir}")
    p = hits[0]
    for k, df, u in lis_frames(p, well):
        if k == frame:
            return df, u, p
    raise KeyError(f"{pattern} {frame} not found")

def build(raw_dir, well, segments, step=0.1524, top=None, base=None):
    tops = [s[3] for s in segments]; bases = [s[4] for s in segments]
    top = top if top is not None else min(tops); base = base if base is not None else max(bases)
    grid = np.round(np.arange(top, base + step / 2, step), 4)
    comp = pd.DataFrame(index=pd.Index(grid, name="DEPTH_M"))
    units, prov = {}, []
    for pattern, frame, cmap, t, b in segments:
        df, u, p = load_frame(raw_dir, well, pattern, frame)
        df = df[~df.index.duplicated()].loc[t:b]
        for src, dst in cmap.items():
            if src not in df: continue
            s = df[src]
            vals = np.interp(grid, s.index.values, s.values, left=np.nan, right=np.nan)
            vals[(grid < max(t, s.index.min())) | (grid > min(b, s.index.max()))] = np.nan
            # do not interpolate across gaps larger than 2 samples
            valid = s.notna()
            gi = pd.Series(valid.values.astype(float), index=s.index)
            vmask = np.interp(grid, gi.index.values, gi.values, left=0, right=0) > 0.99
            vals[~vmask] = np.nan
            if dst not in comp: comp[dst] = np.nan
            fill = comp[dst].isna() & ~np.isnan(vals)
            comp.loc[fill, dst] = vals[fill]
            units[dst] = u.get(src, "")
            prov.append(dict(curve=dst, source_file=os.path.basename(p), frame=frame, source_mnem=src, top=t, base=b))
    return comp, units, pd.DataFrame(prov)
