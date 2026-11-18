import glob, numpy as np, pandas as pd, warnings, logging
warnings.filterwarnings('ignore'); logging.getLogger('dlisio').setLevel(logging.CRITICAL)
from logio import lis_frames, dlis_frames, las_frame
_cache={}
def frame(well, pat, key):
    k=(well,pat,key)
    if k in _cache: return _cache[k]
    p=sorted(glob.glob(f'data/raw/{well}/Raw_logs/{pat}'))[0]
    b=open(p,'rb').read(100)
    it = dlis_frames(p) if b'RECORD' in b[:90] else ([('las',)+las_frame(p)] if b.lstrip()[:2] in (b'~V',b'~v') else lis_frames(p,well))
    for kk,df,u in it:
        if kk.startswith(key):
            df=df[~df.index.duplicated()].sort_index(); _cache[k]=(df,u); return df,u
    raise KeyError(key)
def lag(well, a, b, curve='GR', curve_b=None, maxlag_m=12, step=0.1524):
    A,_=frame(well,*a); B,_=frame(well,*b); cb=curve_b or curve
    lo=max(A.index.min(),B.index.min()); hi=min(A.index.max(),B.index.max())
    g=np.arange(lo,hi,step)
    x=np.interp(g,A.index,A[curve].values); y=np.interp(g,B.index,B[cb].values)
    m=np.isfinite(x)&np.isfinite(y); x=(x-np.nanmean(x[m])); y=(y-np.nanmean(y[m]))
    best=None
    for L in range(-int(maxlag_m/step),int(maxlag_m/step)+1):
        if L>=0: xx,yy=x[L:],y[:len(y)-L]
        else: xx,yy=x[:L],y[-L:]
        mm=np.isfinite(xx)&np.isfinite(yy)
        if mm.sum()<200: continue
        r=np.corrcoef(xx[mm],yy[mm])[0,1]
        if best is None or r>best[1]: best=(L,r)
    r0=np.corrcoef(x[m],y[m])[0,1]
    return dict(overlap=(round(lo,1),round(hi,1)), r_at_0=round(r0,3), best_lag_m=round(best[0]*step,3), r_best=round(best[1],3))
