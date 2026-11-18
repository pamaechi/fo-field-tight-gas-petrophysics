import glob, os, numpy as np, pandas as pd, logging, warnings
warnings.filterwarnings('ignore'); logging.getLogger('dlisio').setLevel(50)
from dlisio import lis
REP=os.path.expanduser('~/.fo_repaired_lis')
def consts(lf):
    out={}
    for wb in lf.wellsite_data():
        try: comps=wb.components()
        except Exception: continue
        cur=None
        for c in comps:
            mn=str(c.mnemonic).strip()
            if mn=='MNEM': cur=str(c.component).strip()
            elif mn=='VALU' and cur and cur not in out: out[cur]=c.component
    return out
def interpret(t, p, hyd_tol=15.0, stable_rate=0.05, min_run_s=15.0, return_trace=False):
    """Formation pressure = median of the longest STABLE run (|dp/dt| < stable_rate psi/s on a 5-sample median-smoothed
    trace, sustained >= min_run_s seconds) that lies in the below-hydrostatic segment and above pmin + 50 % of the drawdown.
    Flags: good; supercharged? (stable run within 30 psi of mud hydrostatic); tight / not stabilised (no such run);
    no drawdown (seal failure)."""
    p=np.asarray(p,float); t=np.asarray(t,float); ok=np.isfinite(p)&np.isfinite(t)&(p>100); p,t=p[ok],t[ok]
    if len(p)<20: return None
    o=np.argsort(t); p,t=p[o],t[o]
    n=max(5,len(p)//20); hyd0=float(np.median(p[:n])); hyd1=float(np.median(p[-n:])); hyd=min(hyd0,hyd1)
    below=np.where(p<hyd-hyd_tol)[0]
    out=dict(hydro_before=hyd0, hydro_after=hyd1, pform=np.nan, pmin=np.nan, stable_s=0.0, quality='no drawdown (seal failure?)')
    if len(below)<5: return out
    i0,i1=below[0],below[-1]; pp=pd.Series(p[i0:i1+1]).rolling(5,center=True,min_periods=1).median().values; tt=t[i0:i1+1]
    pmin=float(pp.min()); dd=hyd-pmin
    rate=np.abs(np.gradient(pp,tt)) if len(tt)>2 else np.full(len(pp),np.inf)
    stable=(rate<stable_rate)&(pp>pmin+0.5*dd)
    best=(0,None,None); k=0
    while k<len(stable):
        if stable[k]:
            j=k
            while j+1<len(stable) and stable[j+1]: j+=1
            dur=tt[j]-tt[k]
            if dur>best[0]: best=(dur,k,j)
            k=j+1
        else: k+=1
    out.update(pmin=pmin)
    if best[1] is None or best[0]<min_run_s:
        out.update(pform=float(pp[-1]) if len(pp) else np.nan, quality='tight / not stabilised'); return out
    pf=float(np.median(pp[best[1]:best[2]+1]))
    out.update(pform=pf, stable_s=float(best[0]), quality='supercharged? (stable ~ mud pressure)' if hyd-pf<30 else 'good')
    return out

# --- readers used to build pressures_all.csv (kept for the record; need the raw tapes) -------------------------------------
# F-O1 (LIS, SGP strain gauge, test depth FDEP per logical file) and F-O3 (LIS, RPQP quartz, depth RFPD): see Notebook 05 notes.
def dlis_station(path, pcols=("QCP", "HPGP", "BQP", "PQP", "SGP", "CSGP")):
    """F-O4 STA_*.DLIS (HPGP psi, ETIM s, FDEP parameter) and F-O8 XPT *LTP*.DLIS (QCP psi, ETIM_XPT s, PDEP_XPT m)."""
    from dlisio import dlis
    out = []
    with dlis.load(path) as F:
        for f in F:
            ps = {p.name: p.values for p in f.parameters}
            for fr in f.frames:
                ch = [c.name for c in fr.channels]; units = {c.name: (c.units or "") for c in fr.channels}
                if not any(p in ch for p in pcols): continue
                a = fr.curves(); tc = "ETIM_XPT" if "ETIM_XPT" in ch else ("ETIM" if "ETIM" in ch else "TIME")
                t = np.asarray(a[tc], float) / (1000.0 if units.get(tc, "").lower() == "ms" else 1.0)
                dep = None
                for dc in ("PDEP_XPT", "RFPD"):
                    if dc in ch:
                        v = np.asarray(a[dc], float); v = v[np.isfinite(v) & (v > 100)]
                        if len(v): dep = float(np.median(v)); break
                if dep is None and "FDEP" in ps and len(ps["FDEP"]): dep = float(ps["FDEP"][0])
                for pc in pcols:
                    if pc in ch and np.nanmax(a[pc]) > 500:
                        r = interpret(t, a[pc]); out.append(dict(depth_md=dep, gauge=pc, gauge_unit=units[pc], **(r or {}))); break
    return out
