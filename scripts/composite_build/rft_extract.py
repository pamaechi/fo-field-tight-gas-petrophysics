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
def lis_tests(well, pattern):
    rows=[]
    for f in sorted(glob.glob(f'{REP}/{pattern}')):
        with lis.load(f) as files:
            for i,lf in enumerate(files):
                c=consts(lf); fdep=c.get('FDEP')
                try:
                    for fs in lf.data_format_specs():
                        a=lis.curves(lf,fs,strict=False)
                        names=[n.strip() for n in a.dtype.names]; a.dtype.names=names
                        tcol='ETIM' if 'ETIM' in names else ('TIME' if 'TIME' in names else names[0])
                        for pc in ('HPGP','HPRE','PQ','SGP','PRES','BQP'):
                            if pc in names and np.nanmax(a[pc])>100:
                                r=interpret(a[tcol],a[pc]); 
                                if r: rows.append(dict(well=well,file=os.path.basename(f),lf=i,depth_md=fdep,gauge=pc,**r))
                                break
                        break
                except Exception as e:
                    rows.append(dict(well=well,file=os.path.basename(f),lf=i,depth_md=fdep,quality=f'read error {e}'))
    return pd.DataFrame(rows)
if __name__=='__main__':
    d=lis_tests('F-O1','*re175i*'); pd.set_option('display.width',220); print(d.round(2).to_string())
