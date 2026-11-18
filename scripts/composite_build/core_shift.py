import pandas as pd, numpy as np
COMP='work/out/composites/{}_composite_v2.csv'
def core_groups(c):
    w=c.well.iloc[0]
    if w in ('F-O2','F-O3','F-O4') and c.core.notna().all(): return c.core.astype(int)
    return c.core_seq.astype(int)
def scan(core_df, log, curve, smin=-8, smax=8, step=0.0762, minn=6):
    z=core_df.depth_m.values; p=core_df.phi_pct.values; ok=np.isfinite(p)
    z,p=z[ok],p[ok]
    if len(z)<minn: return None
    L=log[curve].dropna()
    res=[]
    for s in np.arange(smin,smax+1e-9,step):
        v=np.interp(z+s,L.index.values,L.values,left=np.nan,right=np.nan); m=np.isfinite(v)
        if m.sum()<minn: res.append((s,np.nan)); continue
        res.append((s,np.corrcoef(p[m],-v[m])[0,1]))
    r=pd.DataFrame(res,columns=['shift','r']).dropna()
    if r.empty: return None
    b=r.loc[r.r.idxmax()]
    # second-best peak separated by > 1.5 m (ambiguity check)
    far=r[(r['shift']-b['shift']).abs()>1.5]
    return dict(best_shift=round(b['shift'],2), r_best=round(b.r,3), r_zero=round(r.loc[(r['shift']).abs().idxmin(),'r'],3),
                second_r=round(far.r.max(),3) if len(far) else np.nan, n=int(len(z)), curve=curve)
if __name__=='__main__':
    a=pd.read_csv('work/core_all_raw.csv')
    rows=[]
    for w,c in a.groupby('well'):
        log=pd.read_csv(COMP.format(w),index_col=0); c=c.copy(); c['grp']=core_groups(c)
        curves=['RHOB']+(['RHOB_ALT'] if 'RHOB_ALT' in log else [])
        for g,cc in c.groupby('grp'):
            for cv in curves:
                r=scan(cc,log,cv)
                if r: rows.append(dict(well=w,core=g,top=cc.depth_m.min(),base=cc.depth_m.max(),**r))
    R=pd.DataFrame(rows); pd.set_option('display.width',200); print(R.to_string(index=False)); R.to_csv('work/core_shift_scan.csv',index=False)
