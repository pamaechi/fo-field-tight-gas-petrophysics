import pandas as pd, numpy as np
def bias(rhof=1.0, rows=None):
    c=pd.read_csv('work/core_all_shifted.csv'); out=[]
    for w,g in c.groupby('well'):
        log=pd.read_csv(f'work/out/composites/{w}_composite_v2.csv',index_col=0); g=g[g.phi_pct.notna()]
        for cv in [x for x in ('RHOB','RHOB_ALT') if x in log]:
            L=log[cv].dropna(); rb=np.interp(g.depth_log_m,L.index,L.values,left=np.nan,right=np.nan)
            rg=g.rhog.fillna(g.rhog.median() if g.rhog.notna().any() else 2.67).values; phid=100*(rg-rb)/(rg-rhof); m=np.isfinite(phid)
            d=phid[m]-g.phi_pct.values[m]
            out.append(dict(well=w,curve=cv,n=int(m.sum()),bias_mean_pu=round(d.mean(),2),bias_median_pu=round(float(np.median(d)),2),rmse_pu=round(float(np.sqrt(np.mean(d**2))),2),r=round(np.corrcoef(phid[m],g.phi_pct.values[m])[0,1],2)))
    return pd.DataFrame(out)
if __name__=='__main__': print(bias().to_string(index=False))
