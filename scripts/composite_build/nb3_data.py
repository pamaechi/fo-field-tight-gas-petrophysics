import pandas as pd, numpy as np
V='outputs/nbtest/out/'
WELLS=['F-O1','F-O2','F-O3','F-O4','F-O6','F-O8','F-R1','F-S1']
def logs():
    out={}
    for w in WELLS:
        d=pd.read_csv(V+f'{w}_composite_v3.csv',index_col=0)
        d['RHOB_USE']=d['RHOB_RAW'] if 'RHOB_RAW' in d else d['RHOB']
        if w=='F-O8' and 'POTA' in d: d['POTA']=d['POTA']*100
        if w=='F-O6' and 'HFK' in d: d['POTA']=d['HFK']*100; d['THOR']=d['HTHO']; d['URAN']=d['HURA']
        out[w]=d
    return out
def core_with_logs(L, curves=('GR','RHOB_USE','NPHI','PEF','DT','DTCO','BADHOLE','THOR','POTA','KTIM','KSDR','TCMR','CMFF'), win=1):
    c=pd.read_csv(V+'core_data_shifted.csv'); rows=[]
    for w,g in c.groupby('well'):
        d=L[w]; z=d.index.values
        for _,r in g.iterrows():
            i=np.searchsorted(z,r.depth_log_m); i=min(max(i,win),len(z)-win-1)
            seg=d.iloc[i-win:i+win+1]
            rows.append({**r.to_dict(), **{cv: seg[cv].mean() if cv in d else np.nan for cv in curves}})
    return pd.DataFrame(rows)
