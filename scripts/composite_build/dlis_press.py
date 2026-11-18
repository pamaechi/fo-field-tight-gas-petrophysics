import sys, glob, numpy as np, pandas as pd, logging, warnings
sys.path.insert(0,'work'); warnings.filterwarnings('ignore'); logging.getLogger('dlisio').setLevel(50)
from dlisio import dlis
from rft_extract import interpret
TRACES={}
def stations(well, pattern, pcols=('QCP','HPGP','BQP','PQP','SGP','CSGP')):
    rows=[]
    for x in sorted(glob.glob(pattern)):
        try:
            with dlis.load(x) as F:
                for f in F:
                    ps={p.name: p.values for p in f.parameters}
                    for fr in f.frames:
                        ch=[c.name for c in fr.channels]
                        if 'TIME' not in ch and 'ETIM' not in ch: continue
                        if not any(p in ch for p in pcols): continue
                        a=fr.curves(); units={c.name:(c.units or '') for c in fr.channels}
                        tc='ETIM_XPT' if 'ETIM_XPT' in ch else ('ETIM' if 'ETIM' in ch else 'TIME')
                        t=np.asarray(a[tc],float)/(1000.0 if units.get(tc,'').lower()=='ms' else 1.0)
                        dep=None
                        for dc in ('PDEP_XPT','RFPD'):
                            if dc in ch:
                                v=np.asarray(a[dc],float); v=v[np.isfinite(v)&(v>100)]
                                if len(v): dep=float(np.median(v)); break
                        if dep is None and 'FDEP' in ps and len(ps['FDEP']): dep=float(ps['FDEP'][0])
                        unitfix=1.0
                        if dep and dep>20000: dep=dep*0.00254   # 0.1 in
                        for pc in pcols:
                            if pc in ch:
                                p=np.asarray(a[pc],float)
                                if np.nanmax(p)>500 and np.nanstd(p)>1:
                                    unit=[c.units for c in fr.channels if c.name==pc][0]
                                    r=interpret(t,p)
                                    if r: rows.append(dict(well=well,file=x.split('/')[-1],frame=fr.name,depth_md=dep,gauge=pc,gauge_unit=unit,**r))
                                    TRACES[(well,x.split('/')[-1])]=(t,p,r)
                                    break
        except Exception as e:
            rows.append(dict(well=well,file=x.split('/')[-1],quality=f'read error {e}'))
    return pd.DataFrame(rows)
if __name__=='__main__':
    pd.set_option('display.width',250)
    a=stations('F-O4','data/raw/F-O4/Raw_logs/STA_*'); print(a.drop(columns=['file']).round(2).to_string()); a.to_csv('rft_F-O4.csv',index=False)
    b=stations('F-O8','data/raw/F-O8/Raw_logs/XPT/*LTP*'); print(b.drop(columns=['file']).round(2).to_string()); b.to_csv('rft_F-O8.csv',index=False)
