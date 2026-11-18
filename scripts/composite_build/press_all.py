import sys, os, glob, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'work'); import warnings, logging; warnings.filterwarnings('ignore'); logging.getLogger('dlisio').setLevel(50)
import rft2, dlis_press
from logio import lis_frames
from rft_extract import interpret
# F-O3: quartz RPQP per frame, depth RFPD
def fo3():
    rows=[]
    f=glob.glob('data/raw/F-O3/Raw_logs/*re231e*')[0]
    for k,df,u in lis_frames(f,'F-O3'):
        if 'RPQP' not in df: continue
        tc='ETIM' if 'ETIM' in df and df.ETIM.notna().sum()>10 else 'TIME'
        if tc not in df: continue
        t=df[tc].values/(1000.0 if u.get(tc,'').upper()=='MS' else 1.0)
        dep=float(np.nanmedian(df.RFPD)) if 'RFPD' in df else np.nan
        r=interpret(t,df.RPQP.values)
        if r: rows.append(dict(well='F-O3',file=os.path.basename(f),frame=k,depth_md=dep,gauge='RPQP',gauge_unit=u.get('RPQP',''),**r)); rft2.TRACES[('F-O3',k)]=(t,df.RPQP.values,r)
    return pd.DataFrame(rows)
P=pd.concat([rft2.tests('F-O1','*re175i*'), fo3(),
             dlis_press.stations('F-O4','data/raw/F-O4/Raw_logs/STA_*'),
             dlis_press.stations('F-O8','data/raw/F-O8/Raw_logs/XPT/*LTP*')],ignore_index=True)
# gauge pressure (psig) -> psia
P['pform_psia']=np.where(P.gauge_unit.str.upper().str.contains('PSIG'), P.pform+14.7, P.pform)
P.to_csv('work/pressures_all.csv',index=False)
T={**rft2.TRACES, **dlis_press.TRACES}
n=len(T); cols=8; rows_=int(np.ceil(n/cols))
fig,ax=plt.subplots(rows_,cols,figsize=(cols*2.4,rows_*1.9))
for a,(key,(t,p,r)) in zip(ax.flat,T.items()):
    a.plot(t-np.nanmin(t),p,lw=0.6)
    if r and np.isfinite(r.get('pform',np.nan)): a.axhline(r['pform'],color='#eb6834',lw=0.8)
    a.set_title(f"{key[0]} {str(key[1])[-14:]}\\n{r['quality'][:22] if r else ''}",fontsize=6); a.tick_params(labelsize=5)
for a in ax.flat[n:]: a.axis('off')
plt.tight_layout(); plt.savefig('work/pressure_qc_traces.png',dpi=90)
print(P.groupby(['well','quality']).size())
