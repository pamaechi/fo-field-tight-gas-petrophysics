import glob, os, sys, numpy as np, pandas as pd, warnings, logging
warnings.filterwarnings('ignore'); logging.getLogger('dlisio').setLevel(logging.CRITICAL)
from logio import lis_frames, dlis_frames, las_frame
KEY=['GR','SGR','CGR','ECGR','GR_EDTC','POTA','THOR','URAN','HFK','HTHO','HURA','SP','CALI','CALS','HCAL','DCAL','LCAL','C1','C2','LLD','LLS','MSFL','ILD','ILM','SFLU','RLA5','RLA3','RLA1','RXOZ','AT90','AT60','AT30','AT10','AF90','AO90','RHOB','RHOZ','RHOM','DRHO','HDRA','NPHI','TNPH','NPOR','PEF','PEFZ','PEFL','DT','DTL','DTCO','DTSM','DT4P','DT4S','TCMR','CMFF','BFV','KSDR','KTIM','T2LM']
well=sys.argv[1]; rows=[]
def sniff(p):
    b=open(p,'rb').read(100)
    if b'RECORD' in b[:90]: return 'DLIS'
    if b.lstrip()[:2] in (b'~V',b'~v'): return 'LAS'
    if b[:4]==b'\0\0\0\0' or (b and b[0] in (0x84,0x80,0x82)): return 'LIS'
    return None
for p in sorted(glob.glob(f'data/raw/{well}/Raw_logs/*')):
    if not os.path.isfile(p) or os.path.getsize(p)>150e6: continue
    f=sniff(p)
    try:
        it = lis_frames(p,well) if f=='LIS' else (dlis_frames(p) if f=='DLIS' else ([('las',)+las_frame(p)] if f=='LAS' else []))
        for k,df,u in it:
            ok=np.isfinite(df.index.values)
            if ok.sum()<50 or df.index[ok].max()>6000 or df.index[ok].min()<0: continue
            have=[c for c in KEY if c in df and df[c].notna().sum()>50]
            if not have: continue
            rows.append(dict(file=os.path.basename(p)[:40],fmt=f,frame=k[:22],top=round(df.index[ok].min(),1),base=round(df.index[ok].max(),1),step=round(float(np.median(np.abs(np.diff(df.index[ok])))),4),n=int(ok.sum()),key=' '.join(have)))
    except Exception as e:
        rows.append(dict(file=os.path.basename(p)[:40],fmt=f,frame='ERR',key=str(e)[:80]))
d=pd.DataFrame(rows); d.to_csv(f'inventory/{well}_frames.csv',index=False)
pd.set_option('display.width',260); pd.set_option('display.max_colwidth',130); pd.set_option('display.max_rows',300)
print(d.to_string(index=False))
