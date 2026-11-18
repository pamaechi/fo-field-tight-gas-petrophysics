import sys, os, glob, numpy as np, pandas as pd, logging, warnings
sys.path.insert(0,'work'); warnings.filterwarnings('ignore'); logging.getLogger('dlisio').setLevel(50)
from logio import lis_frames
from rft_extract import consts, interpret
from dlisio import lis
REP=os.path.expanduser('~/.fo_repaired_lis')
def fdeps(rep_path):
    out=[]
    with lis.load(rep_path) as files:
        for lf in files: out.append(consts(lf).get('FDEP'))
    return out
TRACES={}
def tests(well, raw_glob, pcols=('SGP','PRES','RPCO','RPQS','PQ','BQP')):
    rows=[]
    for f in sorted(glob.glob(f'data/raw/{well}/Raw_logs/{raw_glob}')):
        rp=os.path.join(REP,os.path.basename(f)+'.lis')
        frames=list(lis_frames(f,well))
        fd=fdeps(rp) if os.path.exists(rp) else []
        # one data frame per logical file (L<i>D1)
        for k,df,u in frames:
            try: i=int(k.split('D')[0][1:])
            except Exception: continue
            tcol='ETIM' if 'ETIM' in df else ('TIME' if 'TIME' in df else None)
            if tcol is None: continue
            df=df.sort_values(tcol)
            for pc in pcols:
                if pc not in df: continue
                if pc in df and np.nanmax(df[pc])>500:
                    tt=df[tcol].values/(1000.0 if u.get(tcol,'').upper()=='MS' else 1.0)
                    r=interpret(tt, df[pc].values)
                    TRACES[(well,f"{os.path.basename(f)}:{k}")]=(tt,df[pc].values,r)
                    if r: rows.append(dict(well=well,file=os.path.basename(f),lf=i,frame=k,depth_md=fd[i] if i<len(fd) else np.nan,gauge=pc,gauge_unit=u.get(pc,''),**r))
                    break
    return pd.DataFrame(rows)
if __name__=='__main__':
    pd.set_option('display.width',250)
    for w,g in [('F-O1','*re175i*')]:
        d=tests(w,g); print(d.drop(columns=['file']).round(2).to_string()); d.to_csv(f'work/rft_{w}.csv',index=False)
