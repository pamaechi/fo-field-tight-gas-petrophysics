import sys; sys.path.insert(0,'work')
import pandas as pd, glob, os
from qcplot import qc_plot
M={"F-O4":[3700.0,3725.5,3757.0,3771.5,3787.0,3799.0,3821.0,3842.5,3859.0,3876.0],"F-O1":[3695.0,None,3702.5,3715.0,3739.5,3750.5,3771.5,3792.5,3811.5,3815.5],
   "F-O2":[3595.0,None,3600.0,3605.0,3615.5,3627.5,3648.5,3668.0,3685.5,3699.5],"F-O3":[3678.0,3688.5,3702.0,3707.5,3714.5,3727.0,3746.0,3763.5,3780.0,3794.0],
   "F-R1":[3836.0,3882.5,3913.0,3928.0,3940.5,3951.0,3972.5,3992.5,4003.5,4035.0],"F-S1":[3596.0,None,3605.0,3614.0,3620.5,3633.5,3649.0,3664.0,3678.0,3712.5]}
LK={"F-O1":(3709.63,3715.71),"F-O2":(3615.28,3628.02),"F-O4":(3709.83,3729.9),"F-O6":(3708,3723.14),"F-O8":(3708.3,3713.93)}
N=["M10 BLVS","M9","M8","M7","M6","M5","M4","M3","M2","M1"]
for f in sorted(glob.glob('work/out/composites/*_composite_v2.csv')):
    w=os.path.basename(f)[:4]; d=pd.read_csv(f,index_col=0); ren={}
    for tgt,opts in (('CALI',('LCAL','CALS')),('LLD',('RLA5','AT90','ILD')),('LLS',('RLA3','AT30','ILM','SFLU')),('DT',('DTCO',))):
        if tgt not in d:
            for o in opts:
                if o in d: ren[o]=tgt; break
    mk={v:n for v,n in zip(M.get(w,[]),N) if v} or None
    top,base=(3800,4060) if w=='F-R1' else (3560,3960)
    qc_plot(d.rename(columns=ren),f'work/out/composites/{w}_USM_QC.png',f'{w} · reservoir · composite v2 (cleaned)',top,base,mk,[(LK[w][0],LK[w][1],"Low-K unit","#eda100")] if w in LK else None,None)
