import sys; sys.path.insert(0,'work')
import pandas as pd, numpy as np
from core_shift import core_groups, COMP
PUBLISHED = {"F-O1": "RPS 2008: +3.5 m; SOEKOR: +2 m", "F-O2": "RPS 2008: +4.5 m; SOEKOR: cores 1-4 +1 m, core 5 +4.5 m",
             "F-O3": "RPS 2008: +2.75 m", "F-O4": "RPS 2008: +2.75 m"}
OVERRIDE = {  # multi-log consensus where porosity contrast is too low for a single-curve match (evidence in Notebook 02)
    ("F-O4", 1): (2.95, "multi-log consensus: RHOB +3.05 (r 0.51), DTCO +2.90 (r 0.47), DT +2.82 (r 0.44)"),
    ("F-O4", 3): (2.70, "multi-log consensus: RHOB +2.82 (r 0.59), DTCO +2.52 (r 0.57), DT +2.52 (r 0.57)"),
    ("F-O4", 2): (2.83, "mean of F-O4 cores 1 and 3 (core 2 too homogeneous: r < 0.23 on all logs)"),
}
def finalize():
    a=pd.read_csv('work/core_all_raw.csv'); S=pd.read_csv('work/core_shift_scan.csv'); S=S[S.curve=='RHOB']
    S['robust']=(S.r_best>=0.75)&((S.r_best-S.second_r).fillna(1)>=0.10)&(S.n>=6)
    med=S[S.robust].groupby('well').best_shift.median()
    out=[]; table=[]
    for w,c in a.groupby('well'):
        c=c.copy(); c['grp']=core_groups(c)
        for g,cc in c.groupby('grp'):
            s=S[(S.well==w)&(S.core==g)]
            if (w,g) in OVERRIDE: sh,how=OVERRIDE[(w,g)]
            elif len(s) and bool(s.robust.iloc[0]): sh=float(s.best_shift.iloc[0]); how=f"core-log correlation (r={s.r_best.iloc[0]:.2f}, n={int(s.n.iloc[0])})"
            elif len(s) and w not in med.index and s.r_best.iloc[0]>=0.6 and (s.r_best.iloc[0]-s.second_r.iloc[0])>=0.15:
                sh=float(s.best_shift.iloc[0]); how=f"core-log correlation, MODERATE confidence (r={s.r_best.iloc[0]:.2f}, 2nd peak {s.second_r.iloc[0]:.2f}; only core in well)"
            else: sh=float(med.get(w,np.nan)); how="well median of robust cores (this core unconstrained: " + (f"r={s.r_best.iloc[0]:.2f}, 2nd peak {s.second_r.iloc[0]:.2f}" if len(s) else "too few plugs") + ")"
            cc['shift_m']=sh; cc['depth_log_m']=cc.depth_m+sh; out.append(cc)
            table.append(dict(well=w,core=g,core_top_m=cc.depth_m.min(),core_base_m=cc.depth_m.max(),plugs=len(cc),shift_m=round(sh,2),method=how,published=PUBLISHED.get(w,"")))
    return pd.concat(out,ignore_index=True), pd.DataFrame(table)
if __name__=='__main__':
    c,t=finalize(); c.to_csv('work/core_all_shifted.csv',index=False); t.to_csv('work/core_shift_table.csv',index=False)
    pd.set_option('display.width',250); pd.set_option('display.max_colwidth',80); print(t.to_string(index=False))
