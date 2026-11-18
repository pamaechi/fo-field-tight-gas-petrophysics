import sys; sys.path.insert(0,'work')
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from qcplot import C1, C2, INK, MUTED, GRID
def fig(well, out):
    d=pd.read_csv(f'work/out/composites_v3/{well}_composite_v3.csv',index_col=0)
    c=pd.read_csv('work/core_all_shifted.csv'); c=c[c.well==well]
    T=pd.read_csv('work/core_shift_table.csv'); T=T[T.well==well]
    top,base=c.depth_m.min()-8, c.depth_log_m.max()+8; s=d.loc[top:base]
    rhob='RHOB_RAW' if 'RHOB_RAW' in d else 'RHOB'
    phid=(100*(2.67-s[rhob])/1.67).mask(s.get('BADHOLE',0)==1)
    f,ax=plt.subplots(1,3,figsize=(7.2,9),sharey=True,gridspec_kw=dict(width_ratios=[1,1.6,1.4]))
    ax[0].plot(s.GR,s.index,color=C1,lw=1); ax[0].set_xlim(0,150); ax[0].set_title('GR (API)',fontsize=9,color=INK)
    ax[1].plot(phid,s.index,color=C1,lw=1,label=f'φD from {rhob} (2.67/1.0)')
    ax[1].scatter(c.phi_pct,c.depth_m,s=10,facecolors='none',edgecolors=MUTED,lw=0.8,label='Core, core depth')
    ax[1].scatter(c.phi_pct,c.depth_log_m,s=12,color=C2,zorder=3,label='Core, log depth')
    ax[1].set_xlim(25,0); ax[1].set_title('Porosity (%)',fontsize=9,color=INK)
    k=c[c.kkl_md>0]; ax[2].scatter(k.kkl_md,k.depth_log_m,s=12,color=C2)
    if 'KTIM' in s: ax[2].plot(s.KTIM.clip(lower=1e-4),s.index,color=C1,lw=1,label='NMR KTIM')
    ax[2].scatter([],[],s=12,color=C2,label='Core Kkl')
    ax[2].legend(fontsize=6.5,loc='lower left',frameon=False)
    ax[2].set_xscale('log'); ax[2].set_xlim(1e-3,1e3); ax[2].set_title('Core Klinkenberg k (mD)',fontsize=9,color=INK)
    for a in ax:
        a.grid(True,color=GRID,lw=0.6); a.tick_params(labelsize=7,colors=MUTED); [sp.set_color(GRID) for sp in a.spines.values()]
    ax[0].invert_yaxis(); ax[0].set_ylabel('Logger depth (m below KB)',fontsize=8,color=INK)
    ax[1].legend(fontsize=6.5,loc='lower left',frameon=False)
    sh='; '.join(f"core {int(r.core)}: {r.shift_m:+.2f} m" for _,r in T.iterrows())
    f.suptitle(f'{well} · core–log depth match (log = core + shift)\n{sh}',fontsize=9,color=INK,x=0.02,ha='left')
    f.tight_layout(rect=(0,0,1,0.95)); f.savefig(out,dpi=150); plt.close(f)
if __name__=='__main__':
    for w in sys.argv[1:] or ['F-O1','F-O2','F-O3','F-O6','F-S1']:
        fig(w,f'work/out/composites_v3/{w}_core_depth_match.png'); print(w)
