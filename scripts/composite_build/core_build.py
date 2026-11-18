import sys; sys.path.insert(0,'work')
import pandas as pd, numpy as np
from core_parse import *
O='work/ocr/'
def build():
    fo1=parse_fo1_appJ(open(O+'F-O1_Well_Completion_Appendix_J.txt').read()); fo1=fo1[fo1.overburden_psi==800].copy()
    bad=fo1.kkl_md>fo1.kair_md*1.02; fo1.loc[bad,['kair_md','kkl_md']]=np.nan; fo1.loc[bad,'flag']='perm misread (Kkl>Kair) - nulled'
    fo1.loc[fo1.phi_pct<2,'phi_pct']=np.nan
    fo2=apply_manual(parse_layout_full(open(O+'F-O2_Other_Conventional_Core_Analysis_on_Borehole_F-O2.txt').read(),'F-O2'))
    fo3=parse_layout_cms(open(O+'F-O3_Other_Conventional_Core_Analysis_on_Borehole_F-O3.txt').read(),'F-O3')
    fs1=parse_layout_cms(open(O+'F-S1_Other_Conventional_Core_Analysis_on_Borehole_F-S1.txt').read(),'F-S1')
    fo6=pd.read_csv('work/core_F-O6_TerraTek.csv').assign(well='F-O6',kkl_md=lambda x:x.kgas_md)
    parts=[fo1,fo2,fo3,fs1,fo6]
    import os
    if os.path.exists('work/core_F-O4.csv'): parts.append(pd.read_csv('work/core_F-O4.csv'))
    out=[]
    for w,g in pd.concat(parts,ignore_index=True).groupby('well'):
        g=g.sort_values('depth_m').copy(); gap=g.depth_m.diff().fillna(0)>2.0; g['core_seq']=gap.cumsum()+1; out.append(g)
    return pd.concat(out,ignore_index=True)
if __name__=='__main__':
    a=build(); a.to_csv('work/core_all_raw.csv',index=False)
    print(a.groupby(['well','core_seq']).depth_m.agg(['min','max','count']))
