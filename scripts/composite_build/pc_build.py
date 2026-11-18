import re, json, numpy as np, pandas as pd
SIG_AB, SIG_HG, SIG_GW = 72.0, 485.0*abs(np.cos(np.radians(140))), 35.0   # sigma*cos(theta), dyn/cm: air-brine, air-Hg, gas-brine (reservoir)
rows=[]
# F-O6 MICP
m=pd.read_csv('work/micp_F-O6.csv')
for _,r in m.iterrows():
    rows.append(dict(well='F-O6',method='MICP',sample=r['sample'],depth_m=r.depth_m,k_md=r.k_md,phi=r.phi,pc_lab_psi=r.pc_hg_psia,sigcos_lab=SIG_HG,sw=1-r.shg))
# F-O6 porous plate (air-brine)
PP6=json.load(open('data/scal/porous_plate_F-O6.json'))   # {sample: [depth_m, k_md, phi, [Sw at each P6 step]]} from the F-O6 SCAL report (not distributed)
P6=[0,1,2,4,8,15,35,60]
for s,(z,k,phi,sw) in PP6.items():
    for pc,v in zip(P6[1:],sw[1:]): rows.append(dict(well='F-O6',method='porous plate',sample=s,depth_m=z,k_md=k,phi=phi,pc_lab_psi=pc,sigcos_lab=SIG_AB,sw=v))
# F-O2 porous plate (air-brine), Sw % right-aligned to pressures 4..180 psi; porosity from CEC table; k from nearest CCA plug
PP2=json.load(open('data/scal/porous_plate_F-O2.json'))   # {sample: [depth_m, phi_percent, [Sw % right-aligned to P2]]} from the F-O2 SCAL report (not distributed)
P2=[4,8,15,35,50,80,180]
core=pd.read_csv('data/core_data_raw_unshifted.csv'); c2=core[(core.well=='F-O2')&(core.kkl_md>0)]
for s,(z,phi,sw) in PP2.items():
    j=(c2.depth_m-z).abs(); k=float(c2.kkl_md[j.idxmin()]) if j.min()<0.3 else np.nan
    for pc,v in zip(P2[-len(sw):],sw): rows.append(dict(well='F-O2',method='porous plate',sample=s,depth_m=z,k_md=k,phi=phi/100,pc_lab_psi=pc,sigcos_lab=SIG_AB,sw=v/100,k_source='nearest CCA plug (<0.3 m)'))
D=pd.DataFrame(rows)
D['pc_res_psi']=D.pc_lab_psi*SIG_GW/D.sigcos_lab
D['J']=0.2166*D.pc_lab_psi/D.sigcos_lab*np.sqrt(D.k_md/D.phi)
D.to_csv('work/pc_all.csv',index=False)
print(D.groupby(['well','method']).agg(samples=('sample','nunique'),points=('sw','size'),k_min=('k_md','min'),k_max=('k_md','max')).round(3))
