import pandas as pd,numpy as np
V='outputs/nb7test/out/{}_composite_v7.csv'
P={'F-O1':(3695,3815.5,3823),'F-O2':(3595,3699.5,3822),'F-O3':(3678,3794,3819),'F-O4':(3700,3876,3818),'F-O6':(3700,3930,3823.45),'F-O8':(3700,3916,3822),'F-R1':(3836,4035,3827),'F-S1':(3596,3712.5,3819)}
rows=[];allr=[]
for w,(top,base,gwc) in P.items():
    d=pd.read_csv(V.format(w),index_col=0)
    n='NPHI' if d.NPHI.notna().sum()>100 else ('TNPH' if 'TNPH' in d and d.TNPH.notna().sum()>100 else 'NPOR')
    d['DS']=d[n]-d.PHIT
    ok=(d.BADHOLE!=1)&(d.VSH_GR<0.15)&(d.PHIT>0.08)&d.DS.notna()
    win=(d.index>=top)&(d.index<=max(base,gwc+80))
    gas=ok&win&(d.index<gwc-5)&(d.index>=top)
    wat=ok&(d.index>gwc+5)&(d.index<=gwc+150)
    for lab,m in [('gas leg',gas),('water leg',wat)]:
        s=d[m]
        if len(s)>=10:
            rows.append(dict(well=w,zone=lab,n=len(s),phit=s.PHIT.median(),vsh=s.VSH_GR.median(),dsep_med=s.DS.median(),frac_cross=(s.DS<-0.02).mean(),curve=n))
            allr.append(s.assign(well=w,zone=lab)[['well','zone','PHIT','VSH_GR','DS']])
R=pd.DataFrame(rows); print(R.round(3).to_string(index=False))
A=pd.concat(allr); A.to_csv('resgas_samples.csv')
# matched porosity bins
A['pb']=pd.cut(A.PHIT,[0.08,0.10,0.12,0.14,0.2])
print(A.groupby(['zone','pb']).DS.agg(['median','size']).round(3))
R.to_csv('resgas_summary.csv',index=False)
