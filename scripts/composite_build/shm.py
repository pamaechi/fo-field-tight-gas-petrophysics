import pandas as pd, numpy as np
SIG_GW=35.0; DRHO_GRAD=1.024   # psi/m: (rho_w 0.95 - rho_g 0.229) g/cc x 1.422
KB={"F-O1":26.0,"F-O2":25.0,"F-O3":22.0,"F-O4":24.0,"F-O6":26.45,"F-O8":25.0,"F-S1":22.0,"F-R1":22.0}
WIN={"F-O1": (3695.0, 3815.5), "F-O2": (3595.0, 3699.5), "F-O3": (3678.0, 3794.0), "F-O4": (3700.0, 3876.0),"F-O6": (3700.0, 3930.0), "F-O8": (3700.0, 3916.0), "F-R1": (3836.0, 4035.0), "F-S1": (3596.0, 3712.5)}
def fit_J(D):
    D=D[(D.sw<0.98)&(D.sw>0.02)&(D.J>0)].copy(); D['lj']=np.log10(D.J); D['ls']=np.log10(D.sw); D['bin']=np.round(D.lj*4)/4
    B=D.groupby(['sample','bin']).agg(lj=('lj','mean'),ls=('ls','mean')).reset_index()
    p=np.polyfit(B.lj,B.ls,1); sd=float((B.ls-np.polyval(p,B.lj)).std()); return 10**p[1], p[0], sd
def sw_shm(h_m, k_md, phi, a, b):
    pc=np.clip(h_m,0,None)*DRHO_GRAD; J=0.2166*pc/SIG_GW*np.sqrt(np.clip(k_md,1e-4,None)/np.clip(phi,0.01,None))
    with np.errstate(divide='ignore'): sw=np.where(J>0, a*np.power(J,b), 1.0)
    return np.clip(sw,0,1)
def fit_fwl(d, well, a, b, grid=np.arange(3740,3900.1,1.0)):
    t,bw=WIN[well]; s=d.loc[t:bw]; s=s[(s.PHIT>0.05)&(s.VSH_GR<0.4)&s.SW.notna()&s.PERM.notna()]
    z=s.index.values-KB[well]; best=None; curve=[]
    for f in grid:
        pred=sw_shm(f-z, s.PERM.values, s.PHIT.values, a, b); err=np.median(np.abs(pred-s.SW.values)); curve.append((f,err))
        if best is None or err<best[1]: best=(f,err)
    return best, pd.DataFrame(curve,columns=['fwl_tvdss','mae']), s, z
