import sys, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "data"); from petrotools import mineral as mn
S="outputs"; O=f"{S}/paperfigs"; import os; os.makedirs(O,exist_ok=True)
plt.rcParams.update({"font.size":9,"axes.titlesize":9,"axes.labelsize":9,"legend.fontsize":7.5,"savefig.dpi":300})
# ---- Fig 5: XRD vs logs (replicates Notebook 06 cell 5-6)
xrd=pd.read_csv("data/xrd_F-O6.csv"); core=pd.read_csv(f"{S}/nbtest/out/core_data_shifted.csv")
d6=pd.read_csv(f"{S}/nb6test/out/F-O6_composite_v6.csv",index_col=0); d6["DT"]=d6.DTCO; sh=float(core[core.well=="F-O6"].shift_m.iloc[0])
for c in ["RHOB_USE","NPHI","PEF","DT","VSH_GR","POTA"]: xrd[c]=[d6[c].loc[z+sh-0.3:z+sh+0.3].mean() for z in xrd.depth]
c6=core[core.well=="F-O6"]; xrd["phi_core"]=[c6.phi_pct[(c6.depth_m-z).abs()<0.35].mean()/100 if ((c6.depth_m-z).abs()<0.35).any() else np.nan for z in xrd.depth]
xrd["carb"]=xrd.cal+xrd.sid+xrd.ank+xrd.dol; sand=xrd.depth>3707
xrd["phi_cal"]=xrd.phi_core.fillna(0.03); vol=mn.xrd_to_volumes(xrd,xrd.phi_cal); RESP=mn.calibrate(vol,xrd)
nd=mn.invert(xrd[["NPHI","DT","VSH_GR"]],{k:v for k,v in RESP.items() if k!="RHOB_USE"},None); xrd["RHOMAA"]=mn.apparent_matrix_density(xrd.RHOB_USE,nd.PHI_MIN)
cf=(vol.CARB/(1-vol.PHI)); CARB=mn.fit_carbonate_index(xrd.RHOMAA[sand],cf[sand])
loo=[]
for i in range(len(xrd)):
    tr=np.arange(len(xrd))!=i; Ri=mn.calibrate(vol[tr].reset_index(drop=True),xrd[tr].reset_index(drop=True))
    loo.append(mn.invert(xrd.iloc[[i]][["RHOB_USE","NPHI","DT","VSH_GR"]],Ri,None).iloc[0])
loo=pd.DataFrame(loo)
fig,ax=plt.subplots(1,3,figsize=(7.2,2.6))
ax[0].scatter(xrd.clay[sand],xrd.POTA[sand],c="#1f77b4",s=18,label="sandstone"); ax[0].scatter(xrd.clay[~sand],xrd.POTA[~sand],c="#888",marker="s",s=18,label="mudstone")
ax[0].set(xlabel="XRD clay (wt%)",ylabel="K (%)",title="(a) clay vs potassium, r = 0.93"); ax[0].legend()
r=np.corrcoef(vol.CLAY,loo.V_CLAY)[0,1]; ax[1].scatter(100*vol.CLAY,100*loo.V_CLAY,s=18,c=np.where(sand,"#1f77b4","#888")); ax[1].plot([0,55],[0,55],"k-",lw=.7)
ax[1].set(xlabel="XRD clay (vol%)",ylabel="solver clay, sample withheld (vol%)",title=f"(b) clay solver, r = {r:.2f}")
x=np.linspace(2.55,2.80,10); ax[2].scatter(xrd.RHOMAA[sand],100*cf[sand],s=18,c="#1f77b4"); ax[2].plot(x,100*(CARB["a"]+CARB["b"]*x),"k--",lw=.8)
ax[2].set(xlabel="apparent matrix density (g/cm³)",ylabel="XRD carbonate (vol% of solids)",title=f"(c) carbonate cement, r = {CARB['r']:.2f}")
for a in ax: a.grid(alpha=.3)
plt.tight_layout(); plt.savefig(f"{O}/P1_Fig5_xrd_logs.png"); plt.close()
print("Fig5", round(r,2), CARB)
# ---- Fig 7: mineralogy by rock type (training plugs)
T=pd.read_csv("nb7/T.csv"); T["CLAY"]=100*T.V_CLAY; T["CARB"]=100*T.V_CARB
fig,ax=plt.subplots(1,4,figsize=(7.2,2.5))
for a,(c,l) in zip(ax,[("kkl_md","core k (mD)"),("VSH_GR","VSH_GR (v/v)"),("CARB","V_CARB (vol%)"),("DSEP","NPHI − PHIT (v/v)")]):
    a.boxplot([T[T.RT==k][c] for k in (1,2,3,4)],tick_labels=["RT1","RT2","RT3","RT4"],showfliers=False,widths=.6); a.set_title(l); a.grid(alpha=.3,axis="y")
ax[0].set_yscale("log")
plt.tight_layout(); plt.savefig(f"{O}/P1_Fig7_rocktype_logs.png"); plt.close()
print(T.groupby("RT").size().to_dict())
# ---- Fig 8: blind permeability per well
b=pd.read_csv(f"{S}/nb6test/out/blind_perm_rocktype_vs_field.csv"); x=np.arange(len(b)); w=.38
fig,ax=plt.subplots(1,2,figsize=(7.2,2.6))
ax[0].bar(x-w/2,b.field_bias_dex,w,color="#999",label="field φ–k transform"); ax[0].bar(x+w/2,b.rocktype_bias_dex,w,color="#1f77b4",label="rock-type transforms"); ax[0].axhline(0,color="k",lw=.7)
ax[1].bar(x-w/2,b.field_rmse_dex,w,color="#999"); ax[1].bar(x+w/2,b.rocktype_rmse_dex,w,color="#1f77b4")
for a,t in zip(ax,["(a) blind bias (log10 k, dex)","(b) blind RMSE (dex)"]): a.set_xticks(x); a.set_xticklabels(b.well); a.set_title(t); a.grid(alpha=.3,axis="y")
ax[0].legend(loc="lower right"); plt.tight_layout(); plt.savefig(f"{O}/P1_Fig8_blind_perm.png"); plt.close()
# ---- Fig 9: net pay
n=pd.read_csv(f"{S}/nb4test/out/net_pay_sensitivity.csv"); n=n[n.well!="F-R1"]; x=np.arange(len(n))
fig,ax=plt.subplots(figsize=(4.6,2.8))
ax.bar(x,n.net_pay_base_m,.6,color="#1f77b4",label="base (m* 2.44, n* 2.4)")
ax.errorbar(x,n.net_pay_base_m,yerr=[n.net_pay_base_m-n.net_pay_high_m,n.net_pay_low_m-n.net_pay_base_m],fmt="none",ecolor="k",capsize=3,lw=.8,label="m*/n* low–high")
ax.scatter(x,n.net_pay_scal_m,marker="D",s=22,color="#d62728",zorder=3,label="SCAL m* 2.1")
ax.set_xticks(x); ax.set_xticklabels(n.well); ax.set_ylabel("net pay (m)"); ax.grid(alpha=.3,axis="y"); ax.legend(loc="upper left")
plt.tight_layout(); plt.savefig(f"{O}/P1_Fig9_net_pay.png"); plt.close()
