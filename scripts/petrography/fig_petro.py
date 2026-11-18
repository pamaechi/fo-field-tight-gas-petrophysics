import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":9,"axes.titlesize":9,"axes.labelsize":9,"legend.fontsize":7.5,"savefig.dpi":300})
O="outputs/paperfigs"
RTC={1:"#1f77b4",2:"#2ca02c",3:"#ff7f0e",4:"#d62728"}
a=pd.read_csv("fo1_ts_joined.csv").dropna(subset=["k"]); a["k"]=a.k.clip(lower=0.001)
b=pd.read_csv("terratek_ts_rt.csv").sort_values("k_plug")
fig=plt.figure(figsize=(7.2,6.0)); gs=fig.add_gridspec(2,2,height_ratios=[1.15,1],width_ratios=[1.25,1]); ax=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])]
cd=a.por_type.str.contains("cement-dissolution")
for rt in (1,2,3,4):
    s=a[a.RT==rt]
    ax[0].scatter(s.k,s.depth_core,c=RTC[rt],s=22,edgecolor="k",lw=.3,label=f"RT{rt}",zorder=3)
ax[0].scatter(a.k[cd],a.depth_core[cd],s=70,facecolor="none",edgecolor="k",lw=.8,label="cement-dissolution\nporosity",zorder=2)
ax[0].axhline(3711,color="#666",ls="--",lw=.8); ax[0].text(0.0012,3710.6,"3711 m",fontsize=7,color="#444")
ax[0].set(xscale="log",xlabel="plug permeability (mD)",ylabel="core depth (m)",title="(a) F-O1 cores 1–2"); ax[0].invert_yaxis()
ax[0].legend(loc="upper left",bbox_to_anchor=(1.0,1.0),fontsize=6.5,frameon=False); ax[0].grid(alpha=.3)
jit=np.array([-.08,.08])
ax[1].set_position;ax[1].scatter(a.qtz_cmt+np.random.default_rng(1).uniform(-.08,.08,len(a)),a.k,c=[RTC[int(r)] for r in a.RT],s=22,edgecolor="k",lw=.3)
ax[1].set(yscale="log",xlabel="quartz-cement score",ylabel="plug permeability (mD)",title="(b) F-O1 quartz cement",xticks=[1,1.5,2,3])
ax[1].grid(alpha=.3)
x=np.arange(len(b)); cols=[("clay_matrix","clay and matrix","#8c564b"),("sec_cement","authigenic cement","#7f7f7f"),("modal_por","visible porosity","#17becf")]
bot=np.zeros(len(b))
for c,l,col in cols: ax[2].bar(x,b[c],.7,bottom=bot,color=col,label=l); bot+=b[c].values
for i,(rt,k,w) in enumerate(zip(b.RT,b.k_plug,b.well)): ax[2].text(i,bot[i]+0.8,f"RT{rt}",ha="center",fontsize=6.5,color=RTC[rt],fontweight="bold")
ax[2].set_xticks(x); ax[2].set_xticklabels([f"{w}\n{k:.2g}" for w,k in zip(b.well,b.k_plug)],fontsize=7)
ax[2].set(ylabel="modal abundance (%)",title="(c) TerraTek point counts (F-O2, F-O3, F-O4), ordered by plug permeability",ylim=(0,58),xlabel="well and plug permeability (mD)"); ax[2].legend(loc="upper center",fontsize=7,ncol=3,frameon=False); ax[2].grid(alpha=.3,axis="y")
plt.tight_layout(); plt.savefig(f"{O}/P1_Fig11_petrography.png"); plt.close()
print(cd.sum(), (a.RT==1).sum(), ((a.RT==1)&cd).sum(), ((a.RT>1)&cd).sum(), (a.RT>1).sum())
