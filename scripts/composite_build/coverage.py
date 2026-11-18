import re, pandas as pd, numpy as np
f = pd.read_csv("work/inventory/all_frames.csv")
FAM = {
 "GR":        r"^(GR|SGR|CGR|ECGR|GRS|HGR|GR_EDTC|GRD)$",
 "Spectral GR (K/Th/U)": r"^(POTA|THOR|URAN|HFK|HTHO|HURA|K|TH|U|NGS.*)$",
 "SP":        r"^(SP)$",
 "Caliper":   r"^(CALI|CALS|HCAL|C1|C2|CAL|HD|LCAL|DCAL)$",
 "Deep res":  r"^(LLD|ILD|IDPH|RT|HDRS|RLA5|RLA4|AT90|AF90|HLLD|M2R9|RD|ILDC|AO90|AT60)$",
 "Shallow/micro res": r"^(LLS|ILM|SFLU|SFL|MSFL|RXO|RXOZ|LL8|MLL|HMRS|RLA1|RLA2|AT10|AT20|AO10|AO20|HLLS|MCFL|RXO8)$",
 "Density":   r"^(RHOB|RHOZ|ZDEN|DEN|RHO8|HRHO|RHOM)$",
 "Neutron":   r"^(NPHI|TNPH|NPOR|CNC|NPHS|APLC|HNPO|NPRL|CN|TNPB)$",
 "PEF":       r"^(PEF|PEFZ|PE|PEF8|PEFL)$",
 "Sonic DTC": r"^(DT|DTCO|DTC|DTL|DTLN|DTLF|DT4P|AC|DT24|DTCO_FMP)$",
 "Shear DTS": r"^(DTSM|DTS|DTSH|DT4S|DTST)$",
 "NMR":       r"^(TCMR|CMRP|CMFF|BFV|KTIM|KSDR|T2LM|CMRP_3MS|MPHI|MBVI|MFFI|BVI|FFI|PHIE_CMR|T2_DIST)$",
 "Image/FMI": r"^(FCA.*|FMI.*|IMG.*|FBST.*|BS1.*|PAD.*)$",
 "Pressure (RFT/XPT/MDT)": r"^(SGP|HPGP|HPT|PRES|QPF|BQP)$",
}
rows = []
for w, g in f.groupby("well"):
    for fam, rx in FAM.items():
        hits = []
        for _, r in g.iterrows():
            names = str(r.curves).split()
            m = [n for n in names if re.match(rx, n.split(".")[0], re.I)]
            if m: hits.append((r.file, r.top_m, r.base_m, set(m), r.index))
        if not hits:
            rows.append(dict(well=w, family=fam, files=0, mnems="", top=np.nan, base=np.nan)); continue
        tops = [h[1] for h in hits if np.isfinite(h[1]) and h[1] < 6000]; bases = [h[2] for h in hits if np.isfinite(h[2]) and h[2] < 6000]
        rows.append(dict(well=w, family=fam, files=len({h[0] for h in hits}), mnems=" ".join(sorted(set().union(*[h[3] for h in hits])))[:60],
                         top=min(tops) if tops else np.nan, base=max(bases) if bases else np.nan))
cov = pd.DataFrame(rows)
cov.to_csv("work/inventory/coverage.csv", index=False)
piv = cov.assign(cell=lambda d: np.where(d.files > 0, d.apply(lambda r: f"{r.mnems[:22]} [{r.top:.0f}-{r.base:.0f}]" if np.isfinite(r.top) else r.mnems[:30], axis=1), "—")).pivot(index="family", columns="well", values="cell")
pd.set_option("display.width", 400); pd.set_option("display.max_colwidth", 40)
print(piv.reindex(list(FAM)).to_string())
