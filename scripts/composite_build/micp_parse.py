import re, numpy as np, pandas as pd
def parse(path):
    t=open(path).read().splitlines(); samples=[]; cur=None
    for ln in t:
        m=re.search(r"Sample N\s?umber:\s*(\S+)",ln)
        if m:
            if cur and cur['rows']: samples.append(cur)
            cur=dict(sample=m.group(1),rows=[],seen_max=False); continue
        if cur is None: continue
        for key,pat in [('depth_m',r"Sample D\s?epth[^:]*:\s*([\d.]+)"),('k_md',r"Permeability to Air[^:]*:\s*([\d.]+)"),('phi',r"Porosity, fraction:\s*([\d.]+)"),('rhog',r"Grain Density[^:]*:\s*([\d.]+)"),('r50_um',r"Median Pore Throat R\s?adius[^:]*:\s*([\d.]+)")]:
            mm=re.search(pat,ln)
            if mm: cur[key]=float(mm.group(1))
        nums=re.findall(r"-?\d+\.?\d*",ln)
        if len(nums)>=6 and re.match(r"^\s*[\d.]+\s+[\d.]+\s+[\d.]+",ln):
            pc,shg,sw=float(nums[0]),float(nums[1]),float(nums[2])
            if 0<=shg<=1.0001 and abs(shg+sw-1)<0.01: cur['rows'].append((pc,shg))
    if cur and cur['rows']: samples.append(cur)
    out=[]
    for s in samples:
        r=np.array(s['rows']); 
        # drainage branch: up to max pressure (first occurrence)
        imax=int(np.argmax(r[:,0])); dr=r[:imax+1]
        for pc,shg in dr: out.append(dict(sample=s['sample'],depth_m=s.get('depth_m'),k_md=s.get('k_md'),phi=s.get('phi'),rhog=s.get('rhog'),r50_um=s.get('r50_um'),pc_hg_psia=pc,shg=shg))
    return pd.DataFrame(out)
if __name__=='__main__':
    d=parse('work/text/F-O6__Other__Well_F-O6_Special_Core_Analysis_Final_Report.pdf.txt')
    s=d.groupby('sample').agg(depth=('depth_m','first'),k=('k_md','first'),phi=('phi','first'),r50=('r50_um','first'),npts=('pc_hg_psia','size'),pmax=('pc_hg_psia','max'),shg_max=('shg','max'))
    print(s.to_string()); d.to_csv('work/micp_F-O6.csv',index=False)
