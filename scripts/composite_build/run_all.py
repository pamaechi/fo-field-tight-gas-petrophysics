import sys, traceback, pandas as pd
from splice_all import splice, write, sanity
from wells_config import WELLS
import clean as CL
from clean import clean
for w in (sys.argv[1:] or WELLS):
    c = WELLS[w]
    try:
        CL.PASS_LOG.clear()
        comp, units, prov = splice(w, c["segs"], c["top"], c["base"], c["notes"])
        for cv, t0, b0 in c.get("null", []):
            if cv in comp: comp.loc[t0:b0, cv] = float("nan")
        comp, clog = clean(w, comp)
        clog = pd.concat([pd.DataFrame(CL.PASS_LOG).assign(well=w, stage="per-pass"), clog.assign(stage="post-splice")], ignore_index=True)
        clog.to_csv(f"work/out/composites/{w}_cleaning_log.csv", index=False)
        print("  cleaned:", len(clog), "intervals,", int(clog.samples.sum()) if len(clog) else 0, "samples")
        write(w, comp, units, prov, c["notes"], f"work/out/composites")
        sh = prov.drop_duplicates(["source_file","frame"])[["source_file","frame","role","depth_shift_m","gr_lag_check"]]
        print(f"== {w}: {comp.shape[0]} samples {comp.index.min()}-{comp.index.max()} m, curves: {' '.join(comp.columns)}")
        print(sh.to_string(index=False)); print("  sanity (curve, median, % outside range):", sanity(comp))
    except Exception as e:
        print("FAILED", w, e); traceback.print_exc()
