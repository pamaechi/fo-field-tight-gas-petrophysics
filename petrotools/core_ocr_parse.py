"""Parse OCR'd SOEKOR conventional core analysis tables (two layouts) into tidy CSV."""
import re, numpy as np, pandas as pd
DEP = re.compile(r"(3[3-9]\d{2})[.,:](\d{2})")
def num(tok, kind):
    t = tok.strip().replace("o", "0").replace("O", "0").replace(",", ".").replace(":", ".").replace("—", "-").replace("°", "")
    t = re.sub(r"[^0-9.\-]", "", t)
    if t in ("", "-", ".", "-."): return np.nan
    if kind == "perm":                       # leading "." often read as "-" or dropped (".017" -> "-017"/"017")
        t = t.lstrip("-")
        if "." not in t: t = "0." + t if len(t) == 3 else t
        if t.startswith("."): t = "0" + t
    if kind == "rhog":
        t = t.replace("-", "")
        if "." not in t and len(t) == 3: t = t[0] + "." + t[1:]
        if re.fullmatch(r"2\d6\d", t): t = "2." + t[2:]   # "2166" -> "2.66"? keep 2.6x
    if kind == "phi":
        t = t.lstrip("-")
    try: return float(t)
    except ValueError: return np.nan
def parse_layout_full(txt, well):
    """F-O2 style: SAMPLE | DEPTH | BOYLE [SUM] | KAIR KKL | [G O W] | [CAL DOL] | RHOG"""
    rows = []; core = None
    for ln in txt.splitlines():
        mc = re.search(r"CORE\s*([1-9])\b", ln, re.I)
        if mc and not DEP.search(ln): core = int(mc.group(1))
        m = DEP.search(ln)
        if not m: continue
        dep = float(m.group(1) + "." + m.group(2)); rest = ln[m.end():]
        rest = re.sub(r"(\d)\s*[oO]\s*(\d)", r"\1 0 \2", rest)
        rest = re.sub(r"\b2[1\-]([5-9]\d)\b", r"2.\1", rest)          # 2166 / 2-68 -> 2.66 / 2.68
        rest = re.sub(r"(?<=\d)-(?=\d\d\b)", ".", rest).replace(":", ".")
        decs = re.findall(r"\d+[.,:]\d+", rest.replace("o", "0").replace("O", "0"))
        ints = re.findall(r"(?<![\d.])\d{1,3}(?![\d.])", re.sub(r"\d+[.,:]\d+", " ", rest))
        vals = [float(x.replace(",", ".").replace(":", ".")) for x in decs]
        rh = re.search(r"2[.,:1]?6\d|2[.,]\d\d\s*\W*$", rest)
        r = dict(well=well, core=core, depth_m=dep, raw=ln.strip())
        if len(vals) == 3:
            r["phi_pct"], r["kair_md"], r["kkl_md"] = vals; r["flag"] = "no grain density read"
        if len(vals) in (4, 7):
            r["phi_pct"] = vals[0]
            k0 = 2 if len(vals) == 7 else 1
            if len(vals) == 7: r["phi_sum_pct"] = vals[1]
            r["kair_md"], r["kkl_md"] = vals[k0], vals[k0 + 1]; r["rhog"] = vals[-1]
        rows.append(r)
    return pd.DataFrame(rows)

MANUAL = {  # rows re-read by eye from the scanned page where OCR failed: (well, depth) -> values
    ("F-O2", 3616.38): dict(phi_pct=9.0, kair_md=0.13, kkl_md=0.05, rhog=2.66),
    ("F-O2", 3635.11): dict(phi_pct=6.5, phi_sum_pct=8.8, kair_md=0.08, kkl_md=0.03, rhog=2.64),
    ("F-O2", 3649.06): dict(phi_pct=13.7, phi_sum_pct=14.0, kair_md=9.93, kkl_md=8.77),
    ("F-O2", 3729.13): dict(phi_pct=9.8, phi_sum_pct=11.7, kair_md=0.11, kkl_md=0.05),
}
def apply_manual(d):
    for (w, z), v in MANUAL.items():
        m = (d.well == w) & (d.depth_m.round(2) == z)
        for k, x in v.items(): d.loc[m, k] = x
        d.loc[m, "flag"] = "manually re-read from scan"
    return d

TOK = re.compile(r"[-—+.]?\d+(?:[.,]\d+)?")
def _is_phi(t): return re.fullmatch(r"\d{1,2}[.,]\d", t) is not None and 0.5 <= float(t.replace(",", ".")) <= 35
def _is_perm(t):
    s = t.replace("—", "-").replace(",", ".")
    return (re.fullmatch(r"[-+.]?\d{3}", s) is not None) or (re.fullmatch(r"[-+.]?\d+\.\d{2,3}", s) is not None)
def _perm_val(t):
    s = t.replace("—", "-").replace(",", ".").lstrip("-+")
    if s.startswith("."): s = "0" + s
    if "." not in s: s = "0." + s
    return float(s)
def parse_layout_cms(txt, well):
    """F-O3 / F-S1 style (SOEKOR 1990s): CORE | ID | DEPTH | phi(summation, gas exp. 800 psi, 4121 psi, 0 psi) |
    K CMS-200 (KL 800, Ka 800, KL 4121, Ka 4121) | minipermeameter (Ka, Kl) | Sg So Sw | Cal Dol | rhog"""
    rows = []
    for ln in txt.splitlines():
        m = DEP.search(ln)
        if not m: continue
        mc = re.match(r"\W*([1-9])\W", ln); core = int(mc.group(1)) if mc else None
        dep = float(m.group(1) + "." + m.group(2)); rest = ln[m.end():].replace("o", "0").replace("O", "0")
        toks = TOK.findall(rest)
        rh = [t for t in toks if re.fullmatch(r"2[.,][5-9]\d", t)]
        rhog = float(rh[-1].replace(",", ".")) if rh else np.nan
        phis, perms, i = [], [], 0
        while i < len(toks) and _is_phi(toks[i].lstrip("-+.")) and not _is_perm(toks[i]): phis.append(float(toks[i].lstrip("-+.").replace(",", "."))); i += 1
        while i < len(toks) and _is_perm(toks[i]) and not re.fullmatch(r"2[.,][5-9]\d", toks[i]): perms.append(_perm_val(toks[i])); i += 1
        r = dict(well=well, core=core, depth_m=dep, raw=ln.strip(), rhog=rhog)
        if phis:
            r["phi_pct"] = round(float(np.mean(phis)), 2); r["phi_values"] = " ".join(map(str, phis))
        if perms:
            r["kkl_md"] = perms[0]; r["kair_md"] = perms[1] if len(perms) > 1 else np.nan
        rows.append(r)
    return pd.DataFrame(rows)

def parse_fo1_appJ(txt, well="F-O1"):
    """F-O1 Well Completion Appendix J: depth | gas-exp phi | [sum-fluids phi] | air perm | liquid perm | [Sw So Sg Cal Dol] | [rhog] | overburden psi.
    Leading '0.' of permeabilities is lost in the scan ('93' = 0.93)."""
    rows = []
    for ln in txt.splitlines():
        s = ln.replace(",", ".").replace(" .", ".")
        m = re.search(r"\b(3[17]\d{2})\s?\.\s?(\d{2})\b", s)
        if not m: continue
        dep = float(("37" + m.group(1)[2:]) + "." + m.group(2))       # '3107.07' misread -> 3707.07
        ob = re.search(r"(800|6100)\s*$", s.strip())
        if not ob: continue
        rest = s[m.end():ob.start()].replace("l", "1").replace("I", "1").replace("O", "0").replace("o", "0")
        toks = re.findall(r"\d+(?:\.\d+)?", rest)
        if not toks: continue
        rhog = np.nan
        if re.fullmatch(r"2\.[5-9]\d", toks[-1]): rhog = float(toks.pop())
        phi = float(toks[0]) if re.fullmatch(r"\d{1,2}\.\d", toks[0]) else np.nan
        i = 1
        if i < len(toks) and re.fullmatch(r"\d{1,2}\.\d", toks[i]) and len(toks) >= 5: i += 1   # sum-of-fluids porosity
        def k(t):
            if "." in t: return float(t)
            return float("0." + t) if len(t) <= 2 else np.nan
        ka = k(toks[i]) if i < len(toks) else np.nan
        kl = k(toks[i + 1]) if i + 1 < len(toks) else np.nan
        rows.append(dict(well=well, depth_m=dep, phi_pct=phi, kair_md=ka, kkl_md=kl, rhog=rhog, overburden_psi=int(ob.group(1)), raw=ln.strip()))
    return pd.DataFrame(rows)

def _dec(t, kind):
    """Restore lost decimal points in F-O4 Corex Table 5.1 by column type."""
    t = t.strip().replace(",", ".").replace("O", "0").replace("o", "0").replace("I", "1").replace("l", "1").replace("S", "5")
    t = re.sub(r"[^0-9.]", "", t).strip(".")
    if t.count(".") > 1: t = t[:t.index(".") + 1] + t[t.index(".") + 1:].replace(".", "")
    if not t: return np.nan
    if "." in t:
        v = float(t)
    elif kind == "k":
        v = float("0." + t[1:]) if t.startswith("0") and len(t) > 1 else (float(t[0] + "." + t[1:]) if len(t) >= 3 else float(t))
    elif kind == "phi":
        v = float(t[:-1] + "." + t[-1]) if len(t) >= 2 else float(t)
    elif kind == "rhog":
        v = float(t[0] + "." + t[1:]) if len(t) == 3 else np.nan
    elif kind == "sw":
        v = float(t[:-1] + "." + t[-1]) if len(t) == 3 else float(t)
    else: v = float(t)
    ok = {"k": (0, 5000), "phi": (0.5, 30), "rhog": (2.5, 3.0), "sw": (0, 100)}[kind]
    return v if ok[0] <= v <= ok[1] else np.nan
def parse_fo4_table51(pages_txt, well="F-O4"):
    rows = []
    for ln in pages_txt.splitlines():
        ms = list(DEP.finditer(ln))
        if not ms: continue
        m = ms[0]; dep = float(m.group(1) + "." + m.group(2))
        seg = ln[m.end(): ms[1].start()] if len(ms) > 1 else ln[m.end():]
        seg = re.sub(r"G[DE]\w*\s*$", "", seg)                       # drop the next sample id (vertical plug)
        toks = [t for t in re.split(r"[\s|\]\[}{)(;:]+", seg) if re.search(r"\d", t) or "rac" in t or re.fullmatch(r"[IlS0-9.]{2,4}", t)]
        if len(toks) < 3: continue
        frac = "rac" in toks[0]
        k = np.nan if frac else _dec(toks[0], "k")
        phi = _dec(toks[1], "phi"); rh = _dec(toks[2], "rhog"); sw = _dec(toks[3], "sw") if len(toks) > 3 else np.nan
        flag = "fractured plug" if frac else ""
        if phi > 20: phi = np.nan; flag += " phi>20 (OCR column shift) nulled"
        if np.isnan(rh) and re.fullmatch(r"2[.,][5-9]\d", toks[1].strip()):        # grain density slid into the porosity column
            phi = np.nan; flag += " porosity unreadable (column shift) nulled"
            if not re.search(r"\d", toks[0].replace(".", "")) or re.fullmatch(r"\d+\.\d", toks[0]): k = np.nan
        if re.fullmatch(r"1\d", re.sub(r"[^0-9.]", "", toks[1])):                  # '12' = 1.2 or 12.x: ambiguous ('87' = 8.7 is not)
            phi = np.nan; flag += " porosity ambiguous (decimal lost) nulled"
        rows.append(dict(well=well, depth_m=dep, kkl_md=k, phi_pct=phi, rhog=rh, sw_pct=sw, flag=flag.strip(), raw=ln.strip()))
    return pd.DataFrame(rows)
