"""XRD-calibrated linear mineral solver (Notebook 06).

Model: each log  L = sum_i V_i * R_i(L) + PHI * R_f(L),  with  sum_i V_i + PHI = 1,  V_i, PHI >= 0.
Mineral groups: QF (quartz + feldspar), CLAY (illite/mica/glauconite + chlorite + kaolinite), CARB (ankerite/Fe-dolomite +
siderite + calcite + dolomite). Fluid responses fixed (invaded zone): RHOB 1.0 g/cc, NPHI 1.0, DT 189 us/ft, VSH_GR 0.
Mineral responses R_i are calibrated to XRD by bounded least squares, not taken from textbook tables.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear

GROUPS = ["QF", "CLAY", "CARB"]
GRAIN_RHO = {"QF": 2.64, "CLAY": 2.80, "CARB": 3.05}           # for wt% -> vol% conversion
FLUID = {"RHOB_USE": 1.0, "NPHI": 1.0, "DT": 189.0, "VSH_GR": 0.0}
TEXTBOOK = {"RHOB_USE": (2.65, 2.75, 3.00), "NPHI": (-0.02, 0.30, 0.02), "DT": (55.5, 85.0, 45.0), "VSH_GR": (0.0, 1.0, 0.0)}
# QF and CARB held near physical values (sparse calibration set); CLAY (illite/mica/glauconite mixture + bound water) free
BOUNDS = {"RHOB_USE": ((2.63, 2.55, 2.85), (2.68, 3.00, 3.20)), "NPHI": ((-0.04, 0.10, -0.02), (0.00, 0.60, 0.06)),
          "DT": ((53, 60, 43), (57, 130, 50)), "VSH_GR": ((0, 0.3, 0), (0.05, 3.0, 0.15))}

def xrd_to_volumes(xrd: pd.DataFrame, phi: pd.Series) -> pd.DataFrame:
    """wt% groups -> bulk volume fractions using grain densities and porosity."""
    w = pd.DataFrame({"QF": xrd.qtz + xrd.kfs + xrd.plag, "CLAY": xrd.clay, "CARB": xrd.cal + xrd.sid + xrd.ank + xrd.dol})
    v = pd.DataFrame({g: w[g] / GRAIN_RHO[g] for g in GROUPS}); v = v.div(v.sum(axis=1), axis=0).mul(1 - phi.values, axis=0)
    v["PHI"] = phi.values
    return v

def calibrate(vol: pd.DataFrame, logs: pd.DataFrame, tools=("RHOB_USE", "NPHI", "DT", "VSH_GR")) -> dict:
    """Mineral responses per tool: min ||V R - (L - PHI R_f)|| with physical bounds."""
    resp = {}
    for t in tools:
        m = logs[t].notna().to_numpy()
        A = vol.loc[m, GROUPS].to_numpy(); y = logs.loc[m, t].to_numpy() - vol.loc[m, "PHI"].to_numpy() * FLUID[t]
        r = lsq_linear(A, y, bounds=BOUNDS[t]); resp[t] = dict(zip(GROUPS, np.round(r.x, 4)))
    return resp

def invert(logs: pd.DataFrame, resp: dict, phi_prior: pd.Series | None = None, weights: dict | None = None,
           w_unity: float = 10.0, w_phi: float = 3.0) -> pd.DataFrame:
    """Per-sample bounded least squares for V_QF, V_CLAY, V_CARB, PHI. Tool uncertainties set the weights."""
    weights = weights or {"RHOB_USE": 1 / 0.03, "NPHI": 1 / 0.03, "DT": 1 / 4.0, "VSH_GR": 1 / 0.15}
    tools = [t for t in resp if t in logs]
    out = np.full((len(logs), 4), np.nan); resid = np.full(len(logs), np.nan)
    L = logs[tools].to_numpy(); P = phi_prior.to_numpy() if phi_prior is not None else None
    for i in range(len(logs)):
        ok = np.isfinite(L[i])
        if ok.sum() < 2: continue
        rows, y = [], []
        for t, o, v in zip(tools, ok, L[i]):
            if o: rows.append([resp[t][g] * weights[t] for g in GROUPS] + [FLUID[t] * weights[t]]); y.append(v * weights[t])
        rows.append([w_unity] * 4); y.append(w_unity)
        if P is not None and np.isfinite(P[i]): rows.append([0, 0, 0, w_phi / 0.02]); y.append(P[i] * w_phi / 0.02)
        r = lsq_linear(np.array(rows), np.array(y), bounds=(0, 1)); out[i] = r.x; resid[i] = np.sqrt(np.mean(r.fun ** 2))
    df = pd.DataFrame(out, index=logs.index, columns=["V_QF", "V_CLAY", "V_CARB", "PHI_MIN"])
    df["MIN_RESID"] = resid
    return df

def apparent_matrix_density(rhob, phi, rho_f: float = 1.0):
    return (np.asarray(rhob) - np.asarray(phi) * rho_f) / (1 - np.asarray(phi))

def fit_carbonate_index(rhomaa: pd.Series, carb_frac: pd.Series):
    """Linear calibration of carbonate volume fraction of solids against apparent matrix density (XRD sands)."""
    b = np.polyfit(rhomaa, carb_frac, 1); r = float(np.corrcoef(rhomaa, carb_frac)[0, 1])
    rmse = float(np.sqrt(np.mean((np.polyval(b, rhomaa) - carb_frac) ** 2)))
    return dict(a=float(b[1]), b=float(b[0]), r=r, rmse=rmse)

def mineral_model(d: pd.DataFrame, resp: dict, carb_cal: dict) -> pd.DataFrame:
    """Per-depth volumes: PHIT (porosity of record), V_CLAY (solver), V_CARB (apparent matrix density from a density-free
    porosity), V_QF = remainder. Inputs: RHOB_USE, NPHI, DT (or DTCO), VSH_GR, PHIT."""
    x = d.copy()
    if "DT" not in x or x["DT"].notna().sum() < 50: x["DT"] = x.get("DTCO")
    full = invert(x[["RHOB_USE", "NPHI", "DT", "VSH_GR"]], resp, None)
    nd = invert(x[["NPHI", "DT", "VSH_GR"]], {k: v for k, v in resp.items() if k != "RHOB_USE"}, None)
    rhomaa = apparent_matrix_density(x.RHOB_USE, nd.PHI_MIN)
    solid = (1 - x.PHIT)
    v_clay = (full.V_CLAY / (1 - full.PHI_MIN)).clip(0, 1) * solid
    v_carb = (carb_cal["a"] + carb_cal["b"] * rhomaa).clip(0, 0.6) * solid
    v_qf = (solid - v_clay - v_carb).clip(lower=0)
    return pd.DataFrame({"V_QF": v_qf, "V_CLAY": v_clay, "V_CARB": v_carb, "RHOMAA": rhomaa, "PHI_ND": nd.PHI_MIN,
                         "MIN_RESID": full.MIN_RESID}, index=d.index)
