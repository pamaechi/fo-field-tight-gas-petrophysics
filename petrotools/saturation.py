"""Water saturation tools for the F-O field (Notebook 04): temperature, Rw(T), Archie, Waxman-Smits with Qv from CEC."""
from __future__ import annotations
import numpy as np
import pandas as pd

# deep resistivity used as Rt in the reservoir (laterolog where available; HRLA in F-O6/F-O8; induction only in F-R1)
RT_CURVE = {"F-O1": "LLD", "F-O2": "LLD", "F-O3": "LLD", "F-O4": "LLD", "F-O6": "RLA5", "F-O8": "RLA5", "F-R1": "ILD", "F-S1": "LLD"}

def temperature_F(md_m, t_ref_F: float = 309.7, z_ref_m: float = 3644.0, t_surf_F: float = 60.0):
    """Linear geotherm through the F-O2 DST temperature (309.7 degF at 3644 m)."""
    return t_surf_F + (t_ref_F - t_surf_F) * np.asarray(md_m, float) / z_ref_m

def rw_at_T(rw_ref: float, t_ref_F: float, t_F):
    """Arps temperature correction."""
    return rw_ref * (t_ref_F + 6.77) / (np.asarray(t_F, float) + 6.77)

def archie_sw(rt, phi, rw, a: float = 1.0, m: float = 2.0, n: float = 2.0):
    return np.clip(((a * rw) / (np.power(phi, m) * rt)) ** (1.0 / n), 0, 1)

def qv_from_cec(cec_meq100g, phi, rhog: float = 2.67):
    """Qv (meq/ml pore volume) = CEC (1 - phi) rho_g / (100 phi)."""
    return cec_meq100g * (1 - phi) * rhog / (100 * phi)

def b_juhasz(t_F, rw):
    """Equivalent counter-ion conductance B (S cm2/meq, i.e. (1/ohm.m)/(meq/ml)) - Juhasz (1981)."""
    tc = (np.asarray(t_F, float) - 32) / 1.8
    return (-1.28 + 0.225 * tc - 0.0004059 * tc ** 2) / (1 + np.power(rw, 1.23) * (0.045 * tc - 0.27))

def waxman_smits_sw(rt, phi, rw, qv, b, m_star: float = 2.0, n_star: float = 2.0, iters: int = 60):
    """Solve 1/Rt = phi^m* Sw^n* (1/Rw + B Qv / Sw) for Sw by fixed-point iteration."""
    rt, phi, rw, qv, b = (np.asarray(x, float) for x in (rt, phi, rw, qv, b))
    sw = archie_sw(rt, phi, rw, 1.0, m_star, n_star)
    for _ in range(iters):
        sw = np.clip(sw, 1e-3, 1)
        sw = np.clip((1.0 / (rt * np.power(phi, m_star) * (1.0 / rw + b * qv / sw))) ** (1.0 / n_star), 1e-3, 1)
    return sw

def calibrate_m_star(rt, phi, rw, qv, b, n_star: float = 2.4, target: float = 1.0, lo: float = 1.6, hi: float = 3.2):
    """m* such that the MEDIAN Waxman-Smits Sw of water-leg samples equals `target` (Sw is not used above 1, so the
    median of the unclipped solution is matched via bisection on the clipped median approaching 1 from below)."""
    rt, phi, rw, qv, b = (np.asarray(x, float) for x in (rt, phi, rw, qv, b))
    def med_sw(m):  # unclipped: solve WS for Sw allowing >1 by using a large clip
        sw = np.full_like(rt, 0.5)
        for _ in range(80):
            sw = np.clip((1.0 / (rt * np.power(phi, m) * (1.0 / rw + b * qv / np.clip(sw, 1e-3, None)))) ** (1.0 / n_star), 1e-3, 5)
        return float(np.median(sw))
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if med_sw(mid) < target: lo = mid
        else: hi = mid
    return round(0.5 * (lo + hi), 3)
