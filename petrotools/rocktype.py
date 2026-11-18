"""Rock typing (Notebook 06): Winland r35, FZI, rock-type classes and probability-weighted rock-type permeability."""
from __future__ import annotations
import numpy as np
import pandas as pd

RT_EDGES = [0, 0.1, 0.5, 2.0, 1e9]          # r35 (micron): RT4 nano | RT3 micro | RT2 meso | RT1 macro/mega
RT_LABELS = [4, 3, 2, 1]
RT_NAMES = {1: "RT1 macro (r35 >= 2 um) - enhanced porosity", 2: "RT2 meso (0.5-2 um)", 3: "RT3 micro (0.1-0.5 um)", 4: "RT4 nano (< 0.1 um)"}

def winland_r35(k_md, phi_pct):
    return 10 ** (0.732 + 0.588 * np.log10(k_md) - 0.864 * np.log10(phi_pct))

def fzi(k_md, phi_frac):
    rqi = 0.0314 * np.sqrt(np.asarray(k_md) / np.asarray(phi_frac)); return rqi / (np.asarray(phi_frac) / (1 - np.asarray(phi_frac)))

def rock_type(r35):
    return pd.cut(pd.Series(r35), RT_EDGES, labels=RT_LABELS).astype(float)

def rt_transforms(phit_pct, logk, rt, min_n: int = 10):
    df = pd.DataFrame(dict(x=phit_pct, y=logk, rt=rt)).dropna()
    return {int(r): np.polyfit(g.x, g.y, 1) for r, g in df.groupby("rt") if len(g) >= min_n}

def perm_rt(phit_frac, proba: np.ndarray, classes, tf: dict):
    """Probability-weighted log10 k over rock types (classes with a transform only)."""
    num = np.zeros(len(phit_frac)); den = np.zeros(len(phit_frac))
    for j, c in enumerate(classes):
        if int(c) in tf: num += proba[:, j] * np.polyval(tf[int(c)], 100 * np.asarray(phit_frac)); den += proba[:, j]
    return 10 ** (num / np.where(den > 0, den, np.nan))
