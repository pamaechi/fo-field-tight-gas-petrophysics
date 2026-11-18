"""Machine-learning permeability with blind (leave-one-well-out) validation (Notebook 07, Paper 3).

Every score is computed on wells the model has not seen. Random K-fold splits are provided only to show how much they
over-state skill (neighbouring plugs from the same well end up on both sides of the split).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier, HistGradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from . import rocktype as rt

ROCK = ["PHIT", "VSH_GR", "V_CLAY", "V_CARB", "RHOMAA", "NPHI", "DSEP"]      # rock / pore-system features (no resistivity)
RF_REG = dict(n_estimators=400, min_samples_leaf=5, max_features=0.5, random_state=0, n_jobs=-1)
RF_CLF = dict(n_estimators=300, min_samples_leaf=5, class_weight="balanced", random_state=0, n_jobs=-1)
GBM = dict(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=15, l2_regularization=1.0, random_state=0)


def regressors() -> dict:
    return {"MLR": make_pipeline(StandardScaler(), LinearRegression()),
            "RF": RandomForestRegressor(**RF_REG),
            "GBM": HistGradientBoostingRegressor(**GBM)}


def sample_logs(core: pd.DataFrame, logs: dict, feats, extra=("BADHOLE",)) -> pd.DataFrame:
    """3-sample (~0.45 m) log average at each core plug's log depth."""
    rows = []
    for w, g in core.groupby("well"):
        d = logs[w]; z = d.index.to_numpy()
        for _, r in g.iterrows():
            i = int(np.clip(np.searchsorted(z, r.depth_log_m), 1, len(z) - 2)); seg = d.iloc[i - 1:i + 2]
            rows.append({**r.to_dict(), **{f: seg[f].mean() for f in feats}, **{e: seg[e].max() for e in extra if e in seg}})
    return pd.DataFrame(rows)


def folds(T: pd.DataFrame, how: str = "well", k: int = 5):
    if how == "well":
        return [(T.index[T.well != w], T.index[T.well == w]) for w in T.well.unique()]
    return [(T.index[a], T.index[b]) for a, b in KFold(k, shuffle=True, random_state=0).split(T)]


def cross_predict(T: pd.DataFrame, F, how: str = "well") -> pd.DataFrame:
    """Out-of-fold log10 k for: field transform, rock-type transforms (Notebook 06) and the ML regressors."""
    P = pd.DataFrame(index=T.index)
    for tr, te in folds(T, how):
        a, b = T.loc[tr], T.loc[te]
        P.loc[te, "FIELD"] = np.polyval(np.polyfit(100 * a.PHIT, a.lk, 1), 100 * b.PHIT)
        clf = RandomForestClassifier(**RF_CLF).fit(a[ROCK], a.RT)
        P.loc[te, "ROCKTYPE"] = np.log10(rt.perm_rt(b.PHIT.to_numpy(), clf.predict_proba(b[ROCK]), clf.classes_,
                                                    rt.rt_transforms(100 * a.PHIT, a.lk, a.RT)))
        for n, m in regressors().items():
            P.loc[te, n] = m.fit(a[F], a.lk).predict(b[F])
    return P


def score(P: pd.DataFrame, y: pd.Series, well: pd.Series) -> pd.DataFrame:
    out = {}
    for c in P:
        e = P[c] - y; per = e.groupby(well).apply(lambda x: np.sqrt((x ** 2).mean()))
        out[c] = dict(pooled_rmse_dex=np.sqrt((e ** 2).mean()), mean_well_rmse_dex=per.mean(), worst_well_rmse_dex=per.max(),
                      F_O1_bias_dex=e[well == "F-O1"].mean(), r=np.corrcoef(P[c], y)[0, 1])
    return pd.DataFrame(out).T.round(3)


def well_conformal(T: pd.DataFrame, F, q=(0.1, 0.9), n_estimators: int = 200):
    """Split-conformal interval by well: residual quantiles from leave-one-well-out predictions within the training wells.
    Returns (lo, hi) offsets in dex to subtract from a prediction: P10 = p - hi, P90 = p - lo (residual = pred - core)."""
    res = []
    for v in T.well.unique():
        a, b = T[T.well != v], T[T.well == v]
        m = RandomForestRegressor(**{**RF_REG, "n_estimators": n_estimators}).fit(a[F], a.lk)
        res.append(m.predict(b[F]) - b.lk.to_numpy())
    return tuple(np.quantile(np.concatenate(res), q))
