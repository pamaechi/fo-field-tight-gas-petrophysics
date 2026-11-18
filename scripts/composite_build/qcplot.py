import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#222222", "#6b6b66", "#e4e3de"

def qc_plot(comp, out, title, top, base, markers=None, bands=None, cores=None):
    s = comp.loc[top:base]; d = s.index.values
    tracks = [
        ("GR (API)", [("GR", C1)], (0, 150), False),
        ("Caliper (in)", [("CALI", C1), ("CALS", C2)], (6, 16), False),
        ("Resistivity (ohm·m)", [("LLD", C1), ("LLS", C2), ("ILD", C3)], (1, 1000), True),
        ("RHOB (g/cm³)", [("RHOB", C1), ("RHOB_ALT", C2)], (2.0, 2.9), False),
        ("NPHI (v/v, LS)", [("NPHI", C1)], (0.30, -0.05), False),
        ("PEF (b/e)", [("PEF", C1)], (0, 10), False),
        ("DT (µs/ft)", [("DT", C1)], (110, 50), False),
    ]
    fig, axes = plt.subplots(1, len(tracks), figsize=(2.05 * len(tracks), 11), sharey=True)
    for ax, (ttl, curves, lim, logx) in zip(axes, tracks):
        for name, col in curves:
            if name in s and s[name].notna().any():
                ax.plot(s[name], d, color=col, lw=1.0, label=name)
        if logx: ax.set_xscale("log")
        ax.set_xlim(*lim)
        ax.set_title(ttl, fontsize=8.5, color=INK, pad=4)
        ax.tick_params(labelsize=7, colors=MUTED, length=2)
        ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
        for sp in ax.spines.values(): sp.set_color(GRID)
        if len([c for c, _ in curves if c in s]) > 1:
            ax.legend(fontsize=6.5, loc="lower right", frameon=False, labelcolor=INK)
        for t, b, lab, col in (bands or []):
            ax.axhspan(t, b, color=col, alpha=0.12, lw=0)
        for t, b in (cores or []):
            ax.axhspan(t, b, xmin=0.94, xmax=1.0, color="#555555", alpha=0.6, lw=0)
        for z, lab in (markers or {}).items():
            ax.axhline(z, color=MUTED, lw=0.6, ls="--")
    for z, lab in (markers or {}).items():
        axes[0].text(3, z - 0.4, lab, fontsize=6.5, color=MUTED, va="bottom")
    for t, b, lab, col in (bands or []):
        axes[-1].text(108, (t + b) / 2, lab, fontsize=6.5, color=INK, va="center")
    axes[0].set_ylim(base, top); axes[0].set_ylabel("Logger depth (m below KB)", fontsize=8, color=INK)
    fig.suptitle(title, fontsize=10, color=INK, x=0.01, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(out, dpi=160); plt.close(fig)
