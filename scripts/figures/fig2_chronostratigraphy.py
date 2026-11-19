"""Chronostratigraphic chart of the Bredasdorp Basin, redrawn and simplified
(modified after Petroleum Agency SA, 2004/05; Brown et al., 1995; McMillan et al., 1997)."""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
import numpy as np

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5})

# stages, youngest at top; row index 0 = top
stages = ["PL", "MIO", "OL", "EO", "PA", "MA", "CA", "SA", "CO", "TU", "CE", "AL", "AP",
          "BA", "HA", "VA", "BE", "PO", "KIM", "OXF"]
N = len(stages)
row = {s: i for i, s in enumerate(stages)}
def y(s, frac=0.0):            # top of stage s = y(s); frac 1 = base
    return N - row[s] - frac

# unconformities: name, stage whose BASE it sits at, age (Ma, as in source chart)
unconf = [("22At1", "PA", 67), ("17At1", "CA", 80), ("15At1", "TU", 93), ("14At1", "AL", 103),
          ("13At1", "AP", 112), ("6At1", "BA", 117.5), ("1At1", "HA", 126), ("D", "OXF", None)]

# column x-positions (left edge, width)
cols = {"era": (0, 0.55), "stage": (0.55, 0.55), "ma": (1.10, 0.45), "unc": (1.55, 0.85),
        "seq": (2.40, 4.2), "env": (6.60, 1.55), "tect": (8.15, 1.65), "res": (9.80, 0.65),
        "shows": (10.45, 0.65), "trap": (11.10, 1.55)}
TOP = N

C = dict(alluvial="#9ccc65", fluvial="#e57373", lacustrine="#ec407a", lagoon="#8e7cc3",
         deltaic="#f6a23a", shallow="#fff176", shelf="#43a047", slope="#80cbc4",
         basin="#4a8fd4", hiatus="#ffffff", basement="#9e9e9e")

fig, ax = plt.subplots(figsize=(12.0, 7.2), dpi=200)
ax.set_xlim(-0.05, 15.4); ax.set_ylim(-1.0, TOP + 1.55); ax.axis("off")

def box(x0, w, y0, y1, fc="white", ec="k", lw=0.6, **kw):
    ax.add_patch(Rectangle((x0, y0), w, y1 - y0, fc=fc, ec=ec, lw=lw, **kw))

# ---- headers
hdr = {"era": "AGE", "stage": "", "ma": "Ma", "unc": "MAJOR\nUNCON-\nFORMITIES",
       "seq": "SEQUENCE STRATIGRAPHIC FRAMEWORK\n(schematic)", "env": "DEPOSITIONAL\nENVIRONMENT",
       "tect": "TECTONISM", "res": "RESER-\nVOIR", "shows": "GAS/\nOIL\nSHOWS", "trap": "TRAP"}
for k, (x0, w) in cols.items():
    if k == "stage":
        continue
    ww = w + (cols["stage"][1] if k == "era" else 0)
    box(x0, ww, TOP, TOP + 1.45, fc="#eeeeee")
    ax.text(x0 + ww / 2, TOP + 0.72, hdr[k], ha="center", va="center", fontsize=6.2, fontweight="bold")
ax.text(cols["seq"][0] + 0.05, TOP - 0.18, "PROXIMAL", fontsize=6, va="top")
ax.text(cols["seq"][0] + cols["seq"][1] - 0.05, TOP - 0.18, "DISTAL", fontsize=6, va="top", ha="right")

# ---- era/period column
periods = [("TERTIARY", "PL", "PA", "#fdd835"), ("UPPER\nCRETACEOUS", "MA", "CE", "#aed581"),
           ("LOWER\nCRETACEOUS", "AL", "BE", "#8bc34a"), ("JURASSIC", "PO", "OXF", "#4fc3f7")]
x0, w = cols["era"]
for name, a, b, fc in periods:
    box(x0, w, y(b, 1), y(a), fc=fc)
    ax.text(x0 + w / 2, (y(b, 1) + y(a)) / 2, name, rotation=90, ha="center", va="center",
            fontsize=6.5, fontweight="bold")
x0, w = cols["stage"]
for s in stages:
    box(x0, w, y(s, 1), y(s), fc="white")
    ax.text(x0 + w / 2, y(s, 0.5), s, ha="center", va="center", fontsize=6.5)

# ---- Ma and unconformities
xm, wm = cols["ma"]; xu, wu = cols["unc"]
box(xm, wm, 0, TOP); box(xu, wu, 0, TOP)
for name, s, age in unconf:
    yy = y(s, 1) if name != "D" else y("OXF", 0.6)
    ax.plot([xu, xu + wu], [yy, yy], color="k", lw=0.6)
    ax.text(xu + wu - 0.05, yy + 0.08, name, ha="right", va="bottom", fontsize=6.5, fontweight="bold")
    if age:
        ax.text(xm + wm / 2, yy, f"{age:g}", ha="center", va="center", fontsize=6.5,
                bbox=dict(fc="white", ec="none", pad=0.3))
ax.text(xu + 0.12, (y("OXF", 0.6) + y("HA", 1)) / 2, "SYNRIFT", rotation=90, ha="center",
        va="center", fontsize=6.5)
ax.annotate("", xy=(xu + 0.12, y("HA", 1) - 0.05), xytext=(xu + 0.12, y("OXF", 0.6) + 0.05),
            arrowprops=dict(arrowstyle="<->", lw=0.5))
ax.text(xu + 0.12, (y("HA", 1) + TOP) / 2, "DRIFT", rotation=90, ha="center", va="center", fontsize=6.5)

# ---- sequence stratigraphic framework (schematic facies belts, proximal -> distal)
xs, ws = cols["seq"]
box(xs, ws, 0, TOP)
def belts(s0, s1, fracs, fcs, f0=0.0, f1=1.0):
    """fill stages s0 (top) .. s1 (bottom) with vertical facies belts; fracs = belt edges 0..1"""
    yt, yb = y(s0, f0), y(s1, f1)
    for i, fc in enumerate(fcs):
        xa, xb = xs + ws * fracs[i], xs + ws * fracs[i + 1]
        ax.add_patch(Rectangle((xa, yb), xb - xa, yt - yb, fc=fc, ec="none"))

# pre-Mesozoic basement and Devonian Cape Supergroup
ax.add_patch(Rectangle((xs, y("OXF", 1)), ws, 0.35, fc=C["basement"], ec="none"))
ax.text(xs + ws / 2, y("OXF", 1) + 0.17, "PRE-MESOZOIC BASEMENT (incl. Devonian Cape Supergroup)", ha="center", va="center",
        fontsize=4.6, color="white", fontweight="bold")
# synrift: fluvial / lacustrine half-graben fill (Oxfordian-Portlandian)
belts("PO", "OXF", [0, 0.35, 0.7, 1.0], [C["alluvial"], C["fluvial"], C["lacustrine"]], f1=0.75)
# Berriasian-Valanginian: shallow marine to continental (USM sand sheet below 1At1)
belts("VA", "BE", [0, 0.2, 0.42, 0.78, 1.0], [C["fluvial"], C["deltaic"], C["shallow"], C["shelf"]])
# Hauterivian-Albian: deep marine; slope wedge and basin-floor fans
for s, fr in [("HA", [0, 0.18, 0.40, 1.0]), ("BA", [0, 0.16, 0.36, 1.0]), ("AP", [0, 0.2, 0.42, 1.0]),
              ("AL", [0, 0.25, 0.47, 1.0])]:
    belts(s, s, fr, [C["shelf"], C["slope"], C["basin"]])
# turbidite fan lobes (deep marine sandstones)
for s, xc in [("BA", 0.62), ("AP", 0.70), ("AL", 0.80)]:
    yc = y(s, 0.55)
    ax.add_patch(Polygon([[xs + ws * (xc - 0.16), yc], [xs + ws * xc, yc + 0.22],
                          [xs + ws * (xc + 0.16), yc]], fc=C["shallow"], ec="k", lw=0.3))
# Upper Cretaceous: regressive cycles, prograding shelf and coastal plain
ucs = ["CE", "TU", "CO", "SA", "CA", "MA"]
for i, s in enumerate(ucs[::-1]):           # older at bottom: progradation toward distal
    shift = 0.08 * i
    belts(s, s, [0, 0.18 + shift, 0.30 + shift, 0.55 + shift, 0.72 + shift, 1.0],
          [C["alluvial"], C["deltaic"], C["shallow"], C["shelf"], C["slope"]])
# Tertiary: thin shelf succession with non-depositional hiatus
belts("PA", "PA", [0, 0.3, 0.7, 1.0], [C["deltaic"], C["shelf"], C["slope"]])
belts("PL", "EO", [0, 0.25, 0.6, 1.0], [C["hiatus"], C["shelf"], C["slope"]])
ax.add_patch(Rectangle((xs, y("EO", 1)), ws * 0.25, y("PL") - y("EO", 1), fc="white", ec="none",
                       hatch="////", lw=0))
for name, s, age in unconf:
    if name == "D":
        continue
    ax.plot([xs, xs + ws], [y(s, 1)] * 2, color="k", lw=0.9)
ax.plot([xs, xs + ws], [y("OXF", 0.6)] * 2, color="k", lw=0.9, ls="--")
for s in stages:
    ax.plot([xs, xs + ws], [y(s, 1)] * 2, color="k", lw=0.15, alpha=0.4)
box(xs, ws, 0, TOP, fc="none")

# F-O reservoir marker (USM sandstone, Valanginian, immediately below 1At1)
yF = y("VA", 0.35)
ax.add_patch(Rectangle((xs + ws * 0.42, y("VA", 0.6)), ws * 0.36, 0.5, fc="none", ec="#c62828", lw=1.4))
ax.annotate("F-O field reservoir\n(Upper Shallow Marine\nsandstone, Valanginian)",
            xy=(xs + ws * 0.60, y("VA") ), xytext=(xs + ws * 0.03, y("HA", 0.85)),
            fontsize=6.2, color="#c62828", fontweight="bold", ha="left", va="bottom",
            bbox=dict(fc="white", ec="#c62828", lw=0.6, pad=1.5),
            arrowprops=dict(arrowstyle="-|>", color="#c62828", lw=0.9))

# ---- depositional environment, tectonism, trap (text boxes spanning stages)
def span_text(col, s0, s1, text, f0=0.0, f1=1.0, fs=6.2):
    x0, w = cols[col]
    yt, yb = y(s0, f0), y(s1, f1)
    box(x0, w, yb, yt, fc="white")
    ax.text(x0 + w / 2, (yt + yb) / 2, text, ha="center", va="center", fontsize=fs, wrap=True)

for col in ["env", "tect", "res", "shows", "trap"]:
    box(cols[col][0], cols[col][1], 0, TOP)
span_text("env", "PL", "CE", "REGRESSIVE\nCYCLES IN\nRESPONSE TO\nEUSTATIC AND\nTECTONIC\nINFLUENCES")
span_text("env", "AL", "HA", "DEEP\nMARINE")
span_text("env", "VA", "BE", "SHALLOW\nMARINE TO\nCONTINENTAL")
span_text("env", "PO", "OXF", "FLUVIAL,\nPOSSIBLY\nLACUSTRINE")
span_text("tect", "PL", "PA", "LAST MAJOR\nUPLIFT")
span_text("tect", "MA", "TU", "THERMAL\nSUBSIDENCE")
span_text("tect", "CE", "CE", "FINAL SEPARATION\nFROM FALKLANDS", fs=5.6)
span_text("tect", "AL", "AL", "THERMAL SUBSIDENCE/\nEUSTATIC EFFECTS", fs=5.4)
span_text("tect", "AP", "AP", "UPLIFT AND\nEROSION", fs=5.8)
span_text("tect", "BA", "HA", "INITIAL OPENING\nAND MARINE\nFLOODING OF BASIN", fs=5.6)
span_text("tect", "VA", "BE", "MAJOR\nREGIONAL\nUPLIFT")
span_text("tect", "PO", "OXF", "INITIAL\nRIFTING")
span_text("trap", "PL", "CE", "")
span_text("trap", "AL", "AL", "FAULT BOUNDED/\nDRAPE/\nSTRATIGRAPHIC", fs=4.9)
span_text("trap", "AP", "BA", "TURBIDITE FAN\nCHANNEL AND\nFAN SANDSTONES", fs=5.6)
span_text("trap", "HA", "HA", "")
span_text("trap", "VA", "BE", "STRUCTURAL\nDOMES")
span_text("trap", "PO", "OXF", "")

# reservoir intervals and shows
xr, wr = cols["res"]; xsh, wsh = cols["shows"]
for s in ["AL", "AP", "BA", "VA"]:
    ax.add_patch(Rectangle((xr + 0.06, y(s, 0.75)), wr - 0.12, 0.5, fc=C["shallow"], ec="k", lw=0.4,
                           hatch="....."))
def gas(xc, yc):
    ax.plot(xc, yc, marker=(8, 2, 0), ms=7, color="k", mew=0.8)
def oil(xc, yc):
    ax.plot(xc, yc, marker="D", ms=4.2, color="k")
oil(xsh + wsh / 2, y("AL", 0.3)); gas(xsh + wsh / 2, y("AL", 0.75))
oil(xsh + wsh / 2, y("AP", 0.35)); gas(xsh + wsh / 2, y("AP", 0.8))
gas(xsh + wsh / 2, y("BA", 0.5))
gas(xsh + wsh / 2, y("VA", 0.5))

# ---- legend
lx, lw_ = 12.75, 2.6
ax.add_patch(Rectangle((lx, 0), lw_, TOP + 1.45, fc="white", ec="k", lw=0.6))
ax.text(lx + lw_ / 2, TOP + 0.72, "LEGEND", ha="center", va="center", fontweight="bold", fontsize=7)
items = [("NON-MARINE", None), ("Alluvial", C["alluvial"]), ("Fluvial", C["fluvial"]),
         ("Lacustrine", C["lacustrine"]), ("COASTAL / MARGINAL /\nSHALLOW MARINE", None),
         ("Deltaic / coastal plain", C["deltaic"]), ("Shallow marine\n(beach, barrier, shoreface)", C["shallow"]),
         ("SUB-WAVE-BASE MARINE", None), ("Shelf", C["shelf"]), ("Slope", C["slope"]),
         ("Deep marine basin floor", C["basin"]), ("GENERAL", None)]
yy = TOP - 0.25
for lab, fc in items:
    if fc is None:
        ax.text(lx + 0.08, yy, lab, fontsize=5.8, fontweight="bold", va="top")
        yy -= 0.55 if "\n" not in lab else 0.85
    else:
        ax.add_patch(Rectangle((lx + 0.1, yy - 0.42), 0.38, 0.38, fc=fc, ec="k", lw=0.4))
        ax.text(lx + 0.58, yy - 0.23, lab, fontsize=5.8, va="center")
        yy -= 0.55 if "\n" not in lab else 0.8
ax.add_patch(Rectangle((lx + 0.1, yy - 0.42), 0.38, 0.38, fc="white", ec="k", lw=0.4, hatch="////"))
ax.text(lx + 0.58, yy - 0.23, "Non-depositional hiatus", fontsize=5.8, va="center"); yy -= 0.55
ax.plot([lx + 0.1, lx + 0.48], [yy - 0.23] * 2, color="k", lw=0.9)
ax.text(lx + 0.58, yy - 0.23, "Major unconformity", fontsize=5.8, va="center"); yy -= 0.55
ax.add_patch(Rectangle((lx + 0.1, yy - 0.42), 0.38, 0.38, fc=C["shallow"], ec="k", lw=0.4, hatch="....."))
ax.text(lx + 0.58, yy - 0.23, "Reservoir sandstone", fontsize=5.8, va="center"); yy -= 0.55
gas(lx + 0.29, yy - 0.23); ax.text(lx + 0.58, yy - 0.23, "Gas", fontsize=5.8, va="center"); yy -= 0.55
oil(lx + 0.29, yy - 0.23); ax.text(lx + 0.58, yy - 0.23, "Oil", fontsize=5.8, va="center"); yy -= 0.55
ax.add_patch(Rectangle((lx + 0.1, yy - 0.42), 0.38, 0.38, fc=C["basement"], ec="k", lw=0.4))
ax.text(lx + 0.58, yy - 0.23, "Basement", fontsize=5.8, va="center"); yy -= 0.55
ax.add_patch(Rectangle((lx + 0.1, yy - 0.42), 0.38, 0.38, fc="none", ec="#c62828", lw=1.4))
ax.text(lx + 0.58, yy - 0.23, "F-O reservoir (this study)", fontsize=5.8, va="center", color="#c62828")

ax.text(0, -0.55, "Modified after Petroleum Agency SA (2004/05), Brown et al. (1995) and McMillan et al. (1997). "
        "Stage boxes are not to time scale; ages (Ma) are those of the source chart.", fontsize=6, style="italic")
fig.savefig("P1_Fig2_chronostratigraphy.png",
            bbox_inches="tight", dpi=300)
fig.savefig("P1_Fig2_chronostratigraphy.pdf",
            bbox_inches="tight")
print("ok")
