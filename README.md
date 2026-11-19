# F-O field tight-gas petrophysics workflow

Open-source Python workflow behind two companion papers on the F-O gas field, Bredasdorp Basin, offshore South Africa:

1. Amaechi, P.O., Afolayan, B.A., Opuwari, M. *Diagenetic controls on permeability and net pay in a tight gas reservoir: an integrated core–log study of the F-O field, Bredasdorp Basin, offshore South Africa.* (submitted)
2. Amaechi, P.O., Afolayan, B.A., Opuwari, M. *Reconciling log, capillary-pressure and formation-pressure data for gas saturation and volumetric uncertainty in a tight, overpressured reservoir: the F-O field, Bredasdorp Basin, offshore South Africa.* (submitted)

The code regenerates every number, table and figure in both papers from the original field tapes and core reports. **No well, core or pressure data are included** (see [`data/README.md`](data/README.md)). The data belong to the Petroleum Agency SA and are available from the Agency on request.

The deterministic steps (shale volume, porosity, porosity–permeability transforms, Rw checks, Waxman–Smits saturation, net pay, core rock types, pressure–depth plot and saturation–height model) were also reproduced independently in Interactive Petrophysics (IP) 2026. The comparison is given in Supplementary Table S1 of Paper 1.

## Repository layout

| Folder | Content |
|---|---|
| `scripts/composite_build/` | Tape decoding (LIS/DLIS), splicing, log clean-up, core-report OCR parsing, core–log depth matching, pretest extraction from tapes, capillary-pressure assembly. Produces the well composites (v2–v3). Run `run_all.py` after setting the paths at the top of each script. |
| `notebooks/01_lis_dlis_parsing.ipynb` | Worked example of reading the legacy LIS tapes (F-O1). |
| `notebooks/02_…` to `07_…` | Core–log QC; shale volume, porosity and permeability; Rw, m\* calibration, saturation and net pay; pressures, J-function, saturation–height model and Monte Carlo; rock typing and mineral model; machine-learning permeability (Paper 3, in preparation). Edit the paths in Cell 1 only. |
| `petrotools/` | Library used by the notebooks: `saturation.py`, `petromodels.py`, `mineral.py`, `rocktype.py`, `shm.py` (datums, GWC), `pressure_extract.py`, `coretools.py`, `mlperm.py`, readers and QC plots. |
| `scripts/petrography/` | Figure 11 of Paper 1 (thin-section and modal data). |
| `scripts/residual_gas/` | Neutron–density residual-gas test (Paper 2, Section 5.1). |
| `scripts/figures/` | Paper figure scripts, including the QGIS location maps (Figure 1 of both papers; needs the georeferenced PetroSA depth map, which is not distributed) and the chronostratigraphic chart (Paper 1, Figure 2). |

## Running

```bash
pip install -r requirements.txt
# 1. place the PASA data as described in data/README.md
# 2. build the composites
python scripts/composite_build/run_all.py
# 3. run notebooks 02 to 07 in order (each writes the composite version used by the next)
```

Tested with Python 3.12, NumPy 2.5, SciPy 1.18, pandas 3.0, scikit-learn 1.9, Matplotlib 3.11, lasio 0.32 and dlisio 1.0.

Some scripts in `scripts/composite_build/` were written for a single working session and still use relative paths such as `work/` and `data/raw/`; set these to your own folders before running.

## Licence and citation

Code: MIT licence (see `LICENSE`). Please cite the papers above and this repository (`CITATION.cff`; a Zenodo DOI will be added at release).
