"""Splice tables per well (chosen from the frame inventories). Deeper / later main passes listed first take priority."""
STD_LIS = {"GR": "GR", "SP": "SP", "CALS": "CALS", "CALI": "CALI", "LLD": "LLD", "LLS": "LLS", "MSFL": "MSFL",
           "ILD": "ILD", "ILM": "ILM", "SFLU": "SFLU", "RHOB": "RHOB", "DRHO": "DRHO", "NPHI": "NPHI", "TNPH": "TNPH",
           "NPOR": "NPOR", "PEF": "PEF", "DT": "DT", "DTL": "DTL", "DTCO": "DTCO", "DTSM": "DTSM", "DT4P": "DT4P", "DT4S": "DT4S"}
NUC_DLIS = {"GR_EDTC": "GR", "ECGR": "ECGR", "CALI": "CALI", "DCAL": "DCAL", "LCAL": "LCAL", "RHOB": "RHOB", "RHOM": "RHOB",
            "DRHO": "DRHO", "HDRA": "DRHO", "NPHI": "NPHI", "TNPH": "TNPH", "NPOR": "NPOR", "PEF": "PEF", "PEFL": "PEF", "SP": "SP"}
NMR = {"TCMR": "TCMR", "CMFF": "CMFF", "BFV": "BFV", "KSDR": "KSDR", "KTIM": "KTIM", "T2LM": "T2LM"}

WELLS = {
 "F-O1": dict(top=361.6, base=3927.7, notes="Run-3 density exists twice: RHOB (re175j, 17-Sep-1989) and RHOB_ALT (re175n field copy, ~0.085 g/cc lower); resolve with core.", segs=[
    ("re175n.*", "L0D1", {"GR": "GR", "SP": "SP", "CALS": "CALS", "ILD": "ILD", "ILM": "ILM", "SFLU": "SFLU", "DT": "DT", "DTL": "DTL"}, 3088.7, 3927.7, "ref"),
    ("re175j.*", "L0D1", {"LLD": "LLD", "LLS": "LLS", "RHOB": "RHOB", "DRHO": "DRHO", "NPHI": "NPHI", "TNPH": "TNPH", "NPOR": "NPOR", "PEF": "PEF", "CALI": "CALI"}, 3093.6, 3915.8, "same-pass"),
    ("re175n.*", "L3D1", {"RHOB": "RHOB_ALT"}, 3092.0, 3915.8, "same-pass"),
    ("re175g.*", "L0D1", {"GR": "GR", "SP": "SP", "CALS": "CALS", "ILD": "ILD", "ILM": "ILM", "SFLU": "SFLU", "DT": "DT"}, 1494.3, 3088.6, "run2"),
    ("re175q.*", "L1D1", {"RHOB": "RHOB", "DRHO": "DRHO", "PEF": "PEF", "CALI": "CALI"}, 1495.3, 3093.5, "run2"),
    ("re175b.*", "L0D1", {"GR": "GR", "SP": "SP", "CALS": "CALS", "ILD": "ILD", "ILM": "ILM", "SFLU": "SFLU", "DT": "DT"}, 361.6, 1494.2, "run1")]),
 "F-O2": dict(top=1484.7, base=3900.7, notes="9-5/8 in shoe ~3446 m: run-3 data above it (LLS 0.15, DT 59.6 = steel) excluded; run 2b fills 3180-3446 m. Tape re210.fo2_6_7_9 is ~8.4 m deeper than all other F-O2 passes and was NOT used. No deep resistivity 3180-3392 m (not logged on supplied tapes).", segs=[
    ("re210d.*", "L0D1", {"GR": "GR", "SP": "SP", "CALS": "CALS", "LLD": "LLD", "LLS": "LLS", "DT": "DT", "DTL": "DTL"}, 3447.0, 3897.8, "ref"),
    ("re210f.*", "L1D1", {"CALI": "CALI", "RHOB": "RHOB", "DRHO": "DRHO", "NPHI": "NPHI", "TNPH": "TNPH", "NPOR": "NPOR", "PEF": "PEF"}, 3447.0, 3900.7, "run3"),
    ("re210j*", "L0D1", {"DTCO": "DTCO", "DTSM": "DTSM"}, 3447.0, 3888.9, "run3"),
    ("re210c*", "L2D1", {"GR": "GR", "CALI": "CALI", "RHOB": "RHOB", "DRHO": "DRHO", "PEF": "PEF", "DT": "DT", "DTL": "DTL"}, 3180.4, 3466.3, "run2b"),
    ("re210c*", "L0D1", {"GR": "GR", "SP": "SP", "CALI": "CALI", "CALS": "CALS", "LLD": "LLD", "LLS": "LLS", "RHOB": "RHOB", "DRHO": "DRHO", "PEF": "PEF", "DT": "DT", "DTL": "DTL"}, 1484.7, 3180.3, "run2")]),
 "F-O3": dict(top=1469.0, base=3958.9, notes="No deep resistivity 3195-3399 m on supplied tapes. 9-5/8 in shoe 3398 m: run 3 (re231f/c/g) above it was logged inside casing (CALI 8.6 in, RHOB ~2.05, DRHO ~-0.95) and is excluded; open-hole run 2 used to 3195 m; re231c sonic 3195-3399 m kept as cased-hole DT4P/DT4S (formation-like values), its GR as GR_CH.", segs=[
    ("re231g*", "L0D1", {"GR": "GR", "SP": "SP", "CALS": "CALS", "LLD": "LLD", "LLS": "LLS", "DTCO": "DTCO", "DTSM": "DTSM", "DT4P": "DT4P", "DT4S": "DT4S"}, 3399.0, 3958.9, "ref"),
    ("re231f*", "L1D1", {"CALI": "CALI", "RHOB": "RHOB", "DRHO": "DRHO", "NPHI": "NPHI", "TNPH": "TNPH", "NPOR": "NPOR", "PEF": "PEF", "GR": "GR"}, 3399.0, 3958.3, "run3"),
    ("re231c*", "L1D1", {"GR": "GR_CH", "DT4P": "DT4P", "DT4S": "DT4S"}, 3195.3, 3398.9, "run3b-cased"),
    ("re231a*", "L0D1", {"GR": "GR", "SP": "SP", "CALI": "CALI", "CALS": "CALS", "LLD": "LLD", "LLS": "LLS", "RHOB": "RHOB", "DRHO": "DRHO", "PEF": "PEF", "DT": "DT", "DTL": "DTL"}, 1469.0, 3195.2, "run2")]),
 "F-O4": dict(top=160.0, base=3905.4, notes="MSFL resampled from 2-inch frame; NMR (CMR) from PLAY_071 3569-3884 m.", segs=[
    ("PLAY_047*", "L0:60B", {"GR": "GR", "ECGR": "ECGR", "SP": "SP", "CALS": "CALS", "LLD": "LLD", "LLS": "LLS", "DT": "DT", "DTL": "DTL", "DTCO": "DTCO", "DTSM": "DTSM"}, 3085.6, 3905.4, "ref"),
    ("PLAY_047*", "L0:20B", {"MSFL": "MSFL"}, 3085.7, 3905.4, "same-pass"),
    ("PLAY_061*", "L0:60B", NUC_DLIS, 3110.3, 3889.6, "run"),
    ("PLAY_071*", "L0:30B", NMR, 3569.0, 3884.4, "nmr"),
    ("PLAY_049*", "L0:60B", {"GR": "GR", "ECGR": "ECGR", "SP": "SP", "CALS": "CALS", "LLD": "LLD", "LLS": "LLS", "DT": "DT", "DTL": "DTL"}, 160.0, 3085.5, "run1"),
    ("PLAY_049*", "L0:20B", {"MSFL": "MSFL"}, 160.1, 3085.6, "run1")]),
 "F-O6": dict(top=2700.1, base=3929.6, notes="Resistivity: AIT (AT10-AT90, 2700-3842 m) and HRLA (RLA1-RLA5, 3595-3930 m). Spectral GR as HFK (%), HTHO, HURA (ppm). Density mnemonic RHOM renamed RHOB, PEFL renamed PEF.", segs=[
    ("F_O6_MainLogDept*dlis*", "L0:60B", {"GR_EDTC": "GR", "HFK": "HFK", "HTHO": "HTHO", "HURA": "HURA", "SP": "SP", "DCAL": "DCAL", "LCAL": "LCAL", "RHOM": "RHOB", "NPHI": "NPHI", "TNPH": "TNPH", "NPOR": "NPOR", "PEFL": "PEF"}, 2700.1, 3842.3, "ref"),
    ("F_O6_MainLogDept*dlis*", "L0:30B", {"AT90": "AT90", "AT60": "AT60", "AT30": "AT30", "AT10": "AT10", "AF90": "AF90", "AO90": "AO90"}, 2700.1, 3842.3, "same-pass"),
    ("F_O6_ML_TD3920m*DLIS*", "L0:60B", {"GR_EDTC": "GR", "HFK": "HFK", "HTHO": "HTHO", "HURA": "HURA", "SP": "SP", "DCAL": "DCAL", "LCAL": "LCAL", "RHOM": "RHOB", "NPHI": "NPHI", "TNPH": "TNPH", "NPOR": "NPOR", "PEFL": "PEF"}, 3595.0, 3929.6, "run-td"),
    ("F_O6_ML_TD3920m*DLIS*", "L0:20B", {"RLA5": "RLA5", "RLA3": "RLA3", "RLA1": "RLA1"}, 3595.0, 3929.6, "run-td"),
    ("F-O6_SoniScan*", "L11:", {"DTCO": "DTCO", "DT4P": "DT4P"}, 2745.9, 3851.9, "sonic"),
    ("F-O6_SoniScan*", "L10:", {"DTSM": "DTSM"}, 2745.9, 3851.9, "sonic")]),
 "F-O8": dict(null=[("GR",3900.5,3916),("NPHI",3906.6,3916),("RHOB",3901.5,3903.3)], top=1143.3, base=3915.9, notes="Sonic Scanner DTCO/DTSM: Main_1 (073PUC) 3750-3899 m; Main_2 (074PUC) 3430-3750 m with a -0.30 m bulk shift (windowed GR lag vs Run 2A -0.15 to -0.46 m; shifted Main_2 vs Main_1 DTCO r=0.99 over 3750-3762 m). Held (constant) readings removed: Main_1 3737-3750 and >3899 m, Main_2 <3430 and >3762 m, Run 2A GR >3900.5, NPHI >3906.6, RHOB 3901.6-3903.3 m. GR_SS2 = shifted Main_2 gamma ray for QC. 1143-2724 m (I_C_SS run 2B pass 4) was logged INSIDE 9-5/8 in casing (C1/C2 ~5 in, SOBS 8.64 in, DTCO 56.5 us/ft = steel): its DTCO/DTSM are casing arrivals and are excluded; its gamma ray is kept only as GR_CH (through-casing, attenuated). HRLA RLA1-RLA5; spectral GR POTA/THOR/URAN from NGT run 2C; CMR from run 2C (DMRP). RHOM->RHOB, PEFL->PEF.", segs=[
    ("F-O8_Res_Nuclear_Run2A_Main_058PUC.DLIS*", "L0:60B", {"GR_EDTC": "GR", "SP": "SP", "DCAL": "DCAL", "LCAL": "LCAL", "RHOM": "RHOB", "NPHI": "NPHI", "TNPH": "TNPH", "NPOR": "NPOR", "PEFL": "PEF"}, 3385.9, 3915.9, "ref"),
    ("F-O8_Res_Nuclear_Run2A_Main_058PUC.DLIS*", "L0:20B", {"RLA5": "RLA5", "RLA3": "RLA3", "RLA1": "RLA1"}, 3385.9, 3915.9, "same-pass"),
    ("F-O8_NGT_Run2C_Main_096PUC.DLIS*", "L0:60B", {"SGR": "SGR", "CGR": "CGR", "POTA": "POTA", "THOR": "THOR", "URAN": "URAN"}, 3626.4, 3915.2, "run2c"),
    ("F-O8_CMR_Run2C_Main_DMRP_091PUC.DLIS*", "L0:60B", NMR, 3689.9, 3915.9, "run2c"),
    ("F-O8_SS_Run2B_Main_1_OH_073PUC.DLIS*", "L0:60B", {"DTCO": "DTCO", "DTSM": "DTSM"}, 3750.1, 3899.3, "run2b"),
    ("F-O8_SS_Run2B_Main_2_OH_074PUC.DLIS*", "L0:60B", {"DTCO": "DTCO", "DTSM": "DTSM", "GR_EDTC": "GR_SS2"}, 3429.9, 3761.8, "run2b-main2@shift=-0.3048"),
    ("F-O8_I_C_SS_Run2B_4*", "las", {"GR_EDTC": "GR_CH"}, 1143.3, 2723.7, "run2b-cased")]),
 "F-R1": dict(top=232.0, base=4572.0, notes="1980 wireline (DIL/ISF-type: ILD, SFLU, MSFL) and FDC-CNL; no PEF.", segs=[
    ("F-R1-re023e*", "L0D1", {"GR": "GR", "SP": "SP", "CALI": "CALI", "MSFL": "MSFL", "ILD": "ILD", "SFLU": "SFLU", "DT": "DT"}, 232.0, 4572.0, "ref"),
    ("F-R1-re023e*", "L1D1", {"RHOB": "RHOB", "DRHO": "DRHO", "NPHI": "NPHI"}, 232.0, 4572.0, "same-file"),
    ("F-R1-re023d.fr1_prop_tlis_3274", "L0D1", {"GR": "GR", "SP": "SP", "CALI": "CALI", "MSFL": "MSFL", "ILD": "ILD", "SFLU": "SFLU", "DT": "DT"}, 4254.7, 4581.1, "deep-run")]),
 "F-S1": dict(top=2196.8, base=3960.4, notes="9-5/8 in shoe ~3500 m: run-3 data above it (DT ~56 = steel, CALI 8.6, DRHO -0.86) excluded. Dipole sonic (DTCO/DTSM) 3509-3959 m from re225l; nuclear from re225m2 (re225m1 is corrupt: PEF~276, CALI~2.8 in).", segs=[
    ("F-S1-re225i*", "L0D1", {"GR": "GR", "SP": "SP", "CALS": "CALS", "LLD": "LLD", "LLS": "LLS", "MSFL": "MSFL", "DT": "DT"}, 3501.0, 3949.9, "ref"),
    ("F-S1-re225m2*", "L0:60B", NUC_DLIS, 3501.0, 3960.4, "run"),
    ("F-S1-re225l_076*", "L0:60B", {"DTCO": "DTCO", "DTSM": "DTSM", "DT4P": "DT4P", "DT4S": "DT4S"}, 3509.3, 3959.0, "run"),
    ("F-S1-re225b*", "L0D1", {"GR": "GR", "SP": "SP", "CALI": "CALI", "CALS": "CALS", "LLD": "LLD", "LLS": "LLS", "RHOB": "RHOB", "DRHO": "DRHO", "PEF": "PEF", "DT": "DT", "DTL": "DTL"}, 2196.8, 3429.5, "run2")]),
}
