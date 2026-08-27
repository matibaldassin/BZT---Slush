#!/usr/bin/env python3
import re
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import os

# Cambiar al directorio del script para que los archivos relativos se encuentren
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

CONC='0.30'

INPUT_FILE   = f"polmacro-Ca{CONC}"
TABLE_DIA    = f"permittivity-Ca{CONC}-table.dat"

# Unit conversions
U_P_UC_CM2_TO_C_M2 = 1e-2   # μC/cm^2 -> C/m^2
U_E_KV_MM_TO_V_M   = 1e6    # kV/mm   -> V/m
EPS0 = 8.854187817e-12      # F/m

# Regex to parse "5kV_mm"
field_re = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)\s*kV_mm\s*$")

# Data holders by temperature
by_T_Px  = defaultdict(list)   # (E, Px)
by_T_Py  = defaultdict(list)   # (E, Py)
by_T_Pz  = defaultdict(list)   # (E, Pz)

with open(INPUT_FILE, "r") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        # Expect: T, "EkV_mm", Px, Py, Pz, |P|
        if len(parts) < 6:
            continue

        try:
            T_str, E_str = parts[0], parts[1]
            Px_str, Py_str, Pz_str = parts[2], parts[3], parts[4]

            T = int(float(T_str))

            m = field_re.match(E_str)
            if not m:
                continue
            E_kV_mm = float(m.group(1))
            E_V_m = E_kV_mm * U_E_KV_MM_TO_V_M

            Px = float(Px_str) * U_P_UC_CM2_TO_C_M2
            Py = float(Py_str) * U_P_UC_CM2_TO_C_M2
            Pz = float(Pz_str) * U_P_UC_CM2_TO_C_M2

            # Projection onto [111]
            P_parallel = (Px + Py + Pz) / np.sqrt(3.0)

            by_T_Px[T].append((E_V_m, Px))
            by_T_Py[T].append((E_V_m, Py))
            by_T_Pz[T].append((E_V_m, Pz))

        except ValueError:
            continue

def fit_eps_vs_E(data_pairs):
    """
    data_pairs: list of (E, P_component) in SI.
    Returns (eps_r, slope, intercept) where slope = dP/dE [F/m], eps_r = 1 + slope/eps0
    """
    data_pairs = sorted(data_pairs, key=lambda x: x[0])
    E = np.array([d[0] for d in data_pairs], dtype=float)
    P = np.array([d[1] for d in data_pairs], dtype=float)
    slope, intercept = np.polyfit(E, P, 1)

    if slope < 0:
        slope = -slope

    eps_r = 1.0 + slope / EPS0
    return eps_r, slope, intercept

# ---- exx / eyy / ezz from Px, Py, Pz vs E ----
fit_diag = []  # (T, eps_xx, eps_yy, eps_zz, slope_x, slope_y, slope_z)

all_Ts = sorted(set(by_T_Px.keys()) | set(by_T_Py.keys()) | set(by_T_Pz.keys()))
for T in all_Ts:
    # Require at least 2 points per component to fit
    if len(by_T_Px[T]) < 2 or len(by_T_Py[T]) < 2 or len(by_T_Pz[T]) < 2:
        continue
    eps_xx, slope_x, _ = fit_eps_vs_E(by_T_Px[T])
    eps_yy, slope_y, _ = fit_eps_vs_E(by_T_Py[T])
    eps_zz, slope_z, _ = fit_eps_vs_E(by_T_Pz[T])

    fit_diag.append((T, eps_xx, eps_yy, eps_zz, slope_x, slope_y, slope_z))

# Save exx/eyy/ezz table
with open(TABLE_DIA, "w") as fout:
    fout.write("# T(K)    exx(rel)    eyy(rel)    ezz(rel)    slope_dPx_dE_SI    slope_dPy_dE_SI    slope_dPz_dE_SI\n")
    for T, exx, eyy, ezz, sx, sy, sz in sorted(fit_diag):
        fout.write(f"{T:6d}  {exx:10.1f}  {eyy:10.1f}  {ezz:10.1f}    {sx:.6e}    {sy:.6e}    {sz:.6e}\n")