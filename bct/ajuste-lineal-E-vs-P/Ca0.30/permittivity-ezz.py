#!/usr/bin/env python3
"""
Compute ε_zz (relative permittivity along z) from polarization-vs-field data
for each temperature, using a linear fit of P_z vs E_z, and plot ε_zz vs T.

Input file:  "polarizaciones"
  Format (whitespace-separated):
    T(K)   E_label   Px   Py   Pz   |P|
  where E_label looks like "0kV_mm", "5kV_mm", "10kV_mm", etc.
  Units:
    - E_label: kV/mm  (1 kV/mm = 1e6 V/m)
    - P components: microC/cm^2 (μC/cm^2)

Outputs:
  - "permittivity-ezz.dat" with columns:
        #T(K)   ezz(rel)   slope_dP_dE(SI)   R^2   chi^2_unweighted
  - "permittivity-ezz.png" figure of ε_zz vs Temperature
  - "fits/fit-TxxxK.png" figures of P_z vs E_z with linear fit for each temperature

Usage:
    python3 permittivity-ezz.py
"""

import re
import sys
from collections import defaultdict
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

INPUT_FILE = "polmacro-Ca0.30"
OUTPUT_FILE = "permittivity-ezz.dat"
PLOT_FILE = "permittivity-ezz.png"
FIT_PLOT_DIR = "fits"

EPS0 = 8.854187817e-12  # F/m


def parse_field_kVmm(label: str) -> float:
    """
    Extract the numeric value in kV/mm from labels like '0kV_mm', '5kV_mm', '10kV_mm'.
    Returns float (kV/mm). Raises ValueError if it cannot parse.
    """
    m = re.match(r'\s*([+-]?\d+(?:\.\d+)?)\s*kV_mm\s*$', label)
    if not m:
        raise ValueError(f"Unrecognized field label: '{label}'")
    return float(m.group(1))


def uc_cm2_to_C_m2(p_microC_cm2: float) -> float:
    """Convert μC/cm^2 to C/m^2"""
    return p_microC_cm2 * 1e-2  # 1 μC/cm^2 = 0.01 C/m^2


def kVmm_to_Vm(e_kVmm: float) -> float:
    """Convert kV/mm to V/m"""
    return e_kVmm * 1e6  # 1 kV/mm = 1e6 V/m


def load_data(path: str):
    """
    Load rows as (T, E(V/m), Pz(C/m^2)).
    Ignores lines that do not have at least 6 columns.
    """
    rows = []
    with open(path, "r") as f:
        for line_no, line in enumerate(f, start=1):
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                sys.stderr.write(
                    f"Warning: line {line_no} has <6 columns, skipping: {line.strip()}\n"
                )
                continue
            try:
                T = float(parts[0])
                E_label = parts[1]
                Pz_uC_cm2 = float(parts[4])  # <-- use Pz (5th column)
                E_kVmm = parse_field_kVmm(E_label)
                E_SI = kVmm_to_Vm(E_kVmm)
                Pz_SI = uc_cm2_to_C_m2(Pz_uC_cm2)
                rows.append((T, E_SI, Pz_SI))
            except Exception as e:
                sys.stderr.write(
                    f"Warning: line {line_no} parse error: {e}; skipping\n"
                )
                continue
    return rows


def compute_ezz_by_temperature(rows):
    """
    For each temperature, fit Pz vs E to get slope dP/dE (SI), then compute ε_zz.
    Also computes R^2 and the unweighted chi^2 = sum(residual^2).
    Returns dict T -> (ezz_rel, slope, r_squared, chi2_unweighted)
    """
    by_T = defaultdict(lambda: {"E": [], "Pz": []})
    for T, E, Pz in rows:
        by_T[T]["E"].append(E)
        by_T[T]["Pz"].append(Pz)

    results = {}
    for T, data in by_T.items():
        E = np.array(data["E"], dtype=float)
        Pz = np.array(data["Pz"], dtype=float)
        if len(E) < 2:
            continue  # Need at least two points
        # Linear fit: Pz = a*E + b  -> slope a = dP/dE
        a, b = np.polyfit(E, Pz, 1)
        fitted = a * E + b
        residuals = Pz - fitted
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((Pz - np.mean(Pz)) ** 2)
        r_squared = 1.0 if np.isclose(ss_tot, 0.0) else 1.0 - (ss_res / ss_tot)
        chi2_unweighted = ss_res
        ezz_rel = 1.0 + (abs(a) / EPS0)
        results[T] = (ezz_rel, a, r_squared, chi2_unweighted)
    return results


def group_rows_by_temperature(rows):
    by_T = defaultdict(lambda: {"E": [], "Pz": []})
    for T, E, Pz in rows:
        by_T[T]["E"].append(E)
        by_T[T]["Pz"].append(Pz)
    return by_T


def write_output(results, path: str):
    Ts = sorted(results.keys(), key=float)
    with open(path, "w") as out:
        out.write(
            "# T(K)    ezz(rel)        slope_dP_dE_SI(C_per_m2_per_V_per_m)"
            "    R_squared    chi2_unweighted\n"
        )
        for T in Ts:
            ezz_rel, slope, r_squared, chi2_unweighted = results[T]
            T_fmt = f"{int(round(T))}" if abs(T - round(T)) < 1e-6 else f"{T:.2f}"
            out.write(
                f"{T_fmt:>6s}  {ezz_rel:12.1f}   {slope: .6e}"
                f"   {r_squared: .6f}   {chi2_unweighted: .6e}\n"
            )


def plot_results(results, path: str):
    Ts = np.array(sorted(results.keys(), key=float), dtype=float)
    ezz = np.array([results[T][0] for T in Ts], dtype=float)

    plt.figure()
    plt.plot(Ts, ezz, marker="o")
    plt.xlabel("Temperature (K)")
    plt.ylabel(r"$\epsilon_{zz}$ (relative)")
    plt.title(r"$\epsilon_{zz}$ vs Temperature")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)

    # Optional phase markers (adjust/remove as needed)
    plt.axvline(110, color="gray", linestyle="--", label="R–O 110 K")
    plt.axvline(140, color="red", linestyle="--", label="O–T 140 K")
    plt.axvline(280, color="blue", linestyle="--", label="T–C 280 K")
    plt.legend()

    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def plot_fits_by_temperature(rows, out_dir: str):
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    by_T = group_rows_by_temperature(rows)
    for T in sorted(by_T.keys(), key=float):
        E = np.array(by_T[T]["E"], dtype=float)
        Pz = np.array(by_T[T]["Pz"], dtype=float)
        if len(E) < 2:
            continue

        order = np.argsort(E)
        E = E[order]
        Pz = Pz[order]
        a, b = np.polyfit(E, Pz, 1)
        fitted = a * E + b

        plt.figure()
        plt.scatter(E, Pz, label="Data")
        plt.plot(E, fitted, color="red", label="Linear fit")
        plt.xlabel(r"$E_z$ (V/m)")
        plt.ylabel(r"$P_z$ (C/m$^2$)")
        plt.title(f"Pz vs Ez at T = {T:g} K")
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.legend()
        plt.tight_layout()

        t_label = f"{int(round(T)):03d}" if abs(T - round(T)) < 1e-6 else f"{T:06.2f}"
        plt.savefig(out_path / f"fit-T{t_label}K.png", dpi=200)
        plt.close()


def main():
    if not Path(INPUT_FILE).exists():
        sys.stderr.write(f"ERROR: '{INPUT_FILE}' not found in current directory.\n")
        sys.exit(1)

    rows = load_data(INPUT_FILE)
    if not rows:
        sys.stderr.write("ERROR: No valid data rows parsed.\n")
        sys.exit(1)

    results = compute_ezz_by_temperature(rows)
    if not results:
        sys.stderr.write(
            "ERROR: Could not compute any ε_zz values "
            "(insufficient points per temperature?).\n"
        )
        sys.exit(1)

    write_output(results, OUTPUT_FILE)
    plot_results(results, PLOT_FILE)
    plot_fits_by_temperature(rows, FIT_PLOT_DIR)
    print(f"Wrote ε_zz results to '{OUTPUT_FILE}' and figure to '{PLOT_FILE}'.")


if __name__ == "__main__":
    main()
