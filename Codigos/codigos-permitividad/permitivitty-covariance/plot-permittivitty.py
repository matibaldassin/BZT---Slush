#!/usr/bin/env python3
# plot-permittivity.py
import pandas as pd
import matplotlib.pyplot as plt

CSV_PATH = "e-vs-T.csv"  # change to "epsilon_vs_T.csv" if that's your file

# Column sets
cols20 = [
    "file","T_K",
    "eps_xx","eps_yy","eps_zz",
    "eps_xy","eps_xz","eps_yz","eps_iso",
    "eps_xx_err","eps_yy_err","eps_zz_err",
    "eps_xy_err","eps_xz_err","eps_yz_err","eps_iso_err",
    "eps_parallel","eps_perp","eps_parallel_err","eps_perp_err"
]
cols16 = [
    "file","T_K",
    "eps_xx","eps_yy","eps_zz",
    "eps_xy","eps_xz","eps_yz","eps_iso",
    "eps_xx_err","eps_yy_err","eps_zz_err",
    "eps_xy_err","eps_xz_err","eps_yz_err","eps_iso_err"
]

def read_flex(path):
    """
    Try headered CSV first. If parsing fails or column count mismatches,
    fall back to a tolerant reader that accepts ragged rows (old/new mix).
    """
    try:
        df = pd.read_csv(path)  # header present
        return df
    except Exception:
        # Mixed-width fallback: force 20 columns, skip header if present
        try:
            df = pd.read_csv(
                path,
                header=None,
                names=cols20,
                engine="python"
            )
            # If the file had a header row, drop it (it will be the first row with non-numeric T_K)
            if pd.to_numeric(df["T_K"], errors="coerce").isna().iloc[0]:
                df = df.iloc[1:].reset_index(drop=True)
            return df
        except Exception as e:
            raise RuntimeError(f"Could not read CSV '{path}': {e}")

df = read_flex(CSV_PATH)

# If it's a 16-col file, trim to the 16-column schema
if set(cols20).issubset(df.columns):
    schema = cols20
elif set(cols16).issubset(df.columns):
    schema = cols16
else:
    # Try aligning: if df came from fallback, it already has cols20; otherwise coerce
    if len(df.columns) == 20:
        df.columns = cols20
        schema = cols20
    elif len(df.columns) == 16:
        df.columns = cols16
        schema = cols16
    else:
        raise RuntimeError(f"Unexpected columns in CSV ({len(df.columns)}).")

# Extract basics
T    = pd.to_numeric(df["T_K"], errors="coerce")
ex   = pd.to_numeric(df["eps_xx"], errors="coerce")
ey   = pd.to_numeric(df["eps_yy"], errors="coerce")
ez   = pd.to_numeric(df["eps_zz"], errors="coerce")
eiso = pd.to_numeric(df["eps_iso"], errors="coerce")

# --- Figure 1: tensor components + iso ---
plt.figure(figsize=(8,6))
plt.plot(T, ex, "o-", label=r"$\epsilon_{xx}$")
plt.plot(T, ey, "s-", label=r"$\epsilon_{yy}$")
plt.plot(T, ez, "d-", label=r"$\epsilon_{zz}$")
plt.plot(T, eiso, "^-", label=r"$\epsilon_{\mathrm{iso}}$", linewidth=2, markersize=7)

# vertical dashed line (as in your script)
plt.axvline(280, color="blue", linestyle="--", label="T–C 270 K")

plt.xlabel("Temperature (K)")
plt.ylabel(r"Permittivity $\epsilon$")
plt.title("Permittivity vs Temperature")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("permittivity_vs_T.png", dpi=300)
plt.show()
print("Saved figure -> permittivity_vs_T.png")

# --- Figure 2: parallel / perpendicular if present ---
if "eps_parallel" in df.columns and "eps_perp" in df.columns:
    epar     = pd.to_numeric(df["eps_parallel"], errors="coerce")
    eperp    = pd.to_numeric(df["eps_perp"], errors="coerce")
    epar_err = pd.to_numeric(df["eps_parallel_err"], errors="coerce") if "eps_parallel_err" in df.columns else None
    eperp_err= pd.to_numeric(df["eps_perp_err"],  errors="coerce") if "eps_perp_err"  in df.columns else None

    plt.figure(figsize=(8,6))
    if epar_err is not None and eperp_err is not None and epar_err.notna().any() and eperp_err.notna().any():
        plt.errorbar(T, epar,  yerr=epar_err,  fmt="o-", capsize=3, label=r"$\epsilon_{\parallel}$")
        plt.errorbar(T, eperp, yerr=eperp_err, fmt="s-", capsize=3, label=r"$\epsilon_{\perp}$")
    else:
        plt.plot(T, epar,  "o-", label=r"$\epsilon_{\parallel}$")
        plt.plot(T, eperp, "s-", label=r"$\epsilon_{\perp}$")

    plt.axvline(280, color="blue", linestyle="--", label="T–C 270 K")
    plt.xlabel("Temperature (K)")
    plt.ylabel(r"Permittivity $\epsilon$")
    plt.title(r"Longitudinal / Transverse Permittivity vs Temperature")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("permittivity_parallel_perp_vs_T.png", dpi=300)
    plt.show()
    print("Saved figure -> permittivity_parallel_perp_vs_T.png")
else:
    print("Note: columns 'eps_parallel' / 'eps_perp' not found; second figure skipped.")
