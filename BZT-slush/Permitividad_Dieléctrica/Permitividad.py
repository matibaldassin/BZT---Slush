import os
import numpy as np
import matplotlib.pyplot as plt

# ---------------- Parámetros ----------------
Zr = "0.60"
UMBRAL_NO_LINEAL = 0.15     # desviación relativa entre pendientes de tramos
EPS0 = 8.8541878128e-12      # F/m
UC_A_SI = 1e-2               # 1 uC/cm^2 -> C/m^2
KV_MM_A_SI = 1e6             # 1 kV/mm -> V/m

base = os.path.dirname(os.path.abspath(__file__))
dir_datos = os.path.join(base, "..", "Datos-corridas", "Permitividad_Dieléctrica")
dir_tablas = os.path.join(dir_datos, "Tablas_Permitividad")
dir_pe = os.path.join(base, f"PvsE_Zr{Zr}")
os.makedirs(dir_tablas, exist_ok=True)
os.makedirs(dir_pe, exist_ok=True)

archivo = os.path.join(dir_datos, f"Pol_Macro_Permitividad_Zr{Zr}.txt")

# ---------------- Lectura ----------------
# columnas: T, campo, Px, Py, Pz, |P|
datos = {}
with open(archivo) as f:
    for linea in f:
        c = linea.split()
        if len(c) < 6:
            continue
        T = float(c[0])
        E = float(c[1].replace("kV_mm", ""))   # kV/mm
        Pz = float(c[4])                        # uC/cm^2
        datos.setdefault(T, []).append((E, Pz))

# ---------------- Cálculo ----------------
filas = []
for T in sorted(datos):
    E, Pz = np.array(sorted(datos[T])).T
    E_si = E * KV_MM_A_SI
    P_si = Pz * UC_A_SI

    # ajuste lineal P = s*E + P0
    (s, P0), cov = np.polyfit(E_si, P_si, 1, cov=True) if len(E) > 2 else (np.polyfit(E_si, P_si, 1), None)
    err_s = np.sqrt(cov[0, 0]) if cov is not None else np.nan
    res = P_si - (s * E_si + P0)
    ss_tot = np.sum((P_si - P_si.mean()) ** 2)
    r2 = 1 - np.sum(res ** 2) / ss_tot if ss_tot > 0 else np.nan

    chi = s / EPS0
    eps_r = 1 + chi
    err_eps = err_s / EPS0

    # pendientes por tramo y desviación relativa
    s1 = (P_si[1] - P_si[0]) / (E_si[1] - E_si[0])
    s2 = (P_si[2] - P_si[1]) / (E_si[2] - E_si[1])
    s_med = 0.5 * (s1 + s2)
    desv = abs(s2 - s1) / abs(s_med) if s_med != 0 else np.inf

    # curvatura (coef. cuadrático del ajuste exacto de grado 2), en C/m^2 / (kV/mm)^2
    curv = np.polyfit(E, P_si, 2)[0]

    lineal = desv < UMBRAL_NO_LINEAL
    filas.append((T, s, eps_r, err_eps, r2, s1 / EPS0, s2 / EPS0, desv, curv, lineal, E, Pz, P0))

# ---------------- Tabla ----------------
ruta_tabla = os.path.join(dir_tablas, f"Permitividad_Zr{Zr}.txt")
with open(ruta_tabla, "w") as f:
    f.write(f"{'T(K)':>6} {'dP/dE(F/m)':>12} {'eps_r':>9} {'err_eps':>9} {'R2':>7} "
            f"{'chi_0-5':>9} {'chi_5-10':>9} {'desv_tramos':>12} {'curv':>11} {'lineal':>7}\n")
    for T, s, eps, err, r2, c1, c2, desv, curv, lin, *_ in filas:
        f.write(f"{T:6.0f} {s:12.4e} {eps:9.2f} {err:9.2f} {r2:7.4f} "
                f"{c1:9.2f} {c2:9.2f} {desv:12.3f} {curv:11.3e} {'si' if lin else 'NO':>7}\n")
print(f"Tabla guardada en {ruta_tabla}")

# ---------------- Gráficos P vs E ----------------
n = len(filas)
ncol = 4
nrow = int(np.ceil(n / ncol))
fig, axs = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.2 * nrow), sharex=True)
axs = axs.ravel()
Efit = np.linspace(0, 10, 50)
for ax, (T, s, eps, err, r2, c1, c2, desv, curv, lin, E, Pz, P0) in zip(axs, filas):
    ax.plot(E, Pz, "o", color="C0")
    ax.plot(Efit, (s * Efit * KV_MM_A_SI + P0) / UC_A_SI, "-", color="C2" if lin else "C3")
    ax.set_title(f"T = {T:.0f} K  ({'lineal' if lin else 'NO lineal'})", fontsize=10)
    ax.text(0.05, 0.95, f"$R^2$={r2:.3f}\n$\\Delta s/\\bar s$={desv:.2f}\n$\\varepsilon_r$={eps:.1f}",
            transform=ax.transAxes, va="top", fontsize=8)
    ax.set_xticks([0, 5, 10])
    ax.grid(alpha=0.3)
for ax in axs[n:]:
    ax.axis("off")
for ax in axs[::ncol]:
    ax.set_ylabel(r"$P_z$ ($\mu$C/cm$^2$)")
for ax in axs[-ncol:]:
    ax.set_xlabel("E (kV/mm)")
fig.tight_layout()
fig.savefig(os.path.join(dir_pe, f"PvsE_Zr{Zr}_todas.png"), dpi=200)
plt.close(fig)

# un gráfico individual por temperatura
for T, s, eps, err, r2, c1, c2, desv, curv, lin, E, Pz, P0 in filas:
    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    ax.plot(E, Pz, "o", color="C0", label="datos")
    ax.plot(Efit, (s * Efit * KV_MM_A_SI + P0) / UC_A_SI, "-", color="C2" if lin else "C3", label="ajuste lineal")
    ax.set_xlabel("E (kV/mm)")
    ax.set_ylabel(r"$P_z$ ($\mu$C/cm$^2$)")
    ax.set_title(f"T = {T:.0f} K, $\\varepsilon_r$ = {eps:.1f}, $R^2$ = {r2:.3f}, "
                 f"{'lineal' if lin else 'NO lineal'}", fontsize=9)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(dir_pe, f"PvsE_Zr{Zr}_T{T:.0f}K.png"), dpi=200)
    plt.close(fig)

# ---------------- Gráfico eps_r vs T ----------------
T_arr = np.array([f[0] for f in filas])
eps_arr = np.array([f[2] for f in filas])
err_arr = np.array([f[3] for f in filas])
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.errorbar(T_arr, eps_arr, yerr=err_arr, fmt="o-", color="C0", capsize=3)
ax.set_xlabel("T (K)")
ax.set_ylabel(r"$\varepsilon_r$")
ax.set_title(f"Permitividad dieléctrica, Zr = {Zr}")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(base, f"Permitividad_vs_T_Zr{Zr}.png"), dpi=200)
plt.close(fig)
print("Gráficos guardados.")

# ---------------- Gráfico eps_r vs T: todas las composiciones ----------------
ZRS = ["0.40", "0.50", "0.60"]
N_SUAVIZADO = 3              # ventana del promedio móvil; 1 = sin suavizar


def suavizar(y, n=N_SUAVIZADO):
    """Promedio móvil centrado de ventana n; en los bordes usa los puntos disponibles."""
    k = np.ones(n)
    return np.convolve(y, k, mode="same") / np.convolve(np.ones(len(y)), k, mode="same")


fig, ax = plt.subplots(figsize=(6.5, 4.5))
for z, marcador in zip(ZRS, ["o", "s", "^"]):
    tabla = np.loadtxt(os.path.join(dir_tablas, f"Permitividad_Zr{z}.txt"), skiprows=1, usecols=(0, 2))
    ax.plot(tabla[:, 0], suavizar(tabla[:, 1]), f"{marcador}-", label=f"Zr = {z}")
ax.set_xlabel("T (K)")
ax.set_ylabel(r"$\varepsilon_r$")
ax.set_title("Permitividad dieléctrica")
ax.grid(alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(base, "Permitividad_vs_T_todas.png"), dpi=200)
plt.close(fig)
print("Gráfico comparativo guardado.")
