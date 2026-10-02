import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / ".." / "Datos-corridas" / "P-local"

# Composición -> cantidad de celdas de Zr (son las primeras de cada archivo)
SISTEMAS = {
    "0.40": 8780,
    "0.50": 10976,
    "0.60": 13172,
}

OUT_FILE = SCRIPT_DIR / "P-local-promedio.dat"
OUT_IMAGE = SCRIPT_DIR / "P-local-promedio.png"
OUT_IMAGE_TI_ZR = SCRIPT_DIR / "P-local-promedio-Ti-Zr.png"

FACTOR_T = 1.48

COLORES_ZR = ["#2a78d6", "#eb6834", "#1baf7a"]
COLOR_TI = "#2a78d6"
COLOR_ZR = "#eb6834"


def numero_archivo(file_path):
    # P-local-<n>
    return int(file_path.name.split("-")[-1])


def archivos(zr):
    return sorted((DATA_DIR / f"Zr{zr}").glob("P-local-*"), key=numero_archivo)


def p_rms(P_tot):
    # Módulo promedio de la polarización local: sqrt(<|P|^2>)
    return np.sqrt(np.mean(P_tot**2))


def read_data_separeted(file_path, cells_Zr):
    data = np.loadtxt(file_path)
    P_tot = data[:, 4]

    P_tot_Zr_avg = p_rms(P_tot[:cells_Zr])
    P_tot_Ti_avg = p_rms(P_tot[cells_Zr:])
    P_tot_avg = p_rms(P_tot)

    return P_tot_Zr_avg, P_tot_Ti_avg, P_tot_avg


def calcular(zr, cells_Zr):
    files = archivos(zr)
    T = FACTOR_T * np.array([numero_archivo(f) for f in files])
    P_Zr, P_Ti, P_all = np.array([read_data_separeted(f, cells_Zr) for f in files]).T
    return T, P_all, P_Zr, P_Ti


def estilo(ax):
    ax.grid(color="0.9", linewidth=0.8)
    ax.set_axisbelow(True)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)


resultados = {zr: calcular(zr, cells_Zr) for zr, cells_Zr in SISTEMAS.items()}

# Guardar datos: T y, para cada composición, P_total, P_Zr, P_Ti
T = resultados["0.40"][0]
columnas = [T]
encabezado = "T(K)"
for zr, (_, P_all, P_Zr, P_Ti) in resultados.items():
    columnas += [P_all, P_Zr, P_Ti]
    encabezado += f"  P_tot(Zr{zr})  P_Zr(Zr{zr})  P_Ti(Zr{zr})"
np.savetxt(OUT_FILE, np.column_stack(columnas), fmt="%10.4f", header=encabezado)

# Gráfica 1: todos los sitios
fig, ax = plt.subplots(figsize=(6, 4.5))
for (zr, (T, P_all, _, _)), color in zip(resultados.items(), COLORES_ZR):
    ax.plot(T, P_all, "-o", color=color, lw=2, ms=4, label=f"x = {zr}")
ax.set_xlabel("T (K)")
ax.set_ylabel(r"$P_L$")
ax.set_title("Polarización local promedio (todos los sitios)")
ax.legend(title=r"Ba(Zr$_x$Ti$_{1-x}$)O$_3$", frameon=False)
estilo(ax)
fig.tight_layout()
fig.savefig(OUT_IMAGE, dpi=200)

# Gráfica 2: sitios de Ti y de Zr, un panel por composición
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
for ax, (zr, (T, _, P_Zr, P_Ti)) in zip(axes, resultados.items()):
    ax.plot(T, P_Ti, "-o", color=COLOR_TI, lw=2, ms=4, label="Sitios Ti")
    ax.plot(T, P_Zr, "-s", color=COLOR_ZR, lw=2, ms=4, label="Sitios Zr")
    ax.set_title(f"x = {zr}")
    ax.set_xlabel("T (K)")
    estilo(ax)
axes[0].set_ylabel(r"$P_L$")
axes[0].legend(frameon=False)
fig.suptitle("Polarización local promedio por tipo de sitio")
fig.tight_layout()
fig.savefig(OUT_IMAGE_TI_ZR, dpi=200)

plt.show()
