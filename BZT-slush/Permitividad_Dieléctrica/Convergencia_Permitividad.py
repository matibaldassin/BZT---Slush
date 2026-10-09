import os
import re
import numpy as np
import matplotlib.pyplot as plt

# ---------------- Parámetros ----------------
Zr = "0.50"
CAMPOS = [0, 5, 10]          # kV/mm
EPS0 = 8.8541878128e-12      # F/m
UC_A_SI = 1e-2               # 1 uC/cm^2 -> C/m^2
KV_MM_A_SI = 1e6             # 1 kV/mm -> V/m
N_SUAVIZADO = 3              # ventana del suavizado (promedio móvil); 1 = sin suavizar

base = os.path.dirname(os.path.abspath(__file__))
dir_ptot = os.path.join(base, "..", "Datos-corridas", "Permitividad_Dieléctrica",
                        "P-totales", f"P-totales-Zr{Zr}")


def leer_pz(T, E):
    """Serie temporal de Pz (uC/cm^2). Formato: 2 líneas de cabecera + paso Px Py Pz."""
    ruta = os.path.join(dir_ptot, f"{E}kV_mm", f"P-total-{T}")
    return np.loadtxt(ruta, skiprows=2)[:, 3]


def suavizar(y, n=N_SUAVIZADO):
    """Promedio móvil centrado de ventana n; en los bordes usa los puntos disponibles."""
    y = np.asarray(y, dtype=float)
    k = np.ones(n)
    return np.convolve(y, k, mode="same") / np.convolve(np.ones(len(y)), k, mode="same")


def permitividad(Pz_medio):
    """eps_r = 1 + (dP/dE)/eps0 con ajuste lineal de <Pz> vs E."""
    s, _ = np.polyfit(np.array(CAMPOS) * KV_MM_A_SI, np.array(Pz_medio) * UC_A_SI, 1)
    return 1 + s / EPS0


# ---------------- Cálculo ----------------
temps = sorted(int(m.group(1)) for fn in os.listdir(os.path.join(dir_ptot, "0kV_mm"))
               if (m := re.fullmatch(r"P-total-(\d+)", fn)))

eps_todo, eps_mitad, eps_cuarto = [], [], []
eps_inst = []   # eps_inst[i][k]: permitividad de T[i] en el paso k
for T in temps:
    series = [leer_pz(T, E) for E in CAMPOS]
    n = min(len(s) for s in series)
    eps_inst.append([permitividad([s[k] for s in series]) for k in range(n)])
    eps_todo.append(permitividad([s.mean() for s in series]))
    eps_mitad.append(permitividad([s[len(s) // 2:].mean() for s in series]))
    eps_cuarto.append(permitividad([s[3 * len(s) // 4:].mean() for s in series]))

# ---------------- Gráfico ----------------
fig, ax = plt.subplots(figsize=(6.5, 4.5))
ax.plot(temps, suavizar(eps_todo), "o-", color="C0", label="Promedio total")
ax.plot(temps, suavizar(eps_mitad), "s--", color="C3", label="Promedio última mitad")
ax.plot(temps, suavizar(eps_cuarto), "^:", color="C2", label="Promedio último cuarto")
ax.set_xlabel("T (K)")
ax.set_ylabel(r"$\varepsilon_r$")
ax.set_title(f"Convergencia de la permitividad, Zr = {Zr}")
ax.grid(alpha=0.3)
ax.legend()
fig.tight_layout()
salida = os.path.join(base, f"Convergencia_Permitividad_Zr{Zr}.png")
fig.savefig(salida, dpi=200)
print(f"Gráfico guardado en {salida}")

# ---------------- Gráfico instantáneo (sin promediar) ----------------
n_pasos = min(len(e) for e in eps_inst)
eps_inst = np.array([e[:n_pasos] for e in eps_inst])   # (T, paso)
fig2, ax2 = plt.subplots(figsize=(6.5, 4.5))
for k in range(n_pasos):
    alpha = 0.05 + 0.95 * k / max(n_pasos - 1, 1)      # claro = inicial, fuerte = final
    ax2.plot(temps, suavizar(eps_inst[:, k]), "-", color="C0", alpha=alpha, lw=1)
ax2.set_xlabel("T (K)")
ax2.set_ylabel(r"$\varepsilon_r$ instantánea")
ax2.set_title(f"Permitividad paso a paso, Zr = {Zr}")
ax2.grid(alpha=0.3)
fig2.tight_layout()
salida2 = os.path.join(base, f"Permitividad_instantanea_Zr{Zr}.png")
fig2.savefig(salida2, dpi=200)
print(f"Gráfico guardado en {salida2}")
# ---------------- Variantes de color para la convergencia ----------------
pasos = np.loadtxt(os.path.join(dir_ptot, f"{CAMPOS[0]}kV_mm", f"P-total-{temps[0]}"),
                   skiprows=2)[:n_pasos, 0]
T_arr = np.array(temps)

# 1) Colormap secuencial con colorbar (color = paso)
fig3, ax3 = plt.subplots(figsize=(7, 4.5))
cmap = plt.cm.viridis
norm = plt.Normalize(pasos[0], pasos[-1])
for k in range(n_pasos):
    ax3.plot(T_arr, suavizar(eps_inst[:, k]), "-", color=cmap(norm(pasos[k])), lw=1)
fig3.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax3, label="Paso")
ax3.set_xlabel("T (K)")
ax3.set_ylabel(r"$\varepsilon_r$ instantánea")
ax3.set_title(f"Permitividad paso a paso (colormap), Zr = {Zr}")
ax3.grid(alpha=0.3)
fig3.tight_layout()
fig3.savefig(os.path.join(base, f"Permitividad_colormap_Zr{Zr}.png"), dpi=200)

# 2) Bloques temporales: promedio de cada quinto de la serie
n_bloques = 5
fig4, ax4 = plt.subplots(figsize=(6.5, 4.5))
colores = plt.cm.viridis(np.linspace(0, 0.9, n_bloques))
for b, idx in enumerate(np.array_split(np.arange(n_pasos), n_bloques)):
    ax4.plot(T_arr, suavizar(eps_inst[:, idx].mean(axis=1)), "o-", color=colores[b], ms=4,
             label=f"Pasos {pasos[idx[0]]:.0f}-{pasos[idx[-1]]:.0f}")
ax4.set_xlabel("T (K)")
ax4.set_ylabel(r"$\varepsilon_r$")
ax4.set_title(f"Promedio por bloques temporales, Zr = {Zr}")
ax4.grid(alpha=0.3)
ax4.legend(fontsize=8)
fig4.tight_layout()
fig4.savefig(os.path.join(base, f"Permitividad_bloques_Zr{Zr}.png"), dpi=200)

# 3) Promedio acumulado vs paso, un color por temperatura
fig5, ax5 = plt.subplots(figsize=(7, 4.5))
cmap_T = plt.cm.plasma
norm_T = plt.Normalize(T_arr.min(), T_arr.max())
acum = np.cumsum(eps_inst, axis=1) / np.arange(1, n_pasos + 1)
for i, T in enumerate(T_arr):
    ax5.plot(pasos, suavizar(acum[i]), "-", color=cmap_T(norm_T(T)), lw=1.5)
fig5.colorbar(plt.cm.ScalarMappable(norm=norm_T, cmap=cmap_T), ax=ax5, label="T (K)")
ax5.set_xlabel("Paso")
ax5.set_ylabel(r"$\varepsilon_r$ (promedio acumulado)")
ax5.set_title(f"Promedio acumulado vs paso, Zr = {Zr}")
ax5.grid(alpha=0.3)
fig5.tight_layout()
fig5.savefig(os.path.join(base, f"Permitividad_acumulado_Zr{Zr}.png"), dpi=200)

plt.show()
