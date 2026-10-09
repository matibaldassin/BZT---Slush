import os
import numpy as np
import matplotlib.pyplot as plt

Zr = "0.40"
T_MIN_CURIE = 200   # K, se ajusta Curie-Weiss para T >= T_MIN_CURIE

base = os.path.dirname(os.path.abspath(__file__))
dir_tablas = os.path.join(base, "..", "Datos-corridas", "Permitividad_Dieléctrica", "Tablas_Permitividad")
tabla = os.path.join(dir_tablas, f"Permitividad_Zr{Zr}.txt")

# columnas de la tabla: T, dP/dE, eps_r, err_eps, ...
d = np.loadtxt(tabla, skiprows=1, usecols=(0, 2, 3))
T, eps, err_eps = d.T

inv_eps = 1 / eps
err_inv = err_eps / eps ** 2     # propagación: sigma(1/eps) = sigma(eps)/eps^2

# ---------------- Ajuste Curie-Weiss a altas T: 1/eps = (T - T0)/C ----------------
# ajuste sin pesos: con solo 3 puntos por T, los errores individuales no son fiables
m = T >= T_MIN_CURIE
(a1, a0), cov = np.polyfit(T[m], inv_eps[m], 1, cov=True)
C = 1 / a1                     # constante de Curie (K)
T0 = -a0 / a1                  # temperatura de Curie-Weiss (K)
err_a1, err_a0 = np.sqrt(np.diag(cov))
err_C = err_a1 / a1 ** 2
err_T0 = T0 * np.hypot(err_a0 / a0, err_a1 / a1)
res = inv_eps[m] - (a1 * T[m] + a0)
r2 = 1 - np.sum(res ** 2) / np.sum((inv_eps[m] - inv_eps[m].mean()) ** 2)

print(f"Ajuste Curie-Weiss (T >= {T_MIN_CURIE} K, {m.sum()} puntos):")
print(f"  C  = {C:.1f} +/- {err_C:.1f} K")
print(f"  T0 = {T0:.1f} +/- {err_T0:.1f} K")
print(f"  R2 = {r2:.4f}")

with open(os.path.join(dir_tablas, f"CurieWeiss_Zr{Zr}.txt"), "w") as f:
    f.write(f"# 1/eps = (T - T0)/C, ajustado para T >= {T_MIN_CURIE} K\n")
    f.write(f"C  = {C:.2f} +/- {err_C:.2f} K\nT0 = {T0:.2f} +/- {err_T0:.2f} K\nR2 = {r2:.4f}\n")

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.errorbar(T[~m], inv_eps[~m], yerr=err_inv[~m], fmt="o", color="gray", capsize=3)
ax.errorbar(T[m], inv_eps[m], yerr=err_inv[m], fmt="o", color="C0", capsize=3,
            label=f"$T \\geq$ {T_MIN_CURIE} K")
Tf = np.linspace(min(T0, T.min()) if T0 < T.min() else T.min(), T.max(), 100)
ax.plot(Tf, a1 * Tf + a0, "-", color="C3",
        label=f"Curie-Weiss: $C$={C:.0f} K, $T_0$={T0:.0f} K, $R^2$={r2:.3f}")
ax.set_ylim(bottom=min(0, inv_eps.min() - err_inv.max()))
ax.set_xlabel("T (K)")
ax.set_ylabel(r"$1/\varepsilon_r$")
ax.set_title(f"Inversa de la permitividad, Zr = {Zr}")
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(base, f"InvPermitividad_vs_T_Zr{Zr}.png"), dpi=200)
plt.close(fig)
print("Gráfico guardado.")
