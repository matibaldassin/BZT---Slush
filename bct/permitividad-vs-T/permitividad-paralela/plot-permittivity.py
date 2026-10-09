import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import glob
import os

# Busca automáticamente todas las tablas generadas: permittivity-Ca*-table.dat
TABLE_GLOB = "permittivity-*-table.dat"
PLOT_DIA = "e33.png"

# Inset: (x0, y0, ancho, alto) en coordenadas relativas del eje principal
INSET_POSITION = (-0.12, 0.35, 0.6, 0.6)
INSET_XLIM = (260, 320)
INSET_YLIM = (40,200)
INSET_SHOW_GRID = True

mpl.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "axes.labelsize": 26,
        "xtick.labelsize": 22,
        "ytick.labelsize": 22,
        "axes.linewidth": 2.5,
        "lines.linewidth": 3.0,
        "lines.markersize": 8.0,
        "xtick.major.size": 7.0,
        "ytick.major.size": 7.0,
        "xtick.major.width": 3.5,
        "ytick.major.width": 3.5,
        "xtick.minor.size": 3.5,
        "ytick.minor.size": 3.5,
        "xtick.minor.width": 1.2,
        "ytick.minor.width": 1.2,
        "legend.fontsize": 18,
    }
)

def read_table(path):
    """Lee una tabla de permitividad generada por `permittivity.py`.

    Formato esperado (columnas): T exx eyy ezz slope_x slope_y slope_z
    Devuelve: temperatures (list[int]) y comps (dict) con claves 'ex','ey','ez' -> list[float]
    """
    temps = []
    ex = []
    ey = []
    ez = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                T = int(float(parts[0]))
                ex_v = float(parts[1])
                ey_v = float(parts[2])
                ez_v = float(parts[3])
            except ValueError:
                continue
            temps.append(1.45 * T)
            ex.append(ex_v)
            ey.append(ey_v)
            ez.append(ez_v)
    comps = {'ex': ex, 'ey': ey, 'ez': ez}
    return temps, comps

def label_from_filename(path):
    """Extrae una etiqueta legible desde el nombre de archivo.
    Ej: permittivity-Ba0.70-table.dat -> Ba0.70
    """
    base = os.path.basename(path)
    # asume formato: permittivity-<label>-table.dat
    parts = base.split('-')
    if len(parts) >= 3:
        return parts[1]
    return base

def main():
    files = sorted(glob.glob(TABLE_GLOB))
    if not files:
        print(f"No se encontraron archivos con patrón: {TABLE_GLOB}")
        return

    fig, ax = plt.subplots(figsize=(9,6))
    ax.set_ylabel(r"$\varepsilon_{33}$")
    ax.set_xlabel("Temperature (K)")
    ax.set_xlim(250, 510)
    ax.minorticks_on()

    # Ciclo de marcadores para distinguir series
    markers = ['o', 's', 'D', '^', 'v', '>', '<', 'p', 'X']

    for i, fpath in enumerate(files):
        temps, comps = read_table(fpath)
        if not temps:
            continue
        # elegir qué componente graficar: por defecto 'ez', pero Ca0.25 -> 'ex'
        label = label_from_filename(fpath)
        comp_key = 'ez'
        if label == 'Ca0.25' or '0.25' in label:
            comp_key = 'ex'

        values = comps.get({'ex':'ex','ey':'ey','ez':'ez'}[comp_key])
        # ordenar por temperatura por si acaso
        sorted_pairs = sorted(zip(temps, values), key=lambda x: x[0])
        temps_s, vals_s = zip(*sorted_pairs)

        marker = markers[i % len(markers)]
        ax.plot(temps_s, vals_s, marker=marker, label=f"{label}")

    ax.grid(True, linestyle='--', alpha=0.6, linewidth=1.2)
    ax.legend(loc='upper right', frameon=False)

    # Inset con el intervalo solicitado
    inset_width = f"{INSET_POSITION[2] * 100}%"
    inset_height = f"{INSET_POSITION[3] * 100}%"
    axins = inset_axes(
        ax,
        width=inset_width,
        height=inset_height,
        bbox_to_anchor=(INSET_POSITION[0], INSET_POSITION[1], INSET_POSITION[2], INSET_POSITION[3]),
        bbox_transform=ax.transAxes,
    )
    for line in ax.get_lines():
        axins.plot(
            line.get_xdata(),
            line.get_ydata(),
            marker=line.get_marker(),
            linestyle=line.get_linestyle(),
            color=line.get_color(),
            markersize=4,
        )

    axins.set_xlim(*INSET_XLIM)
    if INSET_YLIM is not None:
        axins.set_ylim(*INSET_YLIM)
    if INSET_SHOW_GRID:
        axins.grid(True, linestyle='--', alpha=0.4)
    axins.set_xticks([260, 280, 300, 320])
    # Establecer ticks en el eje Y cada 50 unidades
    ymin, ymax = axins.get_ylim()
    axins.set_yticks(np.arange(np.ceil(ymin / 50) * 50, np.floor(ymax / 50) * 50 + 50, 50))
    axins.tick_params(axis='both', which='both', labelsize=18)

    plt.tight_layout()
    fig.savefig(PLOT_DIA, dpi=300, bbox_inches='tight')
    print(f"Figura guardada en: {PLOT_DIA}")
    plt.show()

if __name__ == '__main__':
    main()