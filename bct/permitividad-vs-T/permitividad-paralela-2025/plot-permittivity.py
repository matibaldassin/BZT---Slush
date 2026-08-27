import matplotlib.pyplot as plt
import numpy as np
import glob
import os

# Busca automáticamente todas las tablas generadas: permittivity-Ca*-table.dat
TABLE_GLOB = "permittivity-Ca*-table.dat"
PLOT_DIA = "permittivity-BCT-all.png"

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

    plt.figure(figsize=(9,6))
    plt.title("Permittivity vs Temperature", fontsize=16, fontweight='bold')
    plt.ylabel(r"$\varepsilon_{r,\parallel}$", fontsize=14, fontweight='bold')
    plt.xlabel("Temperature (K)", fontsize=14, fontweight='bold')
    plt.minorticks_on()

    # Ciclo de marcadores para distinguir series
    markers = ['o', 's', 'D', '^', 'v', '>', '<', 'p', 'X']

    for i, fpath in enumerate(files):
        temps, comps = read_table(fpath)
        if not temps:
            continue
        # elegir qué componente graficar: por defecto 'ez', pero Ba0.75 -> 'ex'
        label = label_from_filename(fpath)
        comp_key = 'ez'
        if label == 'Ca0.25' or '0.25' in label:
            comp_key = 'ex'

        values = comps.get({'ex':'ex','ey':'ey','ez':'ez'}[comp_key])
        # ordenar por temperatura por si acaso
        sorted_pairs = sorted(zip(temps, values), key=lambda x: x[0])
        temps_s, vals_s = zip(*sorted_pairs)

        marker = markers[i]
        plt.plot(temps_s, vals_s, marker=marker, label=f"{label}")

    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(PLOT_DIA, dpi=300, bbox_inches='tight')
    print(f"Figura guardada en: {PLOT_DIA}")
    plt.show()

if __name__ == '__main__':
    main()