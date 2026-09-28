import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from engine import open_and_extract, get_periodic

# ARCHIVO DE ENTRADA (debe tener la columna cluster_id ya calculada)
IN_FILE = '/home/matias/Escritorio/CONICET/bzt/Clusters_BZT/Cluster_sample_cutoff_10/grouped_ss_15.dump'

# CARPETA DE SALIDA
OUT_FOLDER = '/home/matias/Escritorio/CONICET/bzt/Clusters_BZT/Cortes_cutoff_10/grouped_ss_15/'

# PARÁMETROS
N_CELL = 28                # celdas por lado de la supercelda
AXIS = 'z'                 # eje de corte: 'x', 'y' o 'z'
DESIRED_ARROW_FRACTION = 2  # fracción de celda que debe ocupar la flecha más larga

G_COL = 'Group'
GROUP_MAP = {0: 'T', 1: 'R', 2: 'O', 3: 'C', 4: 'Zr'}
DIPOLE_COLS = {'x': 'DipoleOrientationX', 'y': 'DipoleOrientationY', 'z': 'DipoleOrientationZ'}
CLUSTER_COL = 'cluster_id'

# Ejes del plano de corte (horizontal, vertical) para cada eje de corte
PLANE_AXES = {'x': ('y', 'z'), 'y': ('x', 'z'), 'z': ('x', 'y')}

CMAP_PZ = plt.get_cmap('coolwarm')
CLUSTER_CMAP = plt.get_cmap('tab20')


def cluster_color(cid):
    if cid == -1:
        return (0.6, 0.6, 0.6, 1.0)
    return CLUSTER_CMAP(int(cid) % 20)


def cargar_datos(fname):
    df, _, header_lines = open_and_extract(fname, [])
    periodic_axes, periodic_dims = get_periodic(header_lines)

    for i, ax in enumerate('xyz'):
        dims = periodic_dims[i]
        box_min, box_size = min(dims), max(dims) - min(dims)
        a = box_size / N_CELL
        # Los sitios B (Ti/Zr) están desplazados medio parámetro de red
        # respecto al origen de la caja, por eso se usa floor y no round
        # (con round, los valores en X.5 caen en bordes de bin y generan
        # huecos/duplicados en la grilla por bankers' rounding)
        df[f'{ax}_idx'] = np.floor((df[ax] - box_min) / a).astype(int) % N_CELL

    df['GroupName'] = df[G_COL].map(GROUP_MAP)
    return df


def graficar_capa(df, axis, layer, norm, scale, out_path):
    axis_col = f'{axis}_idx'
    h_axis, v_axis = PLANE_AXES[axis]
    h_idx_col, v_idx_col = f'{h_axis}_idx', f'{v_axis}_idx'
    h_dip_col, v_dip_col = DIPOLE_COLS[h_axis], DIPOLE_COLS[v_axis]
    depth_dip_col = DIPOLE_COLS[axis]

    capa = df[df[axis_col] == layer]
    ti = capa[capa['GroupName'] != 'Zr']
    ti_polar = ti[ti['GroupName'] != 'C']
    ti_no_polar = ti[ti['GroupName'] == 'C']
    zr = capa[capa['GroupName'] == 'Zr']

    fig, ax = plt.subplots(figsize=(9, 8))

    # Fondo de celda de cada Ti, coloreado según su cluster_id
    for _, row in ti.iterrows():
        color = cluster_color(row[CLUSTER_COL])
        rect = mpatches.Rectangle((row[h_idx_col] - 0.5, row[v_idx_col] - 0.5), 1, 1,
                                   facecolor=color, edgecolor='none', alpha=0.85, zorder=1)
        ax.add_patch(rect)

    # Flechas de polarización (proyección en el plano de corte) para los Ti
    # polares (T/R/O). Los no polares (C) no llevan flecha, solo un punto.
    if len(ti_polar) > 0:
        ax.quiver(ti_polar[h_idx_col], ti_polar[v_idx_col], ti_polar[h_dip_col], ti_polar[v_dip_col],
                  color='black', angles='xy', scale_units='xy', scale=scale,
                  width=0.006, zorder=3)

    if len(ti_no_polar) > 0:
        ax.scatter(ti_no_polar[h_idx_col], ti_no_polar[v_idx_col], s=12,
                   color='dimgray', zorder=3)

    # Círculos para los Zr (sin flecha, solo posición)
    ax.scatter(zr[h_idx_col], zr[v_idx_col], s=180, facecolors='none',
               edgecolors='black', linewidths=1.3, zorder=2, label='Zr')

    ax.set_xlim(-0.5, N_CELL - 0.5)
    ax.set_ylim(-0.5, N_CELL - 0.5)
    ax.set_aspect('equal')
    ax.set_xlabel(f'{h_axis} (índice de celda)')
    ax.set_ylabel(f'{v_axis} (índice de celda)')
    ax.set_title(f'Corte {axis} = {layer}  (supercelda {N_CELL}x{N_CELL}x{N_CELL})')
    ax.legend(loc='upper right')

    # sm = plt.cm.ScalarMappable(cmap=CMAP_PZ, norm=norm)
    # sm.set_array([])
    # cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    # cbar.set_label(f'Dipolo fuera de plano ({axis})')

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    os.makedirs(OUT_FOLDER, exist_ok=True)

    df = cargar_datos(IN_FILE)

    # Normalización global (misma escala de color en todas las capas)
    depth_dip_col = DIPOLE_COLS[AXIS]
    vmax = df[depth_dip_col].abs().max()
    norm = Normalize(vmin=-vmax, vmax=vmax)

    # Escala global de flechas (mismo tamaño de flecha en todas las capas)
    h_axis, v_axis = PLANE_AXES[AXIS]
    ti_all = df[df['GroupName'] != 'Zr']
    max_mag = np.sqrt(ti_all[DIPOLE_COLS[h_axis]]**2 + ti_all[DIPOLE_COLS[v_axis]]**2).max()
    scale = max_mag / DESIRED_ARROW_FRACTION

    for layer in range(N_CELL):
        out_path = os.path.join(OUT_FOLDER, f'corte_{AXIS}_{layer:02d}.png')
        print(f'Generando corte {AXIS} = {layer} -> {out_path}')
        graficar_capa(df, AXIS, layer, norm, scale, out_path)

    print('Listo.')
