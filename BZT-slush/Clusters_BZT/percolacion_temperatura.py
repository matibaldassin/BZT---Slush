import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from engine import open_and_extract, get_periodic

# ARCHIVOS DE ENTRADA (salida de clasificacion.py, con columna 'Group' ya asignada)
IN_FOLDER = '/home/matias/Escritorio/CONICET/bzt/Clusters_BZT/grouped_files_cutoff_10/'
IN_FILES = [f'grouped_ss_{i}.dump' for i in range(50)]

# CARPETA / ARCHIVOS DE SALIDA
OUT_FOLDER = '/home/matias/Escritorio/CONICET/bzt/Clusters_BZT/Percolacion/'
OUT_FIG = OUT_FOLDER + 'percolacion_vs_temperatura.png'
OUT_CSV = OUT_FOLDER + 'percolacion_vs_temperatura.csv'

# Eje X: se usa el índice de archivo (0..49) como proxy de temperatura.
# Si en algún momento tenés la lista real de temperaturas (mismo orden que
# IN_FILES), reemplazá esta línea por, por ejemplo:
# EJE_X = [10,20,...,590]
EJE_X = list(range(len(IN_FILES)))
EJE_X_LABEL = 'Índice de archivo (proxy de T)'

# PARÁMETROS
RADIO = 5.0
GROUP_MAP = {0: 'T', 1: 'R', 2: 'O', 3: 'C', 4: 'Zr'}
G_COL = 'Group'
PX_COL, PY_COL, PZ_COL = 'DipoleOrientationX', 'DipoleOrientationY', 'DipoleOrientationZ'
AXES = ['x', 'y', 'z']


def reducir_dipolo(atom_g, atom_p):
    '''Reduce el dipolo continuo a componentes {-1,0,1} según simetría local.'''
    reduced_p = [0, 0, 0]
    if atom_g in ('T', 'O', 'R'):
        max_ax = np.argmax(np.abs(atom_p))
        reduced_p[max_ax] = 2 * int(atom_p[max_ax] > 0) - 1
        if atom_g in ('O', 'R'):
            sorted_p = np.sort(np.abs(atom_p))
            second_ax = np.argmax(np.abs(atom_p) == sorted_p[1])
            reduced_p[second_ax] = 2 * int(atom_p[second_ax] > 0) - 1
            if atom_g == 'R':
                third_ax = np.argmin(np.abs(atom_p))
                reduced_p[third_ax] = 2 * int(atom_p[third_ax] > 0) - 1
    return reduced_p


def replicar_con_shift(df, box_sizes):
    '''
    Igual que engine.periodic_replicate, pero además etiqueta cada fila
    replicada con su desplazamiento entero (en unidades de caja) por eje:
    0 (copia original), +1 o -1 (réplica desplazada).
    '''
    base = df.copy()
    base['shift_x'] = 0
    base['shift_y'] = 0
    base['shift_z'] = 0
    replicas = [base]

    for i, ax in enumerate(AXES):
        box_size = box_sizes[i]
        for sign in (1, -1):
            r = df.copy()
            r[ax] = df[ax] + sign * box_size
            r['shift_x'] = 0
            r['shift_y'] = 0
            r['shift_z'] = 0
            r[f'shift_{ax}'] = sign
            replicas.append(r)

    return pd.concat(replicas, ignore_index=True)


def analizar_percolacion(fname):
    '''
    Reconstruye el mismo grafo de conexiones que arma analisis_clusters.py
    (vecinos dentro de RADIO + criterio de alineación de dipolo reducido
    por eje), pero llevando la cuenta de un desplazamiento acumulado por
    átomo dentro de su cluster. Si un átomo es alcanzado dos veces con
    desplazamientos distintos, el cluster se conecta consigo mismo a
    través del borde periódico: percola en el/los eje(s) de la diferencia.

    Devuelve (N_total, lista_de_clusters), donde cada cluster es un dict
    {'atoms': [...], 'percola': set de ejes 'x'/'y'/'z'}.
    '''
    df, _, header_lines = open_and_extract(fname, [])
    periodic_axes, periodic_dims = get_periodic(header_lines)
    box_sizes = [max(d) - min(d) for d in periodic_dims]

    n_atoms = df.shape[0]
    df['GroupName'] = df[G_COL].map(GROUP_MAP)

    positions = df[[*AXES]].to_numpy()
    groups = df['GroupName'].to_numpy()
    dipolos = df[[PX_COL, PY_COL, PZ_COL]].to_numpy()

    dipole_map = np.array([reducir_dipolo(groups[i], dipolos[i]) for i in range(n_atoms)])

    # Réplicas periódicas con shift, y vecinos dentro del radio de corte
    replicas = replicar_con_shift(df, box_sizes)
    replica_ids = np.tile(np.arange(n_atoms), 7)
    replica_shifts = replicas[['shift_x', 'shift_y', 'shift_z']].to_numpy()

    knn = NearestNeighbors(radius=RADIO, n_jobs=-1)
    knn.fit(replicas[[*AXES]])
    _, indices = knn.radius_neighbors(df[[*AXES]])

    # BFS con seguimiento de desplazamiento acumulado (detección de percolación)
    visitado = -np.ones(n_atoms, dtype=int)
    clusters = []

    for atom in range(n_atoms):
        if visitado[atom] != -1:
            continue

        cluster_idx = len(clusters)
        clusters.append({'atoms': [atom], 'percola': set()})
        visitado[atom] = cluster_idx
        desplazamiento = {atom: np.zeros(3, dtype=int)}
        cola = [atom]

        while cola:
            actual = cola.pop(0)
            d_actual = desplazamiento[actual]
            pos_actual = positions[actual]

            for idx in indices[actual]:
                n_id = replica_ids[idx]
                if n_id == actual:
                    continue

                # El eje de enlace se determina con las posiciones reales
                # (sin desplazar), igual que en analisis_clusters.py: la
                # diferencia de coordenada real es enorme (~tamaño de caja)
                # en el eje que efectivamente envuelve la caja periódica,
                # así que sigue identificando el eje correcto.
                n_axis = np.argmax(np.abs(positions[n_id] - pos_actual))

                link_score = dipole_map[n_id][n_axis] * dipole_map[actual][n_axis]
                if link_score <= 0:
                    continue

                n_shift = replica_shifts[idx]
                propuesto = d_actual + n_shift

                if visitado[n_id] == -1:
                    visitado[n_id] = cluster_idx
                    desplazamiento[n_id] = propuesto
                    clusters[cluster_idx]['atoms'].append(n_id)
                    cola.append(n_id)
                elif visitado[n_id] == cluster_idx:
                    diff = propuesto - desplazamiento[n_id]
                    for eje, val in zip(AXES, diff):
                        if val != 0:
                            clusters[cluster_idx]['percola'].add(eje)

    return n_atoms, clusters


if __name__ == '__main__':
    os.makedirs(OUT_FOLDER, exist_ok=True)

    p_inf_list = []
    n_ejes_list = []
    filas_csv = []

    for i, fname in enumerate(IN_FILES):
        print(f'Procesando {fname}...')
        n_atoms, clusters = analizar_percolacion(IN_FOLDER + fname)

        mayor = max(clusters, key=lambda c: len(c['atoms']))
        tam_mayor = len(mayor['atoms'])
        p_inf = tam_mayor / n_atoms
        ejes = sorted(mayor['percola'])
        n_ejes = len(ejes)

        p_inf_list.append(p_inf)
        n_ejes_list.append(n_ejes)
        filas_csv.append({
            'archivo': fname,
            'indice': i,
            'N_atomos': n_atoms,
            'tam_cluster_mayor': tam_mayor,
            'P_inf': p_inf,
            'ejes_percolando': ','.join(ejes),
            'n_ejes_percolando': n_ejes,
        })
        print(f'  cluster mayor: {tam_mayor}/{n_atoms} (P_inf={p_inf:.3f}), percola en: {ejes or "ninguno"}')

    pd.DataFrame(filas_csv).to_csv(OUT_CSV, index=False)
    print(f'CSV guardado en {OUT_CSV}')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    ax1.plot(EJE_X, p_inf_list, marker='o', color='tab:red')
    ax1.set_ylabel(r'$P_\infty$ (cluster mayor / N)')
    ax1.set_ylim(0, 1.05)
    ax1.set_title('Percolación del cluster más grande')

    ax2.step(EJE_X, n_ejes_list, where='mid', marker='o', color='tab:blue')
    ax2.set_ylabel('Ejes que percolan (0-3)')
    ax2.set_yticks([0, 1, 2, 3])
    ax2.set_xlabel(EJE_X_LABEL)

    plt.tight_layout()
    plt.savefig(OUT_FIG, dpi=150)
    print(f'Figura guardada en {OUT_FIG}')
