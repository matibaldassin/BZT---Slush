import os
import numpy as np
import time
from collections import deque
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from engine import open_and_extract, get_periodic, periodic_replicate, export_sample

# ARCHIVOS DE ENTRADA
IN_FOLDER = '/home/matias/Escritorio/CONICET/BZT-slush/Clusters_BZT/Corridas_Nuevas/Zr0.60/grouped_files_cutoff_10/'
IN_FILES = ['grouped_ss-'+str(i)+'.dump' for i in range(60)]

# ARCHIVOS DE SALIDA
OUT_FOLDER = '/home/matias/Escritorio/CONICET/BZT-slush/Clusters_BZT/Corridas_Nuevas/Zr0.60/Cluster_sample_cutoff_10/'

# PARÁMETROS
RADIO = 5.0
GROUP_MAP = {0: 'T', 1: 'R', 2: 'O', 3: 'C', 4: 'Zr'}
GROUP_MAP_INV = {v: k for k, v in GROUP_MAP.items()}
GROUPS_T, GROUPS_R, GROUPS_O = GROUP_MAP_INV['T'], GROUP_MAP_INV['R'], GROUP_MAP_INV['O']
G_COL = 'Group'
PX_COL = 'DipoleOrientationX'
PY_COL = 'DipoleOrientationY'
PZ_COL = 'DipoleOrientationZ'

def contar_clusters(fname):
    """
    Cuenta clusters usando la misma lógica que polar_cluster_analysis.py
    """
    
    print(f"Leyendo archivo: {fname}")
    t1 = time.time()
    
    # Leer archivo (usa la función del engine)
    df, norm_df, header_lines = open_and_extract(fname, [])
    
    # Asignar IDs
    n_atoms = df.shape[0]
    all_ids = np.arange(n_atoms)
    df['Atom_id'] = all_ids
    avail_atoms = list(all_ids)

    # Posiciones indexadas directamente por Atom_id (Atom_id == posición de fila)
    pos_arr = df[['x', 'y', 'z']].values

    # Generar mapa de vecinos CON REPLICAS PERIÓDICAS
    periodic_info = get_periodic(header_lines)
    replicas = periodic_replicate(df, periodic_info[1])
    replica_ids = np.array(replicas['Atom_id'])

    print("Calculando vecinos...")
    knn = NearestNeighbors(radius=RADIO, n_jobs=-1)
    knn.fit(replicas[['x','y','z']])
    distances, indices = knn.radius_neighbors(replicas[['x','y','z']])

    # Guardar vecinos (lista indexada por Atom_id)
    neigh_map = [None] * n_atoms
    for i in range(n_atoms):
        all_indices = indices[i][indices[i] != i]
        neigh_map[i] = replica_ids[all_indices]

    # Generar mapa de dipolos reducidos (según simetría), vectorizado
    print("Generando dipolos reducidos...")
    group_arr = df[G_COL].values
    dipole_arr = df[[PX_COL, PY_COL, PZ_COL]].values
    sign_arr = np.where(dipole_arr > 0, 1, -1)

    # order[:,0]/[:,1]/[:,2] = eje de menor/media/mayor magnitud por fila
    order = np.argsort(np.abs(dipole_arr), axis=1)
    min_ax, mid_ax, max_ax = order[:, 0], order[:, 1], order[:, 2]
    rows = np.arange(n_atoms)

    T, R, O = GROUPS_T, GROUPS_R, GROUPS_O
    mask_tor = (group_arr == T) | (group_arr == O) | (group_arr == R)
    mask_or = (group_arr == O) | (group_arr == R)
    mask_r = group_arr == R

    dipole_map = np.zeros((n_atoms, 3), dtype=int)
    dipole_map[rows[mask_tor], max_ax[mask_tor]] = sign_arr[rows[mask_tor], max_ax[mask_tor]]
    dipole_map[rows[mask_or], mid_ax[mask_or]] = sign_arr[rows[mask_or], mid_ax[mask_or]]
    dipole_map[rows[mask_r], min_ax[mask_r]] = sign_arr[rows[mask_r], min_ax[mask_r]]

    # Clustering
    print("Identificando clusters...")
    cluster_ids = np.full(n_atoms, -1, dtype=int)
    current_cluster = 1

    while len(avail_atoms) > 0:
        # Semilla aleatoria
        current_atom = np.random.choice(avail_atoms)

        cluster_atoms = [current_atom]
        cluster_atoms_set = {current_atom}
        unexplored_atoms = deque([current_atom])

        # BFS
        while len(unexplored_atoms) > 0:
            current_atom = unexplored_atoms.popleft()
            current_pos = pos_arr[current_atom]

            # Organizar vecinos por eje
            current_neighs = [[],[],[]]
            for n in neigh_map[current_atom]:
                n_pos = pos_arr[n]
                n_dist = np.abs(n_pos - current_pos)
                n_axis = np.argmax(n_dist)
                current_neighs[n_axis].append(n)

            # Verificar alineamiento
            for l in range(3):
                for n in current_neighs[l]:
                    link_score = dipole_map[n][l] * dipole_map[current_atom][l]
                    is_linked = link_score > 0
                    if is_linked and n not in cluster_atoms_set:
                        unexplored_atoms.append(n)
                        cluster_atoms.append(n)
                        cluster_atoms_set.add(n)

        # Remover del list
        for c_atom in cluster_atoms:
            avail_atoms.remove(c_atom)

        # Solo contar si > 1 átomo
        if len(cluster_atoms) > 1:
            cluster_ids[cluster_atoms] = current_cluster
            current_cluster += 1

    num_clusters = current_cluster - 1

    t2 = time.time()
    print(f"Tiempo total: {(t2-t1)/60:.1f} minutos")
    print(f"Número de clusters: {num_clusters}")

    # Guardar cluster_id en dataframe
    df['cluster_id'] = cluster_ids

    # Asegurar que la carpeta exista
    os.makedirs(OUT_FOLDER, exist_ok=True)

    # Nombre de salida
    out_name = OUT_FOLDER + fname.split('/')[-1]

    # Exportar
    export_sample(out_name, df, header_lines, 'atom_id cluster_id')
    


# USO
if __name__ == '__main__':
    np.random.seed(1234)  # Seed global
    num_tot_clusters = []
    ids_tot_clusters = []
    for i in range(len(IN_FILES)):
        contar_clusters(IN_FOLDER + IN_FILES[i])
