# ------- SCRIPT PARAMETERS ------
# Input:
#IN_FOLDER = '../Grouped_samples/'
#IN_FILES = ['grouped_' + str(5*i).zfill(2) + '.dump' for i in range(20)]
IN_FOLDER = '/home/matias/Escritorio/CONICET/bzt/Clusters_BZT/0.40-S100/grouped_files_cutoff_10/'
IN_FILES = ['grouped_ss_' + str(i) + '.dump' for i in range(50)]
# Output:
OUT_FOLDER = '/home/matias/Escritorio/CONICET/bzt/Clusters_BZT/0.40-S100/Phase_clustered_samples/'

# Group mapping:
GROUP_MAP = {0: 'T', 1: 'R', 2: 'O', 3: 'C', 4: 'Zr'}

# Column for groups
G_COL = 'Group'

# Columns for dipole 
PX_COL = 'DipoleOrientationX'
PY_COL = 'DipoleOrientationY'
PZ_COL = 'DipoleOrientationZ'

# Max radius for nearest neighbors
RADIUS = 5.

# Maximum atoms per slice
MAX_ATOMS = 40000


# ---------- LIBRARIES ----------
import os
import time
from collections import deque
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

# ----------- ENGINE -------------
from engine import *

# Numeric codes for each group name
GROUP_MAP_INV = {v: k for k, v in GROUP_MAP.items()}
GROUPS_T, GROUPS_R, GROUPS_O, GROUPS_C, GROUPS_ZR = (
    GROUP_MAP_INV['T'], GROUP_MAP_INV['R'], GROUP_MAP_INV['O'],
    GROUP_MAP_INV['C'], GROUP_MAP_INV['Zr'])


# ----------- SCRIPT ------------
# Fix random seed
np.random.seed(1235)

# Make sure the output folder exists
os.makedirs(OUT_FOLDER, exist_ok=True)

# Begin timing
t1 = time.time()

# Operate over each file
for i, f in enumerate(IN_FILES):
    print('\nProcessing file %d of %d' % (i+1, len(IN_FILES)))

    # Fetch file with data
    fname = IN_FOLDER + f
    df, norm_df, header_lines = open_and_extract(fname, [])

    # Assign IDs to the atoms
    n_atoms = df.shape[0]
    all_ids = np.arange(n_atoms)
    df['Atom_id'] = all_ids

    # Initialize list of unexplored atoms
    avail_atoms = list(all_ids)

    # Positions and groups indexed directly by Atom_id
    # (Atom_id == row position)
    pos_arr = df[['x', 'y', 'z']].values
    group_arr = df[G_COL].values

    # Get periodic info and generate replicas
    periodic_info = get_periodic(header_lines)
    replicas = periodic_replicate(df, periodic_info[1])
    replica_ids = np.array(replicas['Atom_id'])

    # Neighbor calculation
    knn = NearestNeighbors(radius=RADIUS, n_jobs=-1)
    knn.fit(replicas[['x','y','z']])
    distances, indices = knn.radius_neighbors(replicas[['x','y','z']])

    # Store neighbors (list indexed by Atom_id)
    neigh_map = [None] * n_atoms
    for a in range(n_atoms):
        all_indices = indices[a][indices[a] != a]
        neigh_map[a] = replica_ids[all_indices]

    # Get smallest dipole component
    min_p_comp = np.min(np.abs(df[[PX_COL, PY_COL, PZ_COL]]))

    # Generate reduced representation of the dipoles, vectorized.
    # The random draws below are consumed in the same atom-major
    # order as the original per-atom loop, so the noise (and thus
    # the tie-breaking it causes) is identical for the same seed.
    dipole_arr = df[[PX_COL, PY_COL, PZ_COL]].values
    noise = np.random.uniform(0, min_p_comp / 100, size=dipole_arr.shape)
    noisy_dipole_arr = dipole_arr + noise * np.sign(dipole_arr)

    # Reproduce the exact tie-breaking rule of the original code:
    # max_ax/third_ax pick the FIRST index achieving the extreme value
    # (np.argmax/np.argmin semantics), and second_ax picks the first
    # index whose magnitude equals the sorted middle value. On an exact
    # tie this can select the same axis twice (a pre-existing quirk of
    # the original algorithm), which is intentionally preserved here.
    abs_noisy = np.abs(noisy_dipole_arr)
    max_ax = np.argmax(abs_noisy, axis=1)
    third_ax = np.argmin(abs_noisy, axis=1)
    sorted_p = np.sort(abs_noisy, axis=1)
    second_ax = np.argmax(abs_noisy == sorted_p[:, 1:2], axis=1)
    rows = np.arange(n_atoms)

    # Sign is unaffected by the noise (it only pushes away from zero)
    sign_arr = np.where(dipole_arr > 0, 1, -1)

    mask_tor = (group_arr == GROUPS_T) | (group_arr == GROUPS_O) | (group_arr == GROUPS_R)
    mask_or = (group_arr == GROUPS_O) | (group_arr == GROUPS_R)
    mask_r = group_arr == GROUPS_R

    dipole_map = np.zeros((n_atoms, 3), dtype=int)
    dipole_map[rows[mask_tor], max_ax[mask_tor]] += sign_arr[rows[mask_tor], max_ax[mask_tor]]
    dipole_map[rows[mask_or], second_ax[mask_or]] += sign_arr[rows[mask_or], second_ax[mask_or]]
    dipole_map[rows[mask_r], third_ax[mask_r]] += sign_arr[rows[mask_r], third_ax[mask_r]]

    # Initialize cluster assignments. Default
    # value for every atom is -1, which means it
    # does not belong to a cluster
    cluster_ids = np.full(n_atoms, -1, dtype=int)

    # Global id for current cluster
    current_cluster = 1

    # Begin cluster identification
    while len(avail_atoms) > 0:
        # Select a random seed atom for the new cluster
        current_atom = np.random.choice(avail_atoms)

        # Initialize cluster info
        cluster_atoms = [current_atom]
        cluster_atoms_set = {current_atom}
        unexplored_atoms = deque([current_atom])

        # Perform cluster growth
        while len(unexplored_atoms) > 0:
            # Update current atom
            current_atom = unexplored_atoms.popleft()
            current_pos = pos_arr[current_atom]
            current_group = group_arr[current_atom]

            # Find neighbors in the three axes
            current_neighs = [[],[],[]]
            for n in neigh_map[current_atom]:
                n_pos = pos_arr[n]
                n_dist = np.abs(n_pos - current_pos)

                # Assign axis of the neighbor
                n_axis = np.argmax(n_dist)
                current_neighs[n_axis].append(n)

            # Operate over axis
            for l in range(3):
                for n in current_neighs[l]:
                    # Only attempt link for dipoles in the same group
                    n_group = group_arr[n]
                    if n_group != current_group or n_group == GROUPS_C or n_group == GROUPS_ZR:
                        continue

                    # Check if neighbor is linked
                    is_linked = np.array_equal(dipole_map[n], dipole_map[current_atom])

                    # If linked, update cluster data
                    if is_linked and (n not in cluster_atoms_set):
                        unexplored_atoms.append(n)
                        cluster_atoms.append(n)
                        cluster_atoms_set.add(n)

        # Remove cluster atoms from available list
        for c_atom in cluster_atoms:
            avail_atoms.remove(c_atom)

        # If a cluster was indeed found,
        # update global cluster id and
        # mappings
        if len(cluster_atoms) > 1:
            cluster_ids[cluster_atoms] = current_cluster

            # Update cluster id
            current_cluster += 1

    # Save cluster assignments to the data
    df['phase_cluster_id'] = cluster_ids

    # Drop atom id
    df.drop(columns=['Atom_id'], inplace=True)

    # Save new file with the results assigned to each atom
    new_path = OUT_FOLDER + 'phase_clusters_' + f
    export_sample(new_path, df, header_lines, 'phase_cluster_id')

# Finish timing and print
t2 = time.time()
print("\nProcess completed.")
print(f"Elapsed time: {round((t2-t1)/60, 1)} minutes.")

# Exit program
exit() 
