#!/usr/bin/env python3
"""
Reescritura en Python de tabpol_B.f90.

Lee un CONFIG (formato DL_POLY-like) y para cada catión de sitio B
(Mg, Zr o Ti) arma la lista de sus 29 vecinos (A-site + O + su propio
shell) dentro de radios de corte por especie, igual que el original.
Escribe tabpol.dat con el mismo formato de columnas.

A diferencia del Fortran original:
 - nat (total de atomos) y ncel (numero de centros) se detectan
   solos leyendo el CONFIG, no hay que tipearlos a mano.
 - la busqueda de vecinos usa un cKDTree periodico en vez de recorrer
   todos los atomos por cada centro (O(N log N) en vez de O(N_centros*N)).
"""

import argparse
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

# especie -> cutoff^2 (mismos valores que el elif de tabpol_B.f90;
# 'Pb-sh' aparece dos veces en el original con el mismo valor, no cambia nada)
CUTOFF2 = {
    'MgA': 19.765, 'MgA-sh': 19.765, 'Pbs-sh': 19.765,
    'Ca': 19.765, 'Ca-sh': 19.765,
    'Sr': 19.765, 'Sr-sh': 19.765,
    'Ba': 19.765, 'Ba-sh': 19.765,
    'Pb': 19.765, 'Pb-sh': 19.765,
    'O': 11.0, 'O-sh': 11.0,
    'Os': 12.0, 'Os-sh': 12.0,
    'Mg-sh': 2.0, 'Zr-sh': 2.0, 'Ti-sh': 2.0,
    'Fe-sh': 2.0,
}

# especie -> peso (igual que el bloque lwei del original)
WEIGHT = {
    'Ba': 8, 'Ba-sh': 8,
    'Pbs': 8, 'Pbs-sh': 8,
    'Pb': 8, 'Pb-sh': 8,
    'Bi': 8, 'Bi-sh': 8,
    'Ca': 8, 'Ca-sh': 8,
    'O': 2, 'O-sh': 2,
    'Os': 2, 'Os-sh': 2,
    'Mg-sh': 1, 'Fe-sh': 1, 'Ti-sh': 1,
}
WEIGHT_DEFAULT = 100

CENTER_SPECIES = {'Mg', 'Zr', 'Ti'}
N_NEIGHBORS = 29


def read_config(path):
    with open(path) as f:
        lines = f.readlines()

    title = lines[0].rstrip('\n')
    lconf = int(lines[1].split()[0])

    cell = np.array(
        [[float(x) for x in lines[2 + j].split()] for j in range(3)]
    )

    names = []
    coords = []
    idx = 5
    block = 3 if lconf > 0 else 2
    while idx < len(lines) and lines[idx].strip():
        name = lines[idx][:10].strip()
        x, y, z = (float(v) for v in lines[idx + 1].split()[:3])
        names.append(name)
        coords.append((x, y, z))
        idx += block

    return title, lconf, cell, names, np.array(coords)


def build_periodic_tree(cell, coords):
    """Devuelve (tree, coords_para_indices, es_diagonal).

    Si la celda es ortogonal usa el soporte periodico nativo de
    cKDTree (rapido, sin replicar atomos). Si es triclinica, replica
    27 imagenes y devuelve un mapa de indices de vuelta a los originales.
    """
    diag = np.diag(np.diag(cell))
    if np.allclose(cell, diag):
        L = np.diag(cell)
        wrapped = coords % L
        tree = cKDTree(wrapped, boxsize=L)
        return tree, None, L

    # celda triclinica: fuerza bruta con 27 imagenes periodicas
    inv = np.linalg.inv(cell)
    frac = coords @ inv
    frac -= np.floor(frac)

    shifts = np.array([[i, j, k] for i in (-1, 0, 1)
                        for j in (-1, 0, 1) for k in (-1, 0, 1)])
    rep_frac = (frac[None, :, :] + shifts[:, None, :]).reshape(-1, 3)
    rep_cart = rep_frac @ cell
    orig_index = np.tile(np.arange(len(coords)), len(shifts))

    tree = cKDTree(rep_cart)
    return tree, orig_index, None


def minimum_image(dx, cell, L):
    if L is not None:
        return dx - L * np.round(dx / L)
    inv = np.linalg.inv(cell)
    frac = dx @ inv
    frac -= np.round(frac)
    return frac @ cell


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('config', nargs='?', default='CONFIG')
    ap.add_argument('out', nargs='?', default='tabpol.dat')
    args = ap.parse_args()

    t0 = time.time()
    title, lconf, cell, names, coords = read_config(args.config)
    n_atoms = len(names)
    print(f"{title}  ({n_atoms} atomos leidos, lconf={lconf})")

    centers = [k for k, nm in enumerate(names) if nm in CENTER_SPECIES]
    ncel = len(centers)
    print(f"{ncel} centros (Mg/Zr/Ti) detectados en {time.time()-t0:.1f}s")

    tree, orig_index, L = build_periodic_tree(cell, coords)
    max_r = max(CUTOFF2.values()) ** 0.5

    out_lines = [f"{ncel:>12d}  Centros"]
    n_warn = 0

    for k in centers:
        pos = coords[k]
        cand = tree.query_ball_point(pos, r=max_r)
        neighbors = []
        for c in cand:
            kp = c if orig_index is None else orig_index[c]
            if kp == k:
                continue
            sp = names[kp]
            cutoff2 = CUTOFF2.get(sp)
            if not cutoff2:
                continue
            dx = coords[kp] - pos
            dx = minimum_image(dx, cell, L)
            rsq = dx @ dx
            if rsq < cutoff2:
                neighbors.append(kp)

        neighbors.sort()

        if len(neighbors) != N_NEIGHBORS:
            n_warn += 1
            print(f"error vecinos {k + 1} {len(neighbors)}", file=sys.stderr)

        out_lines.append(f"{k + 1:8d}{len(neighbors):5d}   {names[k]:<10s}")
        row = neighbors[:N_NEIGHBORS] + [0] * max(0, N_NEIGHBORS - len(neighbors))
        out_lines.append("".join(f"{(n + 1):7d}" for n in row))

    weights = [WEIGHT.get(nm, WEIGHT_DEFAULT) for nm in names]
    for i in range(0, len(weights), 20):
        chunk = weights[i:i + 20]
        out_lines.append("".join(f"{w:5d}" for w in chunk))

    with open(args.out, 'w') as f:
        f.write("\n".join(out_lines) + "\n")

    print(f"listo en {time.time()-t0:.1f}s, {n_warn} centros con vecinos != {N_NEIGHBORS}")


if __name__ == '__main__':
    main()
