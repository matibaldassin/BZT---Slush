#!/usr/bin/env python3
import numpy as np
import random

def read_config_lines(fname):
    """Read entire CONFIG into a list of lines."""
    with open(fname, 'r') as f:
        return f.readlines()

def write_config_lines(lines, fname):
    """Write modified lines back to disk."""
    with open(fname, 'w') as f:
        f.writelines(lines)

def get_lattice(lines):
    """Parse lines 2–4 (0-based) for the three cell-vector rows."""
    vecs = []
    for i in (2, 3, 4):
        vecs.append([float(x) for x in lines[i].split()])
    return np.array(vecs)

def min_image(d, box):
    """Apply minimum-image correction to displacement vector d."""
    return d - np.round(d / box) * box

def main():
    input_file   = 'CONFIG'
    output_file  = 'CONFIG-salida'
    zr_only_file = 'CONFIG-Zr'

    # --- load all lines ---
    lines = read_config_lines(input_file)
    atom_lines = lines[5:]
    if len(atom_lines) % 2 != 0:
        raise ValueError("Unexpected atom-section formatting")
    n_atoms  = len(atom_lines) // 2
    if n_atoms % 10 != 0:
        raise ValueError("CONFIG does not contain 10-atom perovskite blocks")
    n_blocks = n_atoms // 10

    # --- parse box for PBC ---
    lattice = get_lattice(lines)
    box     = np.array([lattice[0,0], lattice[1,1], lattice[2,2]])

    # --- extract Ti (B-site) core coords ---
    coords = np.zeros((n_blocks, 3))
    for b in range(n_blocks):
        species_idx = 5 + (b*10 + 1)*2
        coord_idx   = species_idx + 1
        coords[b]   = np.array([float(x) for x in lines[coord_idx].split()])

    # --- user inputs ---
    num_np = int(input(f"Number of BaZrO₃ nanoparticles (≤ {n_blocks}): "))
    if not (1 <= num_np <= n_blocks):
        raise ValueError("Invalid number of nanoparticles")
    radius = float(input("Sphere radius (Å): "))

    # --- choose random centers ---
    centers_idx = random.sample(range(n_blocks), num_np)
    centers     = coords[centers_idx]

    # --- find blocks to replace ---
    to_replace = set()
    for C in centers:
        for b in range(n_blocks):
            d = coords[b] - C
            d = min_image(d, box)
            if np.linalg.norm(d) <= radius:
                to_replace.add(b)
    print(f"Will replace {len(to_replace)} B-sites (plus their shells) out of {n_blocks}.")

    # --- perform Ti→Zr substitution ---
    for b in to_replace:
        core_line_idx  = 5 + (b*10 + 1)*2
        old_core_label = lines[core_line_idx].split()[0]
        lines[core_line_idx] = lines[core_line_idx].replace(old_core_label, 'Zr', 1)
        shell_line_idx  = 5 + (b*10 + 6)*2
        old_shell_label = lines[shell_line_idx].split()[0]
        lines[shell_line_idx] = lines[shell_line_idx].replace(old_shell_label, 'Zr-sh', 1)

    # --- build ordered CONFIG-salida: Zr-blocks first, then Ti-blocks, renumbered with proper alignment ---
    salida = []
    # header
    salida.extend(lines[:5])
    new_idx = 1
    # helper for aligned species line
    def make_species_line(species, idx):
        # original puts index at column ~20, so pad species accordingly
        pad = max(1, 19 - len(species))
        return f"{species}{' ' * pad}{idx}\n"

    # Zr-containing blocks
    for b in sorted(to_replace):
        for i in range(10):
            li = 5 + (b*10 + i)*2
            species = lines[li].split()[0]
            # write aligned species+new index
            salida.append(make_species_line(species, new_idx))
            # write coord line unchanged
            salida.append(lines[li+1])
            new_idx += 1
    # remaining Ti blocks
    for b in range(n_blocks):
        if b not in to_replace:
            for i in range(10):
                li = 5 + (b*10 + i)*2
                species = lines[li].split()[0]
                salida.append(make_species_line(species, new_idx))
                salida.append(lines[li+1])
                new_idx += 1

    write_config_lines(salida, output_file)
    print(f"Done: reordered and renumbered CONFIG written as '{output_file}'.")

    # --- write Zr-only CONFIG-Zr ---
    zr_lines = []
    zr_lines.extend(lines[:5])
    for b in sorted(to_replace):
        core_line_idx  = 5 + (b*10 + 1)*2
        coord_line_idx = core_line_idx + 1
        zr_lines.append(lines[core_line_idx])
        zr_lines.append(lines[coord_line_idx])
    write_config_lines(zr_lines, zr_only_file)
    print(f"Done: Zr-only CONFIG written as '{zr_only_file}'.")

if __name__ == '__main__':
    main()
