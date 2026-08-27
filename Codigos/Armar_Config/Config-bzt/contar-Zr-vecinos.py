#!/usr/bin/env python3
import numpy as np
import sys
import matplotlib.pyplot as plt

# === User‐modifiable parameters ===
CUTOFF = 7.0  # Ångstroms; change this to any cutoff you like

def read_config(filename):
    """
    Reads a DL_POLY CONFIG file (positions‐only, 2 lines per atom).
    Returns:
      - atoms: list of (label_str, pos_array([x,y,z]))
      - cell:    3×3 numpy array of cell vectors (a, b, c)
      - pbc_flag: integer (0 = no PBC, 1 or 2 = orthorhombic/triclinic)
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    header_line = lines[1].split()
    try:
        pbc_flag = int(header_line[0])
    except (IndexError, ValueError):
        raise RuntimeError(f"Cannot parse PBC flag from line 2: {lines[1]!r}")

    try:
        cell = np.array([
            list(map(float, lines[2].split())),
            list(map(float, lines[3].split())),
            list(map(float, lines[4].split())),
        ])
    except Exception as e:
        raise RuntimeError(f"Error parsing cell vectors (lines 3–5): {e}")

    total_lines = len(lines)
    if total_lines < 6:
        raise RuntimeError("CONFIG file appears too short to contain any atoms.")
    natoms = (total_lines - 5) // 2

    atoms = []
    i = 5
    for _ in range(natoms):
        label = lines[i][:8].strip()
        try:
            pos = np.array(list(map(float, lines[i+1].split())))
        except Exception:
            raise RuntimeError(f"Malformed coordinate line at file line {i+1}: {lines[i+1]!r}")
        atoms.append((label, pos))
        i += 2

    return atoms, cell, pbc_flag


def count_zr_neighbors(atoms, cell, pbc_flag, cutoff):
    """
    Given:
      - atoms:    list of (label_str, pos_array)
      - cell:     3×3 array of cell vectors
      - pbc_flag: 0 or 1 or 2
      - cutoff:   distance in Å
    Returns:
      - zr_positions: N×3 array of all Zr positions
      - counts:       list of length N, where counts[i] = # of other Zr’s within cutoff of Zr i
    """
    zr_positions = [pos for (label, pos) in atoms if label == 'Zr']
    if len(zr_positions) == 0:
        return np.zeros((0,3)), []

    zr_positions = np.array(zr_positions)
    N = zr_positions.shape[0]
    counts = [0] * N

    cell_vecs   = cell
    cell_sqnorm = np.array([np.dot(cell_vecs[k], cell_vecs[k]) for k in range(3)])

    for i in range(N):
        pos_i = zr_positions[i]
        cnt = 0
        for j in range(N):
            if i == j:
                continue
            pos_j = zr_positions[j]
            delta = pos_j - pos_i

            if pbc_flag > 0:
                for k in range(3):
                    frac = np.dot(delta, cell_vecs[k]) / cell_sqnorm[k]
                    delta -= cell_vecs[k] * np.round(frac)

            dist = np.linalg.norm(delta)
            if dist < cutoff:
                cnt += 1

        counts[i] = cnt

    return zr_positions, counts


def main():
    config_file = "CONFIG"
    if len(sys.argv) > 1:
        config_file = sys.argv[1]

    try:
        atoms, cell, pbc_flag = read_config(config_file)
    except Exception as e:
        print(f"ERROR reading CONFIG: {e}", file=sys.stderr)
        sys.exit(1)

    # Use the global CUTOFF defined at top
    zr_positions, counts = count_zr_neighbors(atoms, cell, pbc_flag, cutoff=CUTOFF)
    total_zr = len(zr_positions)

    # 1) Write neighbor counts to vecinos-Zr.dat
    out_filename = "vecinos-Zr.dat"
    try:
        out_f = open(out_filename, 'w')
    except Exception as e:
        print(f"ERROR opening '{out_filename}' for writing: {e}", file=sys.stderr)
        sys.exit(1)

    out_f.write(f"# Zr neighbor counts (cutoff = {CUTOFF} Å)\n")
    out_f.write(f"# Total Zr atoms found: {total_zr}\n")
    out_f.write("# Format: Zr_index   neighbor_count\n")
    for idx, c in enumerate(counts, start=1):
        out_f.write(f"{idx:5d}   {c:3d}\n")
    out_f.close()
    print(f"Found {total_zr} Zr atoms total.")
    print(f"Wrote neighbor counts to '{out_filename}'.")

    # 2) Build and save the *normalized* histogram
    if total_zr == 0:
        print("No Zr atoms found; skipping histogram.")
        return

    counts_arr = np.array(counts)
    max_neighbors = counts_arr.max()
    bin_edges = np.arange(-0.5, max_neighbors + 1.5, 1.0)

    # Normalize so that the sum of bar heights = 1.0
    weights = np.ones_like(counts_arr, dtype=float) / total_zr

    plt.figure(figsize=(6,4))
    n, bins, patches = plt.hist(
        counts_arr,
        bins=bin_edges,
        weights=weights,
        edgecolor='black',
        linewidth=1.0
    )

    plt.xlabel("Number of Zr neighbors within {:.1f} Å".format(CUTOFF), fontweight='bold')
    plt.ylabel("Fraction of Zr atoms",                      fontweight='bold')
    plt.title("Normalized histogram of Zr–Zr neighbor counts", fontweight='bold')
    plt.xticks(np.arange(0, max_neighbors + 1, 1))

    plt.tight_layout()
    hist_filename = "histograma-Zr-normalizado.png"
    try:
        plt.savefig(hist_filename, dpi=300)
        plt.close()
        print(f"Saved normalized histogram to '{hist_filename}'.")
    except Exception as e:
        print(f"ERROR saving histogram: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
