import random
import sys
import re

def format_name_index(elem: str, idx: int) -> str:
    # element left-justified in cols 1–17, index right-justified in 18–20
    return f"{elem:<17}{idx:>3}\n"

def main():
    # 1) read entire file
    with open("CONFIG", "r") as f:
        lines = f.readlines()

    # 2) split header vs body
    header = lines[:5]
    body   = lines[5:]
    if len(body) % 20 != 0:
        print("ERROR: CONFIG body is not a multiple of 20 lines per cell.")
        sys.exit(1)

    num_cells = len(body) // 20
    blocks = [body[i*20:(i+1)*20] for i in range(num_cells)]

    # 3) ask for N
    try:
        N = int(input(f"Enter number of B-site cells to replace (0–{num_cells}): "))
    except ValueError:
        print("ERROR: Please enter an integer.")
        sys.exit(1)
    if not (0 <= N <= num_cells):
        print(f"ERROR: N must be between 0 and {num_cells}.")
        sys.exit(1)

    # 4) select random cells
    zr_indices = set(random.sample(range(num_cells), N))

    # 5) rename selected blocks in-place (element names only)
    for i in zr_indices:
        blk = blocks[i]
        # core        = line 1
        # shell       = line 11
        # replace leading token “Ti”→“Zr” and “Ti-sh”→“Zr-sh”
        blk[2]  = re.sub(r"^Ti\b", "Zr", blk[2])
        blk[12] = re.sub(r"^Ti-sh\b", "Zr-sh", blk[12])

    # 6) reorder: all Zr first, then Ti
    reordered = [blocks[i] for i in range(num_cells) if i in zr_indices] \
              + [blocks[i] for i in range(num_cells) if i not in zr_indices]

    # 7) write output with fresh indices
    with open("CONFIG-salida", "w") as out:
        out.writelines(header)
        new_idx = 1
        for blk in reordered:
            # each block: 10 particles × 2 lines each
            for p in range(10):
                # extract element name from the original name-index line
                elem = blk[2*p].split()[0]
                coords_line = blk[2*p + 1]
                out.write(format_name_index(elem, new_idx))
                out.write(coords_line)
                new_idx += 1

    print(f"Done: wrote CONFIG-salida with {new_idx-1} renumbered particles.")

if __name__ == "__main__":
    main()
