#!/usr/bin/env python3
import random
import sys
import re

LINES_PER_CELL = 20
NAME_LINE_POS_CORE = 2    # 0-based in 20-line block (line 3)
NAME_LINE_POS_SHELL = 12  # 0-based in 20-line block (line 13)
PARTICLES_PER_CELL = 10

# regex: <element> <index>
NAME_IDX_RE = re.compile(r"^\s*(\S+)\s+(\d+)\s*$")

def format_name_index(elem: str, idx: int) -> str:
    # element left-justified in cols 1–17, index right-justified in 18–20
    return f"{elem:<17}{idx:>3}\n"

def first_token(s: str) -> str:
    s = s.strip()
    return s.split()[0] if s else ""

def parse_elem_and_index(name_line: str):
    """
    From a name/index line like 'Ti               1234',
    return ('Ti', 1234). Robust to spacing.
    """
    m = NAME_IDX_RE.match(name_line.rstrip("\n"))
    if not m:
        # fallback: split-based (in case of unexpected spacing)
        parts = name_line.split()
        if len(parts) >= 2 and parts[-1].isdigit():
            return parts[0], int(parts[-1])
        raise ValueError(f"Cannot parse element/index from line: {name_line!r}")
    return m.group(1), int(m.group(2))

def main():
    # 1) read entire file
    with open("CONFIG", "r") as f:
        lines = f.readlines()

    # 2) split header vs body
    if len(lines) < 5:
        print("ERROR: CONFIG seems too short (needs 5-line header).")
        sys.exit(1)

    header = lines[:5]
    body   = lines[5:]
    if len(body) % LINES_PER_CELL != 0:
        print("ERROR: CONFIG body is not a multiple of 20 lines per cell.")
        sys.exit(1)

    num_cells = len(body) // LINES_PER_CELL
    blocks = [body[i*LINES_PER_CELL:(i+1)*LINES_PER_CELL] for i in range(num_cells)]

    # 3) classify cells: existing Zr vs Ti
    zr_existing = []
    ti_cells = []
    for i, blk in enumerate(blocks):
        core_tok = first_token(blk[NAME_LINE_POS_CORE])
        shell_tok = first_token(blk[NAME_LINE_POS_SHELL])
        if core_tok not in ("Ti", "Zr") or shell_tok not in ("Ti-sh", "Zr-sh"):
            print(f"ERROR: Cell {i}: unexpected tokens: core='{core_tok}', shell='{shell_tok}'.")
            print("       Expected 'Ti'/'Zr' and 'Ti-sh'/'Zr-sh'.")
            sys.exit(1)
        if core_tok == "Zr":
            zr_existing.append(i)
        else:
            ti_cells.append(i)

    num_ti = len(ti_cells)

    # 4) ask for N (only among Ti cells)
    try:
        N = int(input(f"Enter number of B-site Ti cells to replace (0–{num_ti}): "))
    except ValueError:
        print("ERROR: Please enter an integer.")
        sys.exit(1)
    if not (0 <= N <= num_ti):
        print(f"ERROR: N must be between 0 and {num_ti}.")
        sys.exit(1)

    # 5) select random Ti cells and replace just those (leave existing Zr untouched)
    selected_ti_to_convert = set(random.sample(ti_cells, N)) if N > 0 else set()

    name_re_core  = re.compile(r"^\s*Ti\b")
    name_re_shell = re.compile(r"^\s*Ti-sh\b")

    # store actual INPUT CONFIG indices for the cores that we replace
    replaced_core_indices_input = []  # 1-based indices from input CONFIG

    for i in selected_ti_to_convert:
        blk = blocks[i]

        # Read the original element/index from the INPUT name line (core)
        elem0, idx0 = parse_elem_and_index(blk[NAME_LINE_POS_CORE])
        # Only record cores that are Ti (they should be)
        if elem0 == "Ti":
            replaced_core_indices_input.append(idx0)

        # convert core if it's Ti
        if name_re_core.search(blk[NAME_LINE_POS_CORE]):
            blk[NAME_LINE_POS_CORE] = name_re_core.sub("Zr", blk[NAME_LINE_POS_CORE], count=1)
        # convert shell if it's Ti-sh
        if name_re_shell.search(blk[NAME_LINE_POS_SHELL]):
            blk[NAME_LINE_POS_SHELL] = name_re_shell.sub("Zr-sh", blk[NAME_LINE_POS_SHELL], count=1)

    # 6) reorder: existing Zr first, then newly converted Zr, then remaining Ti
    converted_set = selected_ti_to_convert
    remaining_ti = [i for i in ti_cells if i not in converted_set]

    # keep deterministic output layout inside each group (optional)
    reordered_indices = zr_existing + sorted(converted_set) + remaining_ti
    reordered_blocks = [blocks[i] for i in reordered_indices]

    # 7) write output with fresh indices
    with open("CONFIG-salida", "w") as out:
        out.writelines(header)
        new_idx = 1
        for blk in reordered_blocks:
            # each block: 10 particles × 2 lines each
            for p in range(PARTICLES_PER_CELL):
                elem = first_token(blk[2*p])
                coords_line = blk[2*p + 1]
                out.write(format_name_index(elem, new_idx))
                out.write(coords_line)
                new_idx += 1

    # 8) write mapping file: ONLY the core indices (from INPUT CONFIG) that were replaced
    replaced_core_indices_input.sort()
    with open("salida-atomos-cambiados", "w") as m:
        m.write("# Random B-site substitutions Ti→Zr (only Ti cells considered)\n")
        m.write("# Input : CONFIG\n")
        m.write("# Output: CONFIG-salida\n")
        m.write(f"# Cells total          : {num_cells}\n")
        m.write(f"# Existing Zr cells    : {len(zr_existing)}\n")
        m.write(f"# Existing Ti cells    : {len(ti_cells)}\n")
        m.write(f"# Replaced Ti cores (N): {len(replaced_core_indices_input)}\n")
        if replaced_core_indices_input:
            m.write("# Replaced core indices in INPUT CONFIG (1-based), one per line:\n")
            for idx in replaced_core_indices_input:
                m.write(f"{idx}\n")

    print(f"Done: wrote CONFIG-salida with {new_idx-1} renumbered particles.")
    print(f"Order: [existing Zr: {len(zr_existing)}] + [new Zr: {len(converted_set)}] + [remaining Ti: {len(remaining_ti)}].")
    print("Mapping written to 'salida-atomos-cambiados' (INPUT core indices only).")

if __name__ == "__main__":
    main()
