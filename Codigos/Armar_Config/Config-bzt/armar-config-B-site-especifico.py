#!/usr/bin/env python3
"""
Script para reemplazar átomos Ti por Zr en posiciones específicas del archivo CONFIG.
Configura los parámetros abajo antes de ejecutar.
"""

import re
import sys

# ============================================================================
# CONFIGURACIÓN: Modifica estos parámetros
# ============================================================================

# Coordenadas de las células donde reemplazar Ti por Zr
# Ejemplo: [(40.0, 40.0, 40.0), (0.0, 0.0, 0.0), (80.0, 80.0, 80.0)]
# El grid es 20×20×20, celda 80×80×80
ZR_COORDS = [(40.0, 40.0, 40.0)]

# ============================================================================
# FIN DE CONFIGURACIÓN
# ============================================================================


def format_name_index(elem: str, idx: int) -> str:
    """Formatea línea de nombre-índice: elemento (cols 1-17) + índice (cols 18-20)"""
    return f"{elem:<17}{idx:>3}\n"


def coords_to_index(x, y, z, cell_x=80.0, cell_y=80.0, cell_z=80.0, 
                    nx=20, ny=20, nz=20):
    """
    Convierte coordenadas reales a índice lineal de célula.
    
    cell_x, cell_y, cell_z: dimensiones de la celda (en unidades físicas)
    nx, ny, nz: número de células en cada dimensión
    """
    # Calcular tamaño de cada célula
    cell_size_x = cell_x / nx
    cell_size_y = cell_y / ny
    cell_size_z = cell_z / nz
    
    # Convertir coordenadas reales a índices de célula (enteros)
    idx_x = int(x / cell_size_x)
    idx_y = int(y / cell_size_y)
    idx_z = int(z / cell_size_z)
    
    # Validar que los índices estén en rango
    idx_x = max(0, min(nx - 1, idx_x))
    idx_y = max(0, min(ny - 1, idx_y))
    idx_z = max(0, min(nz - 1, idx_z))
    
    # Convertir a índice lineal
    linear_idx = idx_x + idx_y * nx + idx_z * (nx * ny)
    
    return linear_idx


def linear_to_coords(idx, cell_x=80.0, cell_y=80.0, cell_z=80.0, 
                     nx=20, ny=20, nz=20):
    """Convierte índice lineal a coordenadas reales (centro de la célula)"""
    cell_size_x = cell_x / nx
    cell_size_y = cell_y / ny
    cell_size_z = cell_z / nz
    
    idx_z = idx // (nx * ny)
    remainder = idx % (nx * ny)
    idx_y = remainder // nx
    idx_x = remainder % nx
    
    x = (idx_x + 0.5) * cell_size_x
    y = (idx_y + 0.5) * cell_size_y
    z = (idx_z + 0.5) * cell_size_z
    
    return x, y, z



def read_config(filepath: str):
    """Lee el archivo CONFIG y devuelve (header, blocks)"""
    with open(filepath, "r") as f:
        lines = f.readlines()
    
    header = lines[:5]
    body = lines[5:]
    
    if len(body) % 20 != 0:
        print("ERROR: El cuerpo de CONFIG no es múltiplo de 20 líneas por célula.")
        sys.exit(1)
    
    num_cells = len(body) // 20
    blocks = [body[i*20:(i+1)*20] for i in range(num_cells)]
    
    return header, blocks, num_cells


def select_zr_indices(coords: list, cell_dims: tuple) -> set:
    """
    Convierte coordenadas reales a índices de células.
    Retorna el conjunto de índices a reemplazar con Zr.
    """
    cell_x, cell_y, cell_z = cell_dims
    zr_indices = set()
    for x, y, z in coords:
        idx = coords_to_index(x, y, z, cell_x, cell_y, cell_z)
        zr_indices.add(idx)
        # Mostrar conversión
        real_x, real_y, real_z = linear_to_coords(idx, cell_x, cell_y, cell_z)
        print(f"  Coord ({x}, {y}, {z}) → Índice {idx} (centro: {real_x:.1f}, {real_y:.1f}, {real_z:.1f})")
    return zr_indices



def main():
    # 1) Leer CONFIG
    header, blocks, num_cells = read_config("CONFIG")
    
    # Extraer dimensiones de la celda de las líneas del header
    # Líneas 2-4 contienen la matriz de la celda
    cell_line_x = [float(x) for x in header[2].split()]
    cell_line_y = [float(x) for x in header[3].split()]
    cell_line_z = [float(x) for x in header[4].split()]
    
    cell_x = cell_line_x[0]  # primera componente
    cell_y = cell_line_y[1]  # segunda componente
    cell_z = cell_line_z[2]  # tercera componente
    
    print(f"✓ Leído CONFIG: {num_cells} células encontradas")
    print(f"✓ Dimensiones de celda: ({cell_x}, {cell_y}, {cell_z})")
    print(f"✓ Coordenadas a reemplazar:")
    
    # 2) Seleccionar índices para Zr
    cell_dims = (cell_x, cell_y, cell_z)
    zr_indices = select_zr_indices(ZR_COORDS, cell_dims)
    
    print(f"✓ Posiciones a reemplazar: {sorted(zr_indices)}")
    
    # 3) Reemplazar Ti por Zr en los bloques seleccionados
    for i in zr_indices:
        blk = blocks[i]
        # Línea 2 (índice 2): core - cambiar "Ti" por "Zr"
        blk[2] = re.sub(r"^Ti\b", "Zr", blk[2])
        # Línea 12 (índice 12): shell - cambiar "Ti-sh" por "Zr-sh"
        blk[12] = re.sub(r"^Ti-sh\b", "Zr-sh", blk[12])
    
    # 4) Reordenar: todos los Zr primero, luego los Ti
    reordered = [blocks[i] for i in range(num_cells) if i in zr_indices] \
              + [blocks[i] for i in range(num_cells) if i not in zr_indices]
    
    # 5) Escribir CONFIG-salida con índices renumerados
    with open("CONFIG-salida", "w") as out:
        out.writelines(header)
        new_idx = 1
        for blk in reordered:
            # Cada bloque: 10 partículas × 2 líneas cada una
            for p in range(10):
                elem = blk[2*p].split()[0]
                coords_line = blk[2*p + 1]
                out.write(format_name_index(elem, new_idx))
                out.write(coords_line)
                new_idx += 1
    
    print(f"✓ Completado: CONFIG-salida escrito con {new_idx-1} partículas renumeradas")
    print(f"  - Zr: {len(zr_indices)} células")
    print(f"  - Ti: {num_cells - len(zr_indices)} células")


if __name__ == "__main__":
    main()
