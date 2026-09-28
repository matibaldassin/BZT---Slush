import os
import numpy as np
import pandas as pd
from engine import open_and_extract, export_sample

# ARCHIVOS DE ENTRADA
IN_FOLDER = '/home/matias/Escritorio/CONICET/bzt/Fracción_de_fase/MultiSOM/ss_xx_files/0.40-S100/'
OUT_FOLDER = '/home/matias/Escritorio/CONICET/bzt/Clusters_BZT/0.40-S100/grouped_files_cutoff_10/'

IN_FILE = [f'ss_{i}.dump' for i in range(50)]

GROUPS = {'T': 0, 'R': 1, 'O': 2, 'C': 3, 'Zr': 4}

#THRESHOLD PARA CLASIFICACIÓN
polar_threshold = 10  # uC/cm^2     

def classify_cells(fname):

    # Leer archivo con nombres de columnas correctos
    df, norm_df, header_lines = open_and_extract(fname, [])
    df['Group'] = 0

    print("Columnas del archivo:")
    print(df.columns.tolist())
    print(f"\nForma del DataFrame: {df.shape}")
    print(f"\nPrimeras filas:")
    print(df.head())

    # CLASIFICACIÓN ENTRE CELDAS CUBICAS y CELDAS POLARES (VECTORIZADA)
    # Inicializar todas como Zr (tipo 2)
    df['Group'] = GROUPS['Zr']

    # Solo clasificar filas desde el índice 8790 en adelante
    mask_base = df.index >= df.index[8790]

    mask_cubic = mask_base & (df['DipoleMagnitude'] < polar_threshold)

    mask_rhombo = mask_base & \
                  (df['DipoleMagnitude'] >= polar_threshold) & \
                  (df['CompMidNorm'] > 0.35) & (df['CompMidNorm'] <= 1.0) & \
                  (df['CompMinNorm'] > 0.35) & (df['CompMinNorm'] <= 1.0)

    mask_ortho = mask_base & \
                 (df['DipoleMagnitude'] >= polar_threshold) & \
                 (df['CompMidNorm'] > 0.5) & (df['CompMidNorm'] <= 1.0) & \
                 (df['CompMinNorm'] < 0.35)

    mask_tetra = mask_base & \
                 (df['DipoleMagnitude'] >= polar_threshold) & \
                 (df['CompMidNorm'] <= 0.5) & (df['CompMinNorm'] < 0.35)

    df.loc[mask_cubic, 'Group'] = GROUPS['C']
    df.loc[mask_rhombo, 'Group'] = GROUPS['R'] 
    df.loc[mask_ortho, 'Group'] = GROUPS['O']
    df.loc[mask_tetra, 'Group'] = GROUPS['T']

    # Exportar
    out_name = OUT_FOLDER + 'grouped_' + fname.split('/')[-1]
    export_sample(out_name, df, header_lines, 'Group')

def main():
    os.makedirs(OUT_FOLDER, exist_ok=True)
    for fname in IN_FILE:
        print(f"Procesando archivo: {fname}")
        classify_cells(IN_FOLDER + fname)

if __name__ == "__main__":
    main()

        