import os
import re
import glob

import pandas as pd

from engine import open_and_extract

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GROUPED_FOLDER = os.path.join(SCRIPT_DIR, 'SOM_files', 'Grouped')
FILE_PATTERN = os.path.join(GROUPED_FOLDER, 'Grouped_SOM_ss_0_trained_s50_0.*.dump')
OUT_TABLE = os.path.join(GROUPED_FOLDER, 'thresholds_vs_temperature.csv')

GROUPS = {'T': 0, 'R': 1, 'O': 2, 'C': 3, 'Zr': 4}

# Indice de archivo -> temperatura: 0 -> 20K, 1 -> 40K, ... 19 -> 400K
T_MIN = 20
T_STEP = 20


def compute_thresholds(df):
    '''
    Determina los thresholds que efectivamente separan los grupos ya
    clasificados en df['Group'] (T/R/O/C/Zr), tomando el punto medio del
    hueco entre el grupo mas cercano de cada lado del corte:

      - DipoleMagnitude_threshold_polar: separa Cubica (no polar) de
        polar (T/O/R).
      - CompMinNorm_threshold_R: separa Rombohedrica (compminnorm alto)
        del resto de celdas polares (T u O, compminnorm bajo).
      - CompMidNorm_threshold_TO: dentro de las celdas no-R, separa
        Tetragonal (compmidnorm bajo) de Ortorrombica (compmidnorm alto).
    '''
    cubic = df.loc[df['Group'] == GROUPS['C']]
    tetra = df.loc[df['Group'] == GROUPS['T']]
    rhombo = df.loc[df['Group'] == GROUPS['R']]
    ortho = df.loc[df['Group'] == GROUPS['O']]

    polar = pd.concat([tetra, rhombo, ortho])
    not_rhombo = pd.concat([tetra, ortho])

    dipole_threshold = (cubic['DipoleMagnitude'].max() +
                         polar['DipoleMagnitude'].min()) / 2

    min_threshold = (rhombo['compminnorm'].min() +
                      not_rhombo['compminnorm'].max()) / 2

    mid_threshold = (tetra['compmidnorm'].max() +
                      ortho['compmidnorm'].min()) / 2

    return {
        'n_cubic': len(cubic),
        'n_tetra': len(tetra),
        'n_rhombo': len(rhombo),
        'n_ortho': len(ortho),
        'DipoleMagnitude_threshold_polar': dipole_threshold,
        'CompMinNorm_threshold_R': min_threshold,
        'CompMidNorm_threshold_TO': mid_threshold,
    }


def main():
    files = glob.glob(FILE_PATTERN)
    files = sorted(files, key=lambda f: int(
        re.search(r's50_0\.(\d+)\.dump', f).group(1)))

    rows = []
    for fname in files:
        idx = int(re.search(r's50_0\.(\d+)\.dump', fname).group(1))
        t = T_MIN + T_STEP * idx
        print('Procesando %s (T = %d K)...' % (os.path.basename(fname), t))

        df, _, _ = open_and_extract(fname, [])
        row = {'index': idx, 'temperature_K': t, 'file': os.path.basename(fname)}
        row.update(compute_thresholds(df))
        rows.append(row)

    table = pd.DataFrame(rows).sort_values('temperature_K').reset_index(drop=True)
    table.to_csv(OUT_TABLE, index=False)

    print()
    print(table.to_string(index=False))
    print('\nGuardado en %s' % OUT_TABLE)


if __name__ == '__main__':
    main()
