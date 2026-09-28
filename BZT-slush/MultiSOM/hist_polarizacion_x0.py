import os

import pandas as pd
import matplotlib.pyplot as plt

from engine import open_and_extract

X_INDEX = 9  # 0 -> 20K, 1 -> 40K, ... 9 -> 200K

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
IN_FILE = os.path.join(SCRIPT_DIR, 'SOM_files', 'Grouped',
                        'Grouped_SOM_ss_0_trained_s50_0.%d.dump' % X_INDEX)
OUT_PATH_MIN = os.path.join(SCRIPT_DIR,
                             'hist_compminnorm_x%d.png' % X_INDEX)
OUT_PATH_MID = os.path.join(SCRIPT_DIR,
                             'hist_compmidnorm_x%d.png' % X_INDEX)
TEMPERATURE_K = 20 + 20 * X_INDEX

GROUPS = {'T': 0, 'R': 1, 'O': 2, 'C': 3, 'Zr': 4}


def plot_hist(values, xlabel, title, out_path):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(values, bins=50, color='#f4a261',
            edgecolor='black', linewidth=0.3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Cantidad de celdas')
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    print(values.describe())
    print('Guardado en %s\n' % out_path)


def main():
    df, _, _ = open_and_extract(IN_FILE, [])
    polar = df.loc[df['Group'].isin([GROUPS['T'], GROUPS['R'], GROUPS['O']])]
    print('n celdas polares:', len(polar))

    plot_hist(polar['compminnorm'], 'CompMinNorm',
              'CompMinNorm de celdas polares (T/R/O) - x=%d (%dK)' %
              (X_INDEX, TEMPERATURE_K),
              OUT_PATH_MIN)
    plot_hist(polar['compmidnorm'], 'CompMidNorm',
              'CompMidNorm de celdas polares (T/R/O) - x=%d (%dK)' %
              (X_INDEX, TEMPERATURE_K),
              OUT_PATH_MID)


if __name__ == '__main__':
    main()
