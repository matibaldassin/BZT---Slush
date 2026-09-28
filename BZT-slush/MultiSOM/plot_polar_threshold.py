import os

import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TABLE_PATH = os.path.join(SCRIPT_DIR, 'SOM_files', 'Grouped',
                           'thresholds_vs_temperature.csv')
OUT_PATH = os.path.join(SCRIPT_DIR, 'threshold_polar_no_polar.png')

NO_POLAR_COLOR = '#a8dadc'
POLAR_COLOR = '#f4a261'

T_COLOR = '#8ecae6'
O_COLOR = '#ffb703'
R_COLOR = '#e63946'

OUT_PATH_MIN = os.path.join(SCRIPT_DIR, 'threshold_R_vs_TO.png')
OUT_PATH_MID = os.path.join(SCRIPT_DIR, 'threshold_T_vs_O.png')


def plot_polar_no_polar(temps, df):
    thr = df['DipoleMagnitude_threshold_polar'].to_numpy()

    fig, ax = plt.subplots(figsize=(7, 5))
    y_top = thr.max() * 1.4
    ax.fill_between(temps, 0, thr, color=NO_POLAR_COLOR, label='No polar')
    ax.fill_between(temps, thr, y_top, color=POLAR_COLOR, label='Polar')
    ax.plot(temps, thr, color='black', marker='o', linewidth=1.5,
            label='Threshold')

    ax.set_xlabel('Temperatura (K)')
    ax.set_ylabel('DipoleMagnitude')
    ax.set_title('Frontera polar / no polar vs. temperatura (celdas de Ti)')
    ax.set_xlim(temps.min(), temps.max())
    ax.set_ylim(0, y_top)
    ax.legend(loc='upper right')
    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150)
    plt.close(fig)
    print('Guardado en %s' % OUT_PATH)


def plot_r_vs_to(temps, df):
    thr = df['CompMinNorm_threshold_R'].to_numpy()

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.fill_between(temps, 0, thr, color=T_COLOR, label='T o O')
    ax.fill_between(temps, thr, 1.05, color=R_COLOR, label='R (rombohedrica)')
    ax.plot(temps, thr, color='black', marker='o', linewidth=1.5)

    ax.set_xlabel('Temperatura (K)')
    ax.set_ylabel('CompMinNorm')
    ax.set_title('Umbral R vs. T/O (CompMinNorm) vs. temperatura')
    ax.set_xlim(temps.min(), temps.max())
    ax.set_ylim(0, 1.05)
    ax.legend(loc='upper right')
    fig.tight_layout()
    fig.savefig(OUT_PATH_MIN, dpi=150)
    plt.close(fig)
    print('Guardado en %s' % OUT_PATH_MIN)


def plot_t_vs_o(temps, df):
    thr = df['CompMidNorm_threshold_TO'].to_numpy()

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.fill_between(temps, 0, thr, color=T_COLOR, label='T (tetragonal)')
    ax.fill_between(temps, thr, 1.05, color=O_COLOR, label='O (ortorrombica)')
    ax.plot(temps, thr, color='black', marker='o', linewidth=1.5)

    ax.set_xlabel('Temperatura (K)')
    ax.set_ylabel('CompMidNorm')
    ax.set_title('Umbral T vs. O (CompMidNorm) vs. temperatura\n'
                  '(entre celdas no-R)')
    ax.set_xlim(temps.min(), temps.max())
    ax.set_ylim(0, 1.05)
    ax.legend(loc='upper right')
    fig.tight_layout()
    fig.savefig(OUT_PATH_MID, dpi=150)
    plt.close(fig)
    print('Guardado en %s' % OUT_PATH_MID)


def main():
    df = pd.read_csv(TABLE_PATH).sort_values('temperature_K')
    temps = df['temperature_K'].to_numpy()

    plot_polar_no_polar(temps, df)
    plot_r_vs_to(temps, df)
    plot_t_vs_o(temps, df)


if __name__ == '__main__':
    main()
