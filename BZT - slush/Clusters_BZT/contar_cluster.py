import pandas as pd
import matplotlib.pyplot as plt

def leer_dump_seleccion(path, columnas_interes=None):
    """
    Lee un archivo .dump de LAMMPS y devuelve un DataFrame con solo las columnas seleccionadas.
    Si hay varios timesteps, devuelve un diccionario {timestep: DataFrame}.
    """
    with open(path, 'r') as f:
        lines = f.readlines()

    data_all = {}
    timestep = None
    reading_atoms = False
    columnas = []

    for line in lines:
        line = line.strip()

        if line.startswith("ITEM: TIMESTEP"):
            reading_atoms = False
            timestep = None
            continue

        elif timestep is None and line.isdigit():
            timestep = int(line)
            data_all[timestep] = {"data": [], "cols": []}
            continue

        elif line.startswith("ITEM: ATOMS"):
            reading_atoms = True
            columnas = line.split()[2:]
            data_all[timestep]["cols"] = columnas
            continue

        elif reading_atoms:
            vals = line.split()
            if len(vals) == len(columnas):
                data_all[timestep]["data"].append(vals)

    # Convertir cada timestep en un DataFrame
    for t in data_all:
        cols = data_all[t]["cols"]
        df = pd.DataFrame(data_all[t]["data"], columns=cols)
        for c in cols:
            df[c] = pd.to_numeric(df[c], errors='ignore')

        if columnas_interes is not None:
            df = df[[c for c in columnas_interes if c in df.columns]]

        data_all[t] = df

    # Si solo hay un timestep, devolver directamente el DataFrame
    if len(data_all) == 1:
        return list(data_all.values())[0]

    return data_all

#******************************************************************

#******************************************************************
# Columnas que queremos conservar

columnas = ['type', 'x', 'y', 'z', 'cluster_id']

camino='/home/matias/Escritorio/CONICET/bzt/Clusters_BZT/0.40-S100/Cluster_sample_cutoff_10/'
IN_FILES=['grouped_ss_' + str(i) + '.dump' for i in range(50)]
Temp = [10*i for i in range(1,41)]+[10*i for i in range(41,60,2)]  # Temperaturas correspondientes a cada archivo
Out_image = '0.40-S100/Fig_Cluster_cutoff_10.png'
Out_file = '0.40-S100/Cluster_summary_cutoff_10.csv'

resultados = []
for i, f in enumerate(IN_FILES):
    fname = camino + f
    df = leer_dump_seleccion(fname, columnas)

    # seleccionar clusters válidos
    df_valid = df.loc[df['cluster_id'] >= 1]

    if df_valid.empty:
        max_cluster = 0
        max_cluster_size = 0
        n_in_clusters = 0
        std_size = 0
    else:
        # contar átomos por cluster
        conteo = df_valid['cluster_id'].value_counts()
        max_cluster = conteo.idxmax()      # cluster_id del más grande
        max_cluster_size = conteo.max()    # tamaño del cluster más grande
        n_in_clusters = df_valid.shape[0]  # total de átomos en clusters
        std_size = conteo.std()            # desviación del tamaño de clusters
        n_clusters = df_valid['cluster_id'].nunique()  # número de clusters

    temp = Temp[i]*1.48
    resultados.append((f, temp, max_cluster, max_cluster_size, n_in_clusters, std_size,
                       n_clusters))

# Crear DataFrame y guardar
resumen = pd.DataFrame(resultados,
                       columns=['file', 'Temperature', 'max_cluster', 'max_size', 'n_in_clusters',
                                'std_size','n_clusters'])

resumen.to_csv(Out_file, index=False)
print("Resumen guardado en " + Out_file)

# === Graficar ===
fig, ax1 = plt.subplots(figsize=(8,5))

# Tamaño del cluster más grande con barra de error
ax1.errorbar(resumen['Temperature'], resumen['max_size'], yerr=resumen['std_size'],
             fmt='o-', capsize=5, label='Tamaño cluster más grande ± std', color='tab:blue')

# Número total de átomos en clusters
ax1.plot(resumen['Temperature'], resumen['n_in_clusters'], 's--', label='Átomos en clusters', color='tab:orange')
ax1.vlines(x=380, ymax=resumen['n_in_clusters'].max(),ymin=0, color='gray', linestyle=':', label='Transición (380K)')
ax1.vlines(x=650, ymax=resumen['n_in_clusters'].max(),ymin=0, color='red', linestyle=':', label='Burns (650K)')
ax1.set_xlabel('Temperatura [K]')
ax1.set_ylabel('Número de átomos')
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.legend(loc='upper left')
ax1.minorticks_on()

# Segundo eje para el número de clusters
ax2 = ax1.twinx()
ax2.plot(resumen['Temperature'], resumen['n_clusters'], 'd-.', label='Número de clusters', color='tab:green')
ax2.set_ylabel('Número de clusters')
ax2.legend(loc='upper right')
ax2.minorticks_on()

plt.title('Evolución de clusters con la temperatura')
fig.tight_layout()
plt.savefig(Out_image, dpi=300)
plt.show()



