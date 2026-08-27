import numpy as np

NUM_CELDAS = 21952
pasos = np.arange(310000, 1000001, 10000)

with open("P-local.txt", "r") as f:
    data = f.read().splitlines()

P_pasos = []
for step in range(len(pasos)):
    inicio = step * (NUM_CELDAS + 1) + 1  # +1 para saltear la cabecera
    bloque = np.array([data[a].split() for a in range(inicio, inicio + NUM_CELDAS)], dtype=float)
    P_pasos.append(bloque)

# promedio temporal: (70, 21952, 3) → (21952, 3)
P_promedio = np.mean(np.array(P_pasos), axis=0)

# módulo de cada celda
P_modulo = np.linalg.norm(P_promedio, axis=1)  # shape (21952,)

with open("P-local-promedio.txt", "w") as f:
    for i in range(NUM_CELDAS):
        f.write(f"{P_promedio[i, 0]:12.5f} {P_promedio[i, 1]:12.5f} {P_promedio[i, 2]:12.5f} {P_modulo[i]:12.5f}\n")
