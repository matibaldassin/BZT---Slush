NUM_CELDAS = 21952
pasos = iter(range(310000, 1000001, 10000))

with open("P-local.txt", "w") as out:
    for i in range(4, 11):
        with open(f"P-local-{i}", "r") as f:
            data = f.read().splitlines()
            for a in range(len(data)):
                if a % (NUM_CELDAS + 1) == 0:
                    out.write(f"{next(pasos)}\n")
                    continue
                values = data[a].split()
                out.write(f"{float(values[1]):12.5f} {float(values[2]):12.5f} {float(values[3]):12.5f}\n")
