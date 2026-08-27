#!/usr/bin/env bash


for ii in `seq 10 10 350`


do


  cp /home/rodrigo/MEGA/Ba1-xCaxTiO3-2025/dl-poly/solucion-solida-con-potenciales-nuevos/BTO-Ca-2025-1/corrida-5/Ba0.70Ca0.30TiO3-001/corridas/temperatura/$ii/P-total borrador-1
  tail -n 100 /home/rodrigo/MEGA/Ba1-xCaxTiO3-2025/dl-poly/solucion-solida-con-potenciales-nuevos/BTO-Ca-2025-1/corrida-5/Ba0.70Ca0.30TiO3-001/corridas/temperatura/mas-pasos-para-permitividad/$ii/P-total > borrador-2

  cat borrador-1 borrador-2 > P-total
  

  cp ../promedio-parametros-de-red/entrada-vectores-de-red-desde-statis/entrada-vectores-de-red-${ii}.dat cell-vectors

python3  permittivity.py P-total -T $ii --lat cell-vectors --block-len 20 --out epsilon_vs_T.csv

tail -n 1 epsilon_vs_T.csv >> e-vs-T.csv 

done
