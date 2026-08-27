#!/bin/bash
for i in $(seq 310 10 400)
do
  cp p-locales/p-locales-0.40-S0/P-local-$i P-local
  gfortran -o ovito3.f90 ov3.x

  ./ov3.x
  mv ovito.xyz ovito-$i.xyz

done

for i in $(seq 410 20 590)
do
  cp p-locales/p-locales-0.40-S0/P-local-$i P-local
  gfortran -o ovito3.f90 ov3.x

  ./ov3.x
  mv ovito.xyz ovito-$i.xyz

done