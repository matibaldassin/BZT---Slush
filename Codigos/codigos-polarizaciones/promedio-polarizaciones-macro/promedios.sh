#!/bin/bash


for ii in `seq 280 10 350`


do


cp ../$ii/P-total entrada.dat

./prom-P

(echo -n $ii "   " ; tail -n 1 pola-total-prom) >> polarizaciones

rm entrada.dat pola-total-prom

done



