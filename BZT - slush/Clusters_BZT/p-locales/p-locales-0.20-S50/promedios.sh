#!/bin/bash


for ii in `seq 410 20 600`


do


cp ../$ii/P-local entrada.in

./loc

mv salida.out P-local-$ii 

rm entrada.in

done



