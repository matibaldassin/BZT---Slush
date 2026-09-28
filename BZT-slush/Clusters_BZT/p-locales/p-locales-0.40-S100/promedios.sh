#!/bin/bash


for ii in `seq 410 20 600`


do


cp ../heating-2-$ii/P-local entrada.in

./loc

mv salida.out P-local-heating-2-$ii 

rm entrada.in

done



