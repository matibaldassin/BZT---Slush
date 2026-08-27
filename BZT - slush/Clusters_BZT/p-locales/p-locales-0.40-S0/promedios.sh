#!/bin/bash


for ii in `seq 10 10 10`


do


cp ../heating-2-$ii/P-local entrada.in

./loc

mv salida.out P-local-$ii 

rm entrada.in

done



