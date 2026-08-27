#!/bin/bash


for ii in `seq 10 10 350`


do


cp ../$ii/P-local entrada.in

./loc

mv salida.out P-local-$ii 

rm entrada.in

done



