#!/bin/bash


for ii in `seq 10 10 310`


do


cp ../P-local-$ii entrada.dat

./loc


if [ $ii -eq 10 ]
then
	head -n 1 salida.dat > polarizaciones
        (echo -n $ii "   " ; tail -n 1 salida.dat) >> polarizaciones
fi



if [ $ii -gt 10 ]
then

      (echo -n $ii "   " ; tail -n 1 salida.dat) >> polarizaciones

fi      

rm entrada.dat salida.dat

done



