#!/bin/bash
# inicio_25_50_25
for ii in `seq 25 25 50`
do
cp ../inicio-${ii}_kV_por_mm/P-total entrada.dat
./prom-P
(echo -n $ii "   " ; tail -n 1 pola-total-prom) >> polarizaciones
rm entrada.dat pola-total-prom
done

# inicio_100_250_50
for ii in `seq 100 50 250`
do
cp ../inicio-${ii}_kV_por_mm/P-total entrada.dat
./prom-P
(echo -n $ii "   " ; tail -n 1 pola-total-prom) >> polarizaciones
rm entrada.dat pola-total-prom
done

# menos_200_100_-50
for ii in `seq 200 -50 100`
do
cp ../menos-${ii}_kV_por_mm/P-total entrada.dat
./prom-P
(echo -n $ii "   " ; tail -n 1 pola-total-prom) >> polarizaciones
rm entrada.dat pola-total-prom
done

# menos_90_-100_-10
for ii in `seq 90 -10 -100`
do
cp ../menos-${ii}_kV_por_mm/P-total entrada.dat
./prom-P
(echo -n $ii "   " ; tail -n 1 pola-total-prom) >> polarizaciones
rm entrada.dat pola-total-prom
done

# menos -150_-250_-50
for ii in `seq -150 -50 -250`
do
cp ../menos-${ii}_kV_por_mm/P-total entrada.dat
./prom-P
(echo -n $ii "   " ; tail -n 1 pola-total-prom) >> polarizaciones
rm entrada.dat pola-total-prom
done

# mas_-200_-100_50
for ii in `seq -200 50 -100`
do
cp ../mas-${ii}_kV_por_mm/P-total entrada.dat
./prom-P
(echo -n $ii "   " ; tail -n 1 pola-total-prom) >> polarizaciones
rm entrada.dat pola-total-prom
done

# mas_-90_10_10
for ii in `seq -90 10 10`
do
cp ../mas-${ii}_kV_por_mm/P-total entrada.dat
./prom-P
(echo -n $ii "   " ; tail -n 1 pola-total-prom) >> polarizaciones
rm entrada.dat pola-total-prom
done

# mas_20_60_10
for ii in `seq 20 10 60`
do
cp ../mas-${ii}_kV_por_mm/P-total entrada.dat
./prom-P-extra
(echo -n $ii "   " ; tail -n 1 pola-total-prom) >> polarizaciones
rm entrada.dat pola-total-prom
done

# mas_70_100_10
for ii in `seq 70 10 100`
do
cp ../mas-${ii}_kV_por_mm/P-total entrada.dat
./prom-P
(echo -n $ii "   " ; tail -n 1 pola-total-prom) >> polarizaciones
rm entrada.dat pola-total-prom
done

# mas_150_250_50
for ii in `seq 150 50 250`
do
cp ../mas-${ii}_kV_por_mm/P-total entrada.dat
./prom-P
(echo -n $ii "   " ; tail -n 1 pola-total-prom) >> polarizaciones
rm entrada.dat pola-total-prom
done
