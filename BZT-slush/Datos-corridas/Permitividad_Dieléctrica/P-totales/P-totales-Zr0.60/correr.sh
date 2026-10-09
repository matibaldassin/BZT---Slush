#!/bin/bash

for i in $(seq 20 20 280); do
	cp ../../Cooling-${i}/P-total 0kV_mm/P-total-${i}
	cp ../${i}/5kV_mm/P-total 5kV_mm/P-total-${i}
	cp ../${i}/10kV_mm/P-total 10kV_mm/P-total-${i}
done

