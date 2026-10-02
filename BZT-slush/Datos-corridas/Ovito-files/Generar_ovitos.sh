#!/bin/bash
set -e   # corta si algún comando falla

zr=Zr0.60

cp ../CONFIGs/CONFIG-$zr CONFIG

for i in $(seq 10 10 600); do
    cp ../P-local/$zr/P-local-$i P-local
    ./ov3.x
    mv ovito.xyz $zr/ovito-$i.xyz
done

rm CONFIG P-local

