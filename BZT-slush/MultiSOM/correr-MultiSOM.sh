#!/bin/bash



for i in $(seq 0 19) 
do
    cp feats_ss50/s50.$i.dump feat_sample.dump
    cp ss_xx_files/0.40-S0/ss_0.dump ss_xx.dump
    python3 main.py 
    mv SOM_ss_xx.dump SOM_files/SOM_ss_0_trained_s50_$i.dump
    python3  clasificar_fases_SOM.py
done

