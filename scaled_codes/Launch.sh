#!/usr/bin/env bash

for Nloop in 6
do
    for theta_loop in 1.07 
    do
        for thetak_loop in 0 0.12 0.22 0.32 0.42 0.45 0.52 0.55 0.62 0.65 0.72 0.82
        do
            

            sed -i "s/N=.*/N=${Nloop}/g" SubmitJob.ll
            sed -i "s/theta=.*/theta=${theta_loop}/g" SubmitJob.ll
            sed -i "s/theta_k=.*/theta_k=${thetak_loop}/g" SubmitJob.ll
            llsubmit SubmitJob.ll
            
        done
    done
done
