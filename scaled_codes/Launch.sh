#!/usr/bin/env bash

for Nloop in 6
do
    for theta_loop in 1.07 0.79
    do
        for thetak_loop in 0.79 0.52
        do
            

            sed -i "s/N=.*/N=${Nloop}/g" SubmitJob.ll
            sed -i "s/theta=.*/theta=${theta_loop}/g" SubmitJob.ll
            sed -i "s/theta_k=.*/theta_k=${thetak_loop}/g" SubmitJob.ll
            llsubmit SubmitJob.ll
            
        done
    done
done
