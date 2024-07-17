#!/usr/bin/env bash

for Nloop in 10
do
    for theta_loop in 1.07
    do
        for thetak_loop in 0.79
        do
            for steps_loop in 100
            do

                sed -i "s/N=.*/N=${Nloop}/g" SubmitJob_ent.ll
                sed -i "s/theta=.*/theta=${theta_loop}/g" SubmitJob_ent.ll
                sed -i "s/theta_k=.*/theta_k=${thetak_loop}/g" SubmitJob_ent.ll
                sed -i "s/max_trotter_steps=.*/max_trotter_steps=${steps_loop}/g" SubmitJob_ent.ll
                llsubmit SubmitJob_ent.ll
            done
        done
    done
done
