#!/usr/bin/env bash

for Aloop in 1 2 3 4 5 6 7 8 9 10 -4.7
do
    for t_primeloop in 1
    do
        for Tloop in 1
        do
            for Dopingloop in 1
            do

                sed -i "s/Ahub=.*/Ahub=${Aloop}/g" SubmitJob.ll
                sed -i "s/Thub=.*/Thub=${Tloop}/g" SubmitJob.ll
                sed -i "s/t_primehub=.*/t_primehub=${t_primeloop}/g" SubmitJob.ll
                sed -i "s/Dopinghub=.*/Dopinghub=${Dopingloop}/g" SubmitJob.ll
                llsubmit SubmitJob.ll
            done
        done
    done
done
