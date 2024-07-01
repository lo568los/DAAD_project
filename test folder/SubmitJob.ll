#!/usr/bin/env bash 

#@ job_name = Doesitwork
#@ job_type = MPICH                        
#@ output = Outputs/$(jobid).output_file  
#@ error  = Errors/$(jobid).error_file      
#@ tasks_per_node = 8
#@ environment = COPY_ALL                   
#@ class = 32core
#@ notification = never                  
#@ queue 

export OMP_NUM_THREADS=8
Ahub=-4.7
Dopinghub=1
t_primehub=1
Thub=1

python3 test.py $Ahub
