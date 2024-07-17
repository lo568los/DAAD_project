#!/usr/bin/env bash 

#@ job_name = KondoH
#@ job_type = MPICH                        
#@ output = Outputs/$(jobid).output_file  
#@ error  = Errors/$(jobid).error_file      
#@ tasks_per_node = 100
#@ environment = COPY_ALL                   
#@ class = 128core
#@ notification = complete
#@ notify_user = e.koenig@fkf.mpg.de
#@ queue 

export OMP_NUM_THREADS=100
N=10
theta=1.07
theta_k=0.79
max_trotter_steps=100


python3 batched_h.py $N $theta $theta_k $max_trotter_steps
