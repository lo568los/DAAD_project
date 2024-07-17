#!/usr/bin/env bash 

#@ job_name = KondoEnt
#@ job_type = MPICH                        
#@ output = Outputs/$(jobid).output_file  
#@ error  = Errors/$(jobid).error_file      
#@ tasks_per_node = 50
#@ environment = COPY_ALL                   
#@ class = 128core
#@ notification = complete
#@ notify_user = e.koenig@fkf.mpg.de
#@ queue 

export OMP_NUM_THREADS=50
N=14
theta=1.07
theta_k=0.79
max_trotter_steps=50


python3 batched_ent.py $N $theta $theta_k $max_trotter_steps