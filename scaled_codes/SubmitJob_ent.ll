#!/usr/bin/env bash 

#@ job_name = KondoEnt
#@ job_type = MPICH                        
#@ output = Outputs/$(jobid).output_file  
#@ error  = Errors/$(jobid).error_file      
#@ tasks_per_node = 1
#@ environment = COPY_ALL                   
#@ class = 128core_new
#@ notification = complete
#@ notify_user = e.koenig@fkf.mpg.de
#@ queue 

export OMP_NUM_THREADS=1
N=10
theta=1.07
theta_k=0.79
max_trotter_steps=100


python3 test_ent.py $N $theta $theta_k $max_trotter_steps
