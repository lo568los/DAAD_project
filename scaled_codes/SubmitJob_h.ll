#!/usr/bin/env bash 

#@ job_name = KondoH
#@ job_type = MPICH                        
#@ output = Outputs/$(jobid).output_file  
#@ error  = Errors/$(jobid).error_file      
#@ tasks_per_node = 1
#@ environment = COPY_ALL                   
#@ class = 32core
#@ notification = complete
#@ notify_user = e.koenig@fkf.mpg.de
#@ queue 

export OMP_NUM_THREADS=1
N=10
theta=0.79
theta_k=0.52
max_trotter_steps=400


python3 test_h.py $N $theta $theta_k $max_trotter_steps
