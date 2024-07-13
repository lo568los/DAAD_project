#!/usr/bin/env bash 

#@ job_name = TestCodes
#@ job_type = MPICH                        
#@ output = Outputs/$(jobid).output_file  
#@ error  = Errors/$(jobid).error_file      
#@ tasks_per_node = 8
#@ environment = COPY_ALL                   
#@ class = 128core
#@ notification = complete
#@ notify_user = soumyadeepsarma3@gmail.com
#@ queue 

export OMP_NUM_THREADS=8
N=14
theta=1.07
theta_k=0.79
max_trotter_steps=50


python3 test_code_multi.py $N $theta $theta_k $max_trotter_steps
