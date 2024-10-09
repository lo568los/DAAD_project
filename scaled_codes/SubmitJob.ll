#!/usr/bin/env bash 

#@ job_name = TestCodes
#@ job_type = MPICH                        
#@ output = Outputs/$(jobid).output_file  
#@ error  = Errors/$(jobid).error_file      
#@ tasks_per_node = 8
#@ environment = COPY_ALL                   
#@ class = 32core
#@ notification = complete
#@ notify_user = ssoumyadeep@iisc.ac.in
#@ queue 

export OMP_NUM_THREADS=8
N=6
theta=1.07
theta_k=0.79


python3 ed.py $N $theta $theta_k 
