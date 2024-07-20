## This is the plotting file which will generate all the plots

## Author: Soumyadeep Sarma

###################    Step 1: Import all the libraries    ###########################

import numpy as np
import matplotlib.pyplot as plt

###################    Step 2: Get data from the text files    ###########################
N = 6
theta= 1.07  #pass the true values here
theta_k = 0.79
max_trotter_steps = 100

def plot_corr_space(pos,corr_super):   # For corr vs time
    vals = corr_super[pos-1]
    #print(vals)
    plt.plot(range(20),vals, label = f"Position (x) = {pos}")
    

def plot_corr_time(t, corr_super): # For corr vs pos
    corr_super = np.array(corr_super)
    vals = corr_super[:,t]
    #print(vals)
    plt.plot(range(1,N+1),vals)
    plt.xlabel("Position of spin")
    plt.ylabel(r"$\langle S_z(x,{t})S_z(0,{t}) \rangle$" +f"at t = {t}")
    plt.title("Correlator expectation as a function of space")
    plt.savefig(f"scaled_codes/plots/Correlator space, N = {N}")
    plt.close()

#i = 0

sz_vals = [0]*max_trotter_steps

for i in range(max_trotter_steps):
    data = np.loadtxt(f"scaled_codes/data/N = {N}, theta = {theta}, theta_k = {theta_k}, t = {i}_sz.txt")
    sz_vals[i] = data[:,1][0]
        




###################    Step 3: Plot the data   ###########################

              

plt.plot(range(max_trotter_steps),sz_vals)
plt.xlabel("Time(trotter steps)")
plt.ylabel(r"$\langle S^z_{imp} \rangle (t)$")
plt.title(f"Impurity Magnetization v/s time for " +  r'FS State, $\theta =$' +  f"{round(theta,2)}, " + r' $\theta_k =$' + f"{round(theta_k,2)}")
plt.savefig(f"scaled_codes/plots/Sz plot_FS , N = {N}", dpi =500)
plt.close()

                







