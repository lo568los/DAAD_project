## This is the plotting file which will generate all the plots

## Author: Soumyadeep Sarma

###################    Step 1: Import all the libraries    ###########################

import numpy as np
import matplotlib.pyplot as plt

###################    Step 2: Get data from the text files    ###########################
N = 10
theta= 1.07  #pass the true values here
theta_k = 0.52
max_trotter_steps = 1000

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
sz_vals2 = [0]*max_trotter_steps
sz_vals3 = [0]*max_trotter_steps
sz_vals4 = [0]*max_trotter_steps
sz_vals5 = [0]*max_trotter_steps

for i in range(max_trotter_steps):
    #data = np.loadtxt(f"scaled_codes/data/N = {N}, theta = {theta}, theta_k = {theta_k}, t = {i}_sz_TS.txt")
    #sz_vals[i] = data[1]
    data = np.loadtxt(f"scaled_codes/data/N = {N}, theta = {theta}, theta_k = {theta_k}, t = {i}_sz.txt")
    sz_vals4[i] = data[1]



"""for i in range(max_trotter_steps):
    data3 = np.loadtxt(f"scaled_codes/data/N = 8, theta = {theta}, theta_k = {theta_k}, t = {i}_sz_TS.txt")
    sz_vals3[i] = data3[1]

for i in range(max_trotter_steps):
    data2 = np.loadtxt(f"scaled_codes/data/N = 6, theta = {theta}, theta_k = {theta_k}, t = {i}_sz.txt")
    sz_vals2[i] = data2[1]  
    #data2 = np.loadtxt(f"scaled_codes/data/N = 6, theta = {theta}, theta_k = {theta_k}, t = {i}_sz.txt")
    #sz_vals5[i] = data2[1]  """  




###################    Step 3: Plot the data   ###########################

              


plt.plot(range(max_trotter_steps),sz_vals4, "b-", label = "N=6,FS")
"""plt.plot(range(max_trotter_steps),sz_vals4, "b-", label = "N=6,FS")
plt.plot(range(max_trotter_steps),sz_vals3,"--",color='orange', label = "N=8,TS")
plt.plot(range(max_trotter_steps),sz_vals,"g--", label = "N=10,TS")
plt.plot(range(max_trotter_steps),sz_vals5,"g-", label = "N=10,FS")"""
plt.xlabel("Time(trotter steps)")
plt.ylabel(r"$\langle S^z_{imp} \rangle (t)$")
plt.legend()
plt.title(f"Impurity Magnetization v/s time for " +  r'$\theta =$' +  f"{round(theta,2)}, " + r' $\theta_k =$' + f"{round(theta_k,2)}")
plt.savefig(f"scaled_codes/plots/Sz plot_FS_rec5 ", dpi =500)
plt.close()

                







