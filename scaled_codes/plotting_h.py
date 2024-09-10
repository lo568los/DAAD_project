## This is the plotting file which will generate all the plots

## Author: Soumyadeep Sarma

###################    Step 1: Import all the libraries    ###########################

import numpy as np
import matplotlib.pyplot as plt

###################    Step 2: Get data from the text files    ###########################
N = 10
theta= 0.79 #pass the true values here
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

h_vals = [0]*max_trotter_steps

for i in range(max_trotter_steps):
    data = np.loadtxt(f"scaled_codes/data/N = {N}, theta = {theta}, theta_k = {theta_k}, t = {i}_h.txt")
    #print(data)
    h_vals[i] = data[1]/N

#data1 = np.loadtxt(f"scaled_codes/data/N = {N}, theta = {theta}, theta_k = {theta_k}_h.txt")
data2 = np.loadtxt(f"scaled_codes/data/N = 6, theta = {theta}, theta_k = {theta_k}_h.txt")

#h_vals1 = data1[:,1]/N
h_vals2 = data2[:,1]/6


        




###################    Step 3: Plot the data   ###########################

              

plt.plot(range(max_trotter_steps),h_vals,"b-", label = "N=10, FS Isotropic")
#plt.plot(range(max_trotter_steps),h_vals1, label = "N=10, old")
plt.plot(range(max_trotter_steps),h_vals2,"r-", label = "N=6, FS Isotropic")
plt.xlabel("Time(trotter steps)")
plt.ylabel(r"$\langle H \rangle (t)$/N")
plt.legend()
plt.grid()
plt.title(f"Hamiltonian expectation v/s time for " +  r'$\theta = \pi/4$ and'  + r' $\theta_k = \pi/4$' )
plt.savefig(f"scaled_codes/plots/H_exp plot_FS_fin1", dpi =500)
plt.close()

                







