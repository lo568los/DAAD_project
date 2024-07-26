## This is the plotting file which will generate all the plots

## Author: Soumyadeep Sarma

###################    Step 1: Import all the libraries    ###########################

import numpy as np
import matplotlib.pyplot as plt

###################    Step 2: Get data from the text files    ###########################
N = 6
theta= 0.79 #pass the true values here
theta_k = 0.52
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

conc_vals = [0]*max_trotter_steps
vne_vals = [0]*max_trotter_steps
ratio_vals = [0]*max_trotter_steps

ratio_vals2 = [0]*max_trotter_steps

conc_vals2 = [0]*max_trotter_steps
vne_vals2 = [0]*max_trotter_steps

for i in range(max_trotter_steps):
    data = np.loadtxt(f"scaled_codes/data/N = {N}, theta = {theta}, theta_k = {theta_k}, t = {i}_ent.txt")
    conc_vals[i] = data[1]
    vne_vals[i] = data[2]
    ratio_vals[i] = vne_vals[i]/conc_vals[i]

for i in range(max_trotter_steps):
    data2 = np.loadtxt(f"scaled_codes/data/N = 6, theta = {theta}, theta_k = {theta_k}, t = {i}_ent.txt")
    conc_vals2[i] = data2[1]
    vne_vals2[i] = data2[2]

for i in range(max_trotter_steps):
    data3 = np.loadtxt(f"scaled_codes/data/N = {N}, theta = {theta}, theta_k = {theta_k}, t = {i}_sz.txt")
    ratio_vals2[i] = (1+data3[1])/conc_vals[i]




###################    Step 3: Plot the data   ###########################

              

plt.plot(range(max_trotter_steps),conc_vals,"r-", label = "Concurrence, N = 6")
plt.plot(range(max_trotter_steps),vne_vals,"b--", label = "Von Neumann, N = 6")
plt.plot(range(max_trotter_steps),ratio_vals2,"r.", label = "Ratio of 1+Sz and Conc., N = 6")

#plt.plot(range(max_trotter_steps),conc_vals2,"b-", label = "Concurrence, N = 6")
#plt.plot(range(max_trotter_steps),vne_vals2,"b--", label = "Von Neumann, N = 6")
plt.xlabel("Time(trotter steps)")
plt.ylabel(r"Entanglement between subsystems")
plt.legend()
plt.title(f"Entanglement measure v/s time for " +  r'FS State, $\theta =$' +  f"{round(theta,2)}, " + r' $\theta_k =$' + f"{round(theta_k,2)}")
plt.savefig(f"scaled_codes/plots/Ent plot_FS_ratio3", dpi =500)
plt.close()

                







