## This is the plotting file which will generate all the plots

## Author: Soumyadeep Sarma

###################    Step 1: Import all the libraries    ###########################

import numpy as np
import matplotlib.pyplot as plt

###################    Step 2: Get data from the text files    ###########################
N = 6
theta= 1.07 #pass the true values here
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

conc_vals = [0]*max_trotter_steps
vne_vals = [0]*max_trotter_steps
ratio_vals = [0]*max_trotter_steps

ratio_vals2 = [0]*max_trotter_steps

conc_vals2 = [0]*max_trotter_steps
vne_vals2 = [0]*max_trotter_steps

conc_vals3 = [0]*max_trotter_steps
vne_vals3 = [0]*max_trotter_steps

for i in range(max_trotter_steps):
    data = np.loadtxt(f"../scaled_codes/data/N = {N}, theta = {theta}, theta_k = {theta_k}, t = {i}_sz_TS.txt")
    conc_vals[i] = np.sqrt(1-data[1]**2)
    vne_vals[i] = -0.5*(np.log((1-data[1]**2)/4) + np.abs(data[1])*np.log((1+np.abs(data[1]))/(1-np.abs(data[1]))))
    #ratio_vals[i] = vne_vals[i]/conc_vals[i]

for i in range(max_trotter_steps):
    data2 = np.loadtxt(f"../scaled_codes/data/N = 10, theta = {theta}, theta_k = {theta_k}, t = {i}_sz_TS.txt")
    conc_vals2[i] = np.sqrt(1-data2[1]**2)
    vne_vals2[i] = -0.5*(np.log((1-data2[1]**2)/4) + np.abs(data2[1])*np.log((1+np.abs(data2[1]))/(1-np.abs(data2[1]))))

for i in range(max_trotter_steps):
    data2 = np.loadtxt(f"../scaled_codes/data/N = 8, theta = {theta}, theta_k = {theta_k}, t = {i}_sz_TS.txt")
    conc_vals3[i] = np.sqrt(1-data2[1]**2)
    vne_vals3[i] = -0.5*(np.log((1-data2[1]**2)/4) + np.abs(data2[1])*np.log((1+np.abs(data2[1]))/(1-np.abs(data2[1]))))

"""for i in range(max_trotter_steps):
    data3 = np.loadtxt(f"scaled_codes/data/N = {N}, theta = {theta}, theta_k = {theta_k}, t = {i}_sz.txt")
    ratio_vals2[i] = (1+data3[1])/conc_vals[i]"""




###################    Step 3: Plot the data   ###########################

              

plt.plot(range(max_trotter_steps),conc_vals,"r-", label = "Concurrence, N = 6")
#plt.plot(range(max_trotter_steps),vne_vals,"r--", label = "Von Neumann, N = 6")
plt.plot(range(max_trotter_steps),conc_vals2,"b-", label = "Concurrence, N = 10")
#plt.plot(range(max_trotter_steps),vne_vals2,"b--", label = "Von Neumann, N = 10")
plt.plot(range(max_trotter_steps),conc_vals3,"--",color = "#ee9190", label = "Concurrence, N = 8")

#plt.plot(range(max_trotter_steps),ratio_vals2,"r.", label = "Ratio of 1+Sz and Conc., N = 6")

#plt.plot(range(max_trotter_steps),conc_vals2,"b-", label = "Concurrence, N = 6")
#plt.plot(range(max_trotter_steps),vne_vals2,"b--", label = "Von Neumann, N = 6")
"""plt.xlabel("Time(trotter steps)")
plt.ylabel(r"Entanglement between subsystems")
plt.legend()
plt.title(f"Entanglement measure v/s time for " +  r'FS State, $\theta =$' +  f"{round(theta,2)}, " + r' $\theta_k =$' + f"{round(theta_k,2)}")"""
plt.savefig(f"Ent_plot_FS_panel_sum,theta=pi3,theta_k=pi6_tol", dpi =3000)
plt.close()

                







