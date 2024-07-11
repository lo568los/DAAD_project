## This is the plotting file which will generate all the plots

## Author: Soumyadeep Sarma

###################    Step 1: Import all the libraries    ###########################

import numpy as np
import matplotlib.pyplot as plt

###################    Step 2: Get data from the text files    ###########################
N_list = [6]
theta_list = [1.07]  #pass the true values here
thetak_list = [0.79]
max_trotter_steps = 5

def plot_corr_space(pos,corr_super):   # For corr vs time
    vals = corr_super[pos-1]
    #print(vals)
    plt.plot(range(max_trotter_steps),vals)
    plt.xlabel("Time(trotter steps)")
    plt.ylabel(r"$\langle S_z(x,t)S_z(0,t) \rangle$")
    plt.title("Correlator expectation as a function of time")
    plt.savefig(f"scaled_codes/Correlator time, N = {N}")
    plt.close()

def plot_corr_time(t, corr_super): # For corr vs pos
    corr_super = np.array(corr_super)
    vals = corr_super[:,t]
    #print(vals)
    plt.plot(range(1,N+1),vals)
    plt.xlabel("Position of spin")
    plt.ylabel(r"$\langle S_z(x,{t})S_z(0,{t}) \rangle$" +f"at t = {t}")
    plt.title("Correlator expectation as a function of space")
    plt.savefig(f"scaled_codes/Correlator space, N = {N}")
    plt.close()


for N in N_list:
    for theta in theta_list:
        for theta_k in thetak_list:
            if theta_k > theta:
                continue
            else:


                data1 = np.loadtxt(f"scaled_codes/N = {N}, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}_sz.txt")
                data2 = np.loadtxt(f"scaled_codes/N = {N}, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}_h.txt")
                #data3 = np.loadtxt(f"scaled_codes/N = {N}, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}_corr.txt")
                data4 = np.loadtxt(f"scaled_codes/N = {N}, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}_ent.txt")
                sz_vals = data1[:,1]
                h_vals = data2[:,1]
                conc_vals = data4[:,1]
                vne_vals = data4[:,2]
                #conc_vals = [x.real for x in conc_vals]
                #vne_vals = [y.real for y in conc_vals]

                """corr_super = []
                for i in range(1,N+1):
                    corr_super.append(data3[:,i])"""




###################    Step 3: Plot the data   ###########################

                plt.plot(range(max_trotter_steps),sz_vals)
                plt.xlabel("Time(trotter steps)")
                plt.ylabel(r"$\langle S^z_{imp}(t)\rangle$")
                plt.title(f"Impurity magnetization v/s time for N = {N}, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}")
                plt.savefig(f"scaled_codes/Sz plot, N = {N}", dpi =500)
                plt.close()

                plt.plot(range(max_trotter_steps),h_vals)
                plt.xlabel("Time(trotter steps)")
                plt.ylabel(r"$\langle H \rangle (t)$")
                plt.title(f"Impurity magnetization v/s time for N = {N}, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}")
                plt.savefig(f"scaled_codes/H_exp plot, N = {N}", dpi =500)
                plt.close()

                plt.plot(range(max_trotter_steps),conc_vals, label = "Concurrence")
                plt.plot(range(max_trotter_steps),vne_vals, label = "Von Neumann")
                plt.xlabel("Time(trotter steps)")
                plt.ylabel(r"Entanglement between subsystems")
                plt.title(f"Entanglement measures between impurity and leads for N = {N}, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}")
                plt.legend()
                plt.savefig(f"scaled_codes/Entanglement plot, N = {N}", dpi = 500)
                plt.close()



            #plot_corr_space(1,corr_super)
            #plot_corr_time(2,corr_super)







