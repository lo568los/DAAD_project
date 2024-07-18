

import numpy as np
import matplotlib.pyplot as plt




N_list = [10]
theta_list = [1.07]  #pass the true values here
thetak_list = [0.79]
max_trotter_steps = 25

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


for N in N_list:
    for theta in theta_list:
        for theta_k in thetak_list:
            if theta_k > theta:
                continue
            else:
                data3 = np.loadtxt(f"scaled_codes/data/N = {N}, theta = 0.79, theta_k = 0.53, t = 25_corr.txt")
                corr_super = []
                for i in range(1,N+1):
                    corr_super.append(data3[:,i])






                pos_list = [1,2,3]
                for pos in pos_list:
                    plot_corr_space(pos,corr_super)
                    plt.xlabel("Time(trotter steps)")
                    plt.ylabel(r"$\langle S_z(x,t)S_z(0,t) \rangle$")
                    plt.title("Correlator expectation as a function of time")
                    plt.legend()
                    plt.savefig(f"scaled_codes/plots/Correlator time, N = {N}")
                    plt.close()
                plot_corr_time(2,corr_super)