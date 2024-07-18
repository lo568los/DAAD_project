import numpy as np
import matplotlib.pyplot as plt




N_list = [10]
theta = 1.07  #pass the true values here
theta_k = 0.79
t = 25

data1 = np.loadtxt(f"scaled_codes/data/N = 10, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}, t= {t}_corr.txt")

vals = data1[:,1]

plt.plot(range(N_list[0]),vals)
plt.xlabel("Position")
plt.ylabel(f"Correlator value at t = {t}")
plt.savefig(f"scaled_codes/plots/Corr_plot_N = {N_list[0]}", dpi =500)

