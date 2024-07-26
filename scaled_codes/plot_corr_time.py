import numpy as np
import matplotlib.pyplot as plt




N_list = [10]
theta = 1.07  #pass the true values here
theta_k = 0.79
t = 50

data1 = np.loadtxt(f"scaled_codes/data/N = 10, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}, t= {t}_corr.txt")
data2 = np.loadtxt(f"scaled_codes/data/N = 10, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}, t= 25_corr.txt")
data3 = np.loadtxt(f"scaled_codes/data/N = 10, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}, t= 10_corr.txt")
data4 = np.loadtxt(f"scaled_codes/data/N = 10, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}, t= 30_corr.txt")
data5 = np.loadtxt(f"scaled_codes/data/N = 10, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}, t= 40_corr.txt")

vals = data1[:,1]
vals2 = data2[:,1]
vals3 = data3[:,1]
vals4 = data4[:,1]
vals5 = data5[:,1]

plt.plot(range(1,N_list[0]+1),vals, label = "t = 50")
plt.plot(range(1,N_list[0]+1),vals2, label = "t = 25")
plt.plot(range(1,N_list[0]+1),vals3, label = "t = 10")
plt.plot(range(1,N_list[0]+1),vals4, label = "t = 30")
plt.plot(range(1,N_list[0]+1),vals5, label = "t = 40")
plt.xlabel("Position")
plt.ylabel(f"Correlator value")
plt.title("Correlator expectation at fixed time")
plt.legend()
plt.savefig(f"scaled_codes/plots/Corr_plot34_N = {N_list[0]}", dpi =500)

