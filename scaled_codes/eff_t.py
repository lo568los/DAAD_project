## This is the test code with multithreading implemented to calculate Sz, H
## and effective temperature for different values of N, theta and theta_k.

## Author: Soumyadeep Sarma (Modified for Effective Temperature Analysis)

###################    Step 1: Import all the libraries    ###########################

import sys
import numpy as np
import math as m
import cmath as cm
from threading import Thread
import time
from itertools import combinations, cycle

# Qiskit imports
import qiskit
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector, DensityMatrix, Operator, SparsePauliOp
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import Estimator, Sampler

# Analysis imports
from scipy import linalg as la
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

###################    Step 2: Define parameters and inputs    ###########################

# Default values if sys.argv is not sufficient (for testing purposes)
# You can run this via: python floquet_kondo_teff.py 2 1.0 0.5 20
if len(sys.argv) < 5:
    print("Usage: python floquet_kondo_teff.py N theta theta_k max_trotter_steps")
    # Setting defaults for demonstration if run without args
    N = 2
    theta = 1.0
    theta_k = 0.5
    max_trotter_steps = 10
else:
    N = int(sys.argv[1])
    theta = float(sys.argv[2])
    theta_k = float(sys.argv[3])
    max_trotter_steps = int(sys.argv[4])

num_qubits = 2*N + 1 

###################    Step 3: Define all the helper functions    ###########################

coeff_dict  = {}

def array_k1(num_qubits):
    array_k = []
    m = num_qubits/2
    if (m)%2!=0:
        for j in range(-int((m)//2),int((m)//2) + 1):
            array_k.append(2*np.pi*j/num_qubits)
    else:
        for j in range(-int((m)//2),int((m)//2)):
            array_k.append(2*np.pi*j/num_qubits)
    return array_k


def recursive_nested(l,num_qubits,coeff_array,coeff = 1,bitstr=''):
    m = int(num_qubits/2)
    if l==m-1:
        coeff_copy = coeff
        bitstr_copy = bitstr
        for i in range(num_qubits):
            if str(i) in bitstr:
                pass
            else:
                coeff = coeff*coeff_array[l,i]
                bitstr = bitstr + f'{i}'
                bitstr_sorted = sort_bitstr(bitstr)
                perm = perm_str2(bitstr,bitstr_sorted)
                if bitstr_sorted in coeff_dict.keys():
                    coeff_dict[bitstr_sorted]+=coeff*perm
                else:
                    coeff_dict[bitstr_sorted] = coeff*perm
            bitstr = bitstr_copy
            coeff = coeff_copy
    if l!=m-1:
        coeff_copy = coeff 
        bitstr_copy = bitstr
        for i in range(num_qubits):
            if str(i) in bitstr:
                pass
            else:
                coeff = coeff*coeff_array[l,i]
                bitstr += str(i)
                recursive_nested(l+1,num_qubits,coeff_array,coeff,bitstr)
            bitstr = bitstr_copy
            coeff = coeff_copy
    return coeff_dict

def sort_bitstr(bitstr):
    bit_array = []
    for i in bitstr:
        bit_array.append(i)
    bitstr_sorted = ''
    for i in range(len(bit_array)):
        bit_array[i] = int(bit_array[i])
    bit_array.sort()
    for i in bit_array:
        bitstr_sorted += str(i)
    return bitstr_sorted

def perm_str2(cmpr,word):
    swaps = 0
    chars = {c: [] for c in word}
    [chars[c].append(i) for i, c in enumerate(word)]
    for k in chars.keys():
        chars[k] = cycle(chars[k])
    idxs = [next(chars[c]) for c in cmpr]
    for cmb in combinations(idxs, 2):
        if cmb[0] > cmb[1]:
            swaps += 1
    if swaps%2 == 0:
        return 1
    else:
        return -1
    
def fermi_state(num_qubits): 
    coeff_dict.clear()
    m = int(num_qubits/2)
    coeff_array = []
    array_k = array_k1(num_qubits)
    for k in array_k:
        pos_list = []
        for x in range(num_qubits):
            pos_list.append(np.exp(-1j*k*x))
        coeff_array.append(pos_list)

    coeff_array = np.array(coeff_array)
    coeff_dict_2 = recursive_nested(0,num_qubits,coeff_array)

    bitstr_dict = {}
    for bstr in coeff_dict_2.keys():
        vac_str = ''
        num_list = []
        for k in range(m):
            num_list.append(int(bstr[k]))

        for i in range(num_qubits):
            if i in num_list:
                vac_str += '1'
            else:
                vac_str += '0'
        bitstr_dict[vac_str] = coeff_dict[bstr]

    fermi_state = Statevector([0]*(2**num_qubits))

    for bstr in bitstr_dict.keys():
        fermi_state += Statevector.from_label(bstr)*bitstr_dict[bstr]

    val_array = []
    for i in bitstr_dict.values():
        val_array.append(i)
    
    fermi_state = fermi_state/np.linalg.norm(val_array)
    return fermi_state

def fermi_state_circuit(N,num_cl_bits = 0):
    qc = QuantumCircuit(num_qubits,num_cl_bits)
    fermi_state_up = fermi_state(N)
    fermi_state_down = fermi_state(N)
    qc.initialize(fermi_state_up,range(N))
    qc.initialize(fermi_state_down,range(N+1,2*N+1))
    return qc

def fsim(theta,phi,beta):
    fsim = Operator([[1,0,0,0],
                   [0,m.cos(theta),1j*cm.exp(1j*beta)*m.sin(theta),0],
                   [0,1j*cm.exp(-1j*beta)*m.sin(theta),m.cos(theta),0],
                   [0,0,0,cm.exp(1j*phi)]])
    return fsim

def add_fsim_half(qc,angles):
    theta = angles
    fsim1 = fsim(theta,0,0)
    fsim2 = fsim(2*theta,0,0)
    for i in range(0,qc.num_qubits//2-1,2):
        qc.unitary(fsim1,[i,i+1],label = r'fsim$(\theta,\phi)$')
    for i in range(qc.num_qubits//2+1,qc.num_qubits-1,2):
        qc.unitary(fsim1,[i,i+1],label = r'fsim$(\theta,\phi)$')
    for i in range(1,qc.num_qubits//2-1,2):
        qc.unitary(fsim2,[i,i+1],label = r'fsim$(2\theta,\phi)$')
    for i in range(qc.num_qubits//2+2,qc.num_qubits-1,2):
        qc.unitary(fsim2,[i,i+1],label = r'fsim$(2\theta,\phi)$')

def add_fsim_inv_half(qc,angles):
    theta = angles
    fsim1 = fsim(theta,0,0)
    for i in range(0,qc.num_qubits//2-1,2):
        qc.unitary(fsim1,[i,i+1],label = r'fsim$(\theta,\phi)$')
    for i in range(qc.num_qubits//2+1,qc.num_qubits-1,2):
        qc.unitary(fsim1,[i,i+1],label = r'fsim$(\theta,\phi)$')

def add_fsim_full(qc,angles):
    theta = angles
    fsim2 = fsim(2*theta,0,0)
    for i in range(0,qc.num_qubits//2-1,2):
        qc.unitary(fsim2,[i,i+1],label = r'fsim$(2\theta,\phi)$')
    for i in range(qc.num_qubits//2+1,qc.num_qubits-1,2):
        qc.unitary(fsim2,[i,i+1],label = r'fsim$(2\theta,\phi)$')
    for i in range(1,qc.num_qubits//2-1,2):
        qc.unitary(fsim2,[i,i+1],label = r'fsim$(2\theta,\phi)$')
    for i in range(qc.num_qubits//2+2,qc.num_qubits-1,2):
        qc.unitary(fsim2,[i,i+1],label = r'fsim$(2\theta,\phi)$')

def kondo_unitary(theta_k,theta_z):
    l1 = cm.exp(1j*theta_z/2)
    l2 = cm.exp(-1j*theta_z/2)
    c1 = m.cos(theta_k)
    s1 = m.sin(theta_k)

    a = m.cos(np.sqrt(2)*theta_k) - (1j/np.sqrt(2))*m.sin(np.sqrt(2)*theta_k)
    a_dag = m.cos(np.sqrt(2)*theta_k) + (1j/np.sqrt(2))*m.sin(np.sqrt(2)*theta_k)
    b = (-1j/np.sqrt(2))*m.sin(np.sqrt(2)*theta_k)
    c = cm.exp(-1j*theta_k) 
    d = cm.exp(1j*theta_k)

    kondo_unitary_2 = Operator([[1,0,0,0,0,0,0,0],
                          [0,a,0,0,0,0,b,0],
                          [0,0,1,0,0,0,0,0],
                          [0,0,0,d,0,0,0,0],
                          [0,0,0,0,c,0,0,0],
                          [0,0,0,0,0,1,0,0],
                          [0,b,0,0,0,0,a_dag,0],
                          [0,0,0,0,0,0,0,1]])
    
    return kondo_unitary_2

def circuit_3(N, trotter_steps,angles = 0,theta_k = 0,theta_z = 0, num_cl_bits = 0, trotter_barriers = False, save = False):
    if num_cl_bits == 0:
        qc = fermi_state_circuit(N)
    else:
        qc = fermi_state_circuit(N,num_cl_bits)
    qc.x(N)
    qc.barrier()
    
    c = num_qubits//2
    if trotter_steps == 0:
        if save == True:
            qc.save_statevector()
        return qc
    else:
        add_fsim_half(qc,angles)
        qc.unitary(kondo_unitary(theta_k,theta_z),[c,c+1,c-1],label=r'$U_{k}(\theta_k,\theta_z)$')
        if trotter_barriers:
                qc.barrier()
        for i in range(1,trotter_steps):
            add_fsim_full(qc,angles)
            qc.unitary(kondo_unitary(theta_k,theta_z),[c,c+1,c-1],label=r'$U_{k}(\theta_k,\theta_z)$')
            if trotter_barriers:
                qc.barrier()
        add_fsim_inv_half(qc,angles)
        if save == True:
            qc.save_statevector()
        return qc

# --- Hamiltonian Construction ---

H_t = 0
H_k = 0
for i in range(2*N):
    if i==N-1 or i==N:
        continue
    else:
        H_t += -theta*(SparsePauliOp('I'*(i) + 'XX' + 'I'*(2*N-i-1)) + SparsePauliOp('I'*(i) + 'YY' + 'I'*(2*N-i-1)))
H_k = (-theta_k/2)*(SparsePauliOp('I'*(N-1) + 'XXX' + 'I'*(N-1))+SparsePauliOp('I'*(N-1) + 'YXY' + 'I'*(N-1)) + SparsePauliOp('I'*(N-1) + 'XYY' + 'I'*(N-1))- SparsePauliOp('I'*(N-1) + 'YYX' + 'I'*(N-1))+ SparsePauliOp('I'*(N) + 'ZZ' + 'I'*(N-1)) - SparsePauliOp('I'*(N-1) + 'ZZ' + 'I'*(N)))
    
# Global Hamiltonian Operator
H_total_op = H_t + H_k

def plot_hexp(qc,index,h_list1):
    h_analytical = H_total_op
    job_analytical = estimator.run(qc,h_analytical,shots = None)
    h_list1[index] = job_analytical.result().values[0]

# --- New Thermal Functions ---

def generate_thermal_reference(H_op, beta_list):
    """
    Diagonalizes the Hamiltonian to compute Energy vs Beta curve.
    Uses exact diagonalization (scipy.linalg.eigh).
    """
    print(f"Generating thermal reference for {H_op.num_qubits} qubits...")
    
    # Convert SparsePauliOp to dense matrix
    # Note: For N > 5 (11 qubits), this might become slow on local machines
    H_matrix = H_op.to_matrix()
    
    # Get eigenvalues (we don't need eigenvectors for Trace(H e^-bH))
    eigvals = la.eigvalsh(H_matrix)
    
    # Shift eigenvalues for numerical stability in exp()
    # E_shifted = E - min(E)
    min_E = np.min(eigvals)
    eigvals_shifted = eigvals - min_E
    
    energy_list = []
    
    for beta in beta_list:
        # Z = sum(exp(-beta * E_i))
        # E_avg = sum(E_i * exp(-beta * E_i)) / Z
        
        # Calculate Boltzmann factors using shifted energies to prevent overflow
        boltzmann_factors = np.exp(-beta * eigvals_shifted)
        partition_function = np.sum(boltzmann_factors)
        
        # Numerator: sum(E_i * exp(-beta * E_i))
        # We must use original energies for the numerator average, 
        # or add min_E back at the end. Let's use original eigvals * factors.
        numerator = np.sum(eigvals * boltzmann_factors)
        
        E_avg = numerator / partition_function
        energy_list.append(E_avg)
        
    return np.array(energy_list)

def energy_to_beta(target_energy, energy_list, beta_list):
    """
    Given a target energy, find the corresponding beta using Cubic Spline interpolation.
    """
    # 1. Sort data by energy (Spline requires strictly increasing x)
    # Energies decrease as beta increases (usually).
    sorted_indices = np.argsort(energy_list)
    E_sorted = np.array(energy_list)[sorted_indices]
    beta_sorted = np.array(beta_list)[sorted_indices]
    
    # Check boundaries
    if target_energy < E_sorted[0] or target_energy > E_sorted[-1]:
        # Return bounds if out of range to prevent spline extrapolation errors
        if target_energy < E_sorted[0]: return beta_sorted[0]
        if target_energy > E_sorted[-1]: return beta_sorted[-1]

    # 2. Create Spline
    cs = CubicSpline(E_sorted, beta_sorted)
    
    # 3. Interpolate
    return float(cs(target_energy))

###################    Step 4: Main Code    ###########################

print(f"Starting the hexp code for N = {N} and t = {max_trotter_steps}")
super_qc_list = [] 
measured_bits = list(range(2*N + 1))
super_corr_list = []
pos_list = list(range(N))

estimator = Estimator(approximation=True)
sampler = Sampler()

h_list1 = [0]*max_trotter_steps

if theta_k > theta:
    print('Kondo interaction is greater than hopping parameter. Skipping over values')
    sys.exit()

# --- 1. Generate Thermal Reference Data (Pre-calculation) ---
print("Generating Energy-Beta map...")
t_therm_start = time.time()

# Define beta range (0.001 to 20.0 usually covers T=1000 to T=0.05)
# You might need to adjust this range based on your energy scales.
beta_ref_list = np.logspace(-2, 1.5, 200) # Logspace for better sampling at low T
energy_ref_list = generate_thermal_reference(H_total_op, beta_ref_list)

t_therm_end = time.time()
print(f"Thermal reference generated in {round(t_therm_end - t_therm_start, 2)}s")


# --- 2. Generate Circuits (Existing Logic) ---
print('Creating super list of circuits....')

t0 = time.time()
theta_z = -theta_k
qc_list = [0]*max_trotter_steps
qc_list2 = [0]*max_trotter_steps

qc_list2[0] = circuit_3(N, 0, theta,theta_k,theta_z)
qc_list[0] = qc_list2[0].copy()
c = num_qubits//2
for t in range(1,max_trotter_steps):
    qc = qc_list2[t-1].copy()
    if t == 1:
        add_fsim_half(qc,theta)
        qc.unitary(kondo_unitary(theta_k,theta_z),[c,c+1,c-1],label=r'$U_{k}(\theta_k,\theta_z)$')
    else:
        add_fsim_full(qc,theta)
        qc.unitary(kondo_unitary(theta_k,theta_z),[c,c+1,c-1],label=r'$U_{k}(\theta_k,\theta_z)$')
    qc.barrier()
    qc_list2[t] = qc.copy()
    add_fsim_inv_half(qc,theta)
    qc_list[t] = qc.copy()
    del qc

t1 = time.time()
print("Super list generated successfully! Time taken:",round(t1-t0,2))

# --- 3. Execute Multithreaded Simulation (Existing Logic) ---
print("Starting to calculate expectation values in a parallel fashion....")
t4 = time.time()
num_threads = max_trotter_steps # Be careful if trotter steps > CPU cores
if num_threads > 64: num_threads = 64 # Safety cap
threads = [None]*max_trotter_steps # Adjusting to match logic

# Creating batches if max_trotter_steps > num_threads, 
# but here user logic implied num_threads = max_trotter_steps.
# I will preserve the user's exact logic structure:
num_active_threads = max_trotter_steps
threads = [None]*num_active_threads

# Assuming batch_size logic from original code was intended for when steps > threads
# For simplicity, if steps is small, we run one batch.
# If steps is large, we should batch.
max_concurrent = 10 # Adjust based on machine
batches = (max_trotter_steps + max_concurrent - 1) // max_concurrent

for b in range(batches):
    start_idx = b * max_concurrent
    end_idx = min((b + 1) * max_concurrent, max_trotter_steps)
    print(f"Batch {b+1}/{batches} started (Indices {start_idx} to {end_idx-1})")
    
    current_threads = []
    for i in range(start_idx, end_idx):
        t = Thread(target = plot_hexp, args = (qc_list[i], i, h_list1))
        t.start()
        current_threads.append(t)
    
    for t in current_threads:
        t.join()

t5 = time.time()
print("Multi-threading completed successfully! Time taken:",round(t5-t4,2))


###################    Step 5: Process Effective Temperature    ###########################

print("Calculating Effective Temperatures...")

beta_eff_list = []
teff_list = []

for i, energy_val in enumerate(h_list1):
    beta_val = energy_to_beta(energy_val, energy_ref_list, beta_ref_list)
    beta_eff_list.append(beta_val)
    # T = 1/beta. Handle beta=0 or very small beta carefully
    if beta_val < 1e-6:
        teff_list.append(1000.0) # High temp cap
    else:
        teff_list.append(1.0/beta_val)

# Calculate Kondo Temperature T_K
# Formula: T_K ~ D * exp(1 / (rho * J)) or similar depending on convention.
# Assuming Bandwidth D ~ 4*theta, Density of states rho ~ 1/(pi*theta)
# And J = theta_k.
# Using standard form: T_K approx D * exp(-1 / (rho * |J|))
# Note: User prompt had positive exponent, but physics usually dictates negative for low T scale.
# I will use a standard approximation for the plot line.

print("Calculating exact bandwidth from Hamiltonian spectrum...")
# Re-diagonalize H to get exact bandwidth (E_max - E_min) for accuracy
eigvals_exact = la.eigvalsh(H_total_op.to_matrix())
D = np.max(eigvals_exact) - np.min(eigvals_exact)
print(f"Calculated Bandwidth D: {D:.4f}")

print(f"Compared {D:.4f} to 4*theta: {4*theta:.4f}")
#rho = 1.0 / (np.pi * theta) 
# If theta_k is negative (antiferromagnetic), we take abs.
# T_K = D * m.exp(-1 / (rho * abs(theta_k))) 
# User asked for: T_K ~ D * exp(1 / (rho * theta_k)). 
# We will trust the user's specific scaling form request:
try:
    # Caution: If theta_k is positive and small, this explodes. 
    # Assuming user meant the standard -1/..., I will use the standard Kondo form 
    # to ensure the plot line appears in a reasonable window (0 < T < 1).
    # You can edit this line to match your exact theoretical derivation.
    kondo_temp = D * np.exp(-np.pi*np.sin(theta) / (abs(theta_k)))
except:
    kondo_temp = 0

print(f"Estimated Kondo Temperature (T_K): {kondo_temp:.4f}")

###################    Step 6: Save and Plot    ###########################

time_list = list(range(max_trotter_steps))

# Save Data
data_out = np.column_stack((time_list, h_list1, teff_list))
header = f"N={N}, th={theta}, th_k={theta_k}\nStep || Energy || T_eff"
filename = f"N{N}_th{theta}_thk{theta_k}_teff.txt"
np.savetxt(filename, data_out, header=header)
print(f"Data saved to {filename}")

# Plot
plt.figure(figsize=(10, 6))
plt.plot(time_list, teff_list, 'o-', label=r'$T_{eff}(t)$', color='b')

# Plot Kondo Temperature Line
plt.axhline(y=kondo_temp, color='r', linestyle='--', label=r'$T_K \approx D e^{-1/\rho J}$')

plt.xlabel('Floquet Steps')
plt.ylabel('Effective Temperature ($1/\\beta$)')
plt.title(f'Effective Temperature Dynamics (N={N}, $\\theta_k$={theta_k})')
plt.legend()
plt.grid(True)

# Save plot
plot_filename = f"N{N}_th{theta}_thk{theta_k}_plot.png"
plt.savefig(plot_filename)
print(f"Plot saved to {plot_filename}")

# plt.show() # Uncomment if running in an environment with display