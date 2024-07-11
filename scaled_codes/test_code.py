## This is the test version of the python file that is passed on to the cluster

## Author: Soumyadeep Sarma


###################    Step 1: Import all the libraries    ###########################

import qiskit
import sys
from qiskit.quantum_info import Statevector
from qiskit import QuantumCircuit
import qiskit_aer 
from qiskit.quantum_info import state_fidelity
from qiskit_aer import AerSimulator
#from qiskit import transpile
#from qiskit.quantum_info.states.random import random_statevector
#from qiskit.circuit.library import Initialize
#from qiskit.visualization import plot_bloch_multivector
import numpy as np
from qiskit.quantum_info import partial_trace # To check later whether our derived density matrix is correct
from qiskit.quantum_info import DensityMatrix
from qiskit.quantum_info import purity
from scipy import linalg as la


from qiskit_aer.primitives import Sampler
from qiskit_aer.primitives import Estimator
from qiskit.quantum_info import SparsePauliOp
from qiskit.quantum_info import Operator

from itertools import combinations, cycle #Used for fermi state
import math as m
import cmath as cm


###################    Step 2: Initialize the system parameters (this format utilizes paralellization)   ###########################

N = int(sys.argv[1])  #Number of fermionic sites
theta = float(sys.argv[2]) #hoping parameter for free fermions
theta_k = float(sys.argv[3]) #Kondo interaction
max_trotter_steps = int(sys.argv[4]) #number of time steps
#time_corr = int(sys.argv[5]) #time for correlator functions

num_qubits = 2*N + 1  #In split side configuration

###################    Step 3: Define all the helper functions    ###########################

coeff_dict  = {}

def array_k1(num_qubits):  #assuming num_qubits is even
    array_k = []
    m = num_qubits/2
    if (m)%2!=0: #replace with m-1 for previous results
        for j in range(-int((m)//2),int((m)//2) + 1): #replace with m-1 for previous results
            array_k.append(2*np.pi*j/num_qubits)
    else:
        for j in range(-int((m)//2),int((m)//2)):
            array_k.append(2*np.pi*j/num_qubits)
    return array_k


def recursive_nested(l,num_qubits,coeff_array,coeff = 1,bitstr=''):

    m = int(num_qubits/2) # taking always even number of qubits

    #coeff_dict_2 = {}
    if l==m-1: # m-1, but we start with 0 indexing
        coeff_copy = coeff  #to ensure multiplied coeffs in previous rounds is preserved and reused
        bitstr_copy = bitstr
        for i in range(num_qubits):
            if str(i) in bitstr:
                pass
            else:
                coeff = coeff*coeff_array[l,i]
                bitstr = bitstr + f'{i}'
                bitstr_sorted = sort_bitstr(bitstr) #to sort the string first
                perm = perm_str2(bitstr,bitstr_sorted) # to compare hamming distance
                if bitstr_sorted in coeff_dict.keys():
                    coeff_dict[bitstr_sorted]+=coeff*perm
                else:
                    coeff_dict[bitstr_sorted] = coeff*perm
            bitstr = bitstr_copy
            coeff = coeff_copy
    if l!=m-1:
        coeff_copy = coeff  #to ensure multiplied coeffs in previous rounds is preserved and reused
        bitstr_copy = bitstr
        for i in range(num_qubits):
            if str(i) in bitstr:
                pass
            else:
                #print(coeff)
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
    #print(bit_array)
    bitstr_sorted = ''
    for i in range(len(bit_array)):
        bit_array[i] = int(bit_array[i])
    bit_array.sort()
    #print(bit_array)
    for i in bit_array:
        #print(i,str(i))
        bitstr_sorted += str(i)
        #print(bitstr_sorted)

    return bitstr_sorted

def perm_str2(cmpr,word):

   #  word = 'eyssaasse' base string
   # cmpr = 'seasysaes'  a string to find number of swaps from the base string
    swaps = 0

    # 1)
    chars = {c: [] for c in word}
    [chars[c].append(i) for i, c in enumerate(word)]
    for k in chars.keys():
        chars[k] = cycle(chars[k])

    # 2)
    idxs = [next(chars[c]) for c in cmpr]

    # 3)
    for cmb in combinations(idxs, 2):
        if cmb[0] > cmb[1]:
            swaps += 1

    #print(swaps)
    if swaps%2 == 0:
        return 1
    else:
        return -1
    
def fermi_state(num_qubits): #of the form num_qubits = 2*odd number. This function created the fermi state (FS)


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

    #print(coeff_array)

    coeff_dict_2 = recursive_nested(0,num_qubits,coeff_array)
    #print(coeff_dict_2)
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

    #print(bitstr_dict)
    

    fermi_state = Statevector([0]*(2**num_qubits))

    for bstr in bitstr_dict.keys():
        fermi_state += Statevector.from_label(bstr)*bitstr_dict[bstr]

    val_array = []
    for i in bitstr_dict.values():
        val_array.append(i)
    np_array = np.array(val_array)
    #print(np_array)
    #print(np.linalg.norm(np_array))
    
    fermi_state = fermi_state/np.linalg.norm(val_array)
    #print(fermi_state)
    fermi_state.is_valid()

    return fermi_state

def fermi_state_circuit(N,num_cl_bits = 0):  #Initialize circuit with FS with both spin chains having it's own FS
    qc = QuantumCircuit(num_qubits,num_cl_bits)
    fermi_state_up = fermi_state(N)
    fermi_state_down = fermi_state(N)
    qc.initialize(fermi_state_up,range(N))
    qc.initialize(fermi_state_down,range(N+1,2*N+1))
    return qc

def fermi_state_circuit2(num_qubits,num_cl_bits = 0):  #Initialize circuit with FS over all non-impurity sites
    qc = QuantumCircuit(num_qubits,num_cl_bits)
    fermi_state_up = fermi_state(int(num_qubits))
    qc.initialize(fermi_state_up,range(int(num_qubits)))
    return qc


def fermion_state(N,pos_list,num_cl_bits = 0):   #Initialize state with fixed number of paricles in up or down chain
    qc = QuantumCircuit(2*N+1,num_cl_bits)
    for i in range(2*N+1):
        if i in pos_list or i==N:
            continue
        else:
            qc.x(i)
    return qc

def fsim(theta,phi,beta):  #Block for free fermions
    fsim = Operator([[1,0,0,0],
                   [0,m.cos(theta),1j*cm.exp(1j*beta)*m.sin(theta),0],
                   [0,1j*cm.exp(-1j*beta)*m.sin(theta),m.cos(theta),0],
                   [0,0,0,cm.exp(1j*phi)]])
    return fsim

def add_fsim_half(qc,angles):
    theta = angles[0]
    phi = angles[1]
    beta = angles[2]

    fsim1 = fsim(theta,phi,beta)
    #Adding fsim in even layers
    for i in range(0,qc.num_qubits//2-1,2):
        qc.unitary(fsim1,[i,i+1],label = r'fsim$(\theta,\phi)$')

    for i in range(qc.num_qubits//2+1,qc.num_qubits-1,2):
        qc.unitary(fsim1,[i,i+1],label = r'fsim$(\theta,\phi)$')
        
    #Adding fsim in odd layers
    for i in range(1,qc.num_qubits//2-1,2):
        qc.unitary(fsim1,[i,i+1],label = r'fsim$(\theta,\phi)$')

    for i in range(qc.num_qubits//2+2,qc.num_qubits-1,2):
        qc.unitary(fsim1,[i,i+1],label = r'fsim$(\theta,\phi)$')

def add_fsim_inv_half(qc,angles):
    theta = angles[0]
    phi = angles[1]
    beta = angles[2]

    fsim1 = fsim(theta,phi,beta)

    #Adding fsim in odd layers
    for i in range(1,qc.num_qubits//2-1,2):
        qc.unitary(fsim1,[i,i+1],label = r'fsim$(\theta,\phi)$')

    for i in range(qc.num_qubits//2+2,qc.num_qubits-1,2):
        qc.unitary(fsim1,[i,i+1],label = r'fsim$(\theta,\phi)$')
        
    #Adding fsim in even layers
    for i in range(0,qc.num_qubits//2-1,2):
        qc.unitary(fsim1,[i,i+1],label = r'fsim$(\theta,\phi)$')

    for i in range(qc.num_qubits//2+1,qc.num_qubits-1,2):
        qc.unitary(fsim1,[i,i+1],label = r'fsim$(\theta,\phi)$')


    
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

    """kondo_unitary = Operator([[1,0,0,0,0,0,0,0],
                          [0,1,0,0,0,0,0,0],
                          [0,0,l1,0,0,0,0,0],
                          [0,0,0,c1*l2,1j*l2*s1,0,0,0],
                          [0,0,0,1j*l2*s1,c1*l2,0,0,0],
                          [0,0,0,0,0,l1,0,0],
                          [0,0,0,0,0,0,1,0],
                          [0,0,0,0,0,0,0,1]])"""
    kondo_unitary = Operator([[1,0,0,0,0,0,0,0],
                          [0,a,0,0,0,0,b,0],
                          [0,0,1,0,0,0,0,0],
                          [0,0,0,d,0,0,0,0],
                          [0,0,0,0,c,0,0,0],
                          [0,0,0,0,0,1,0,0],
                          [0,b,0,0,0,0,a_dag,0],
                          [0,0,0,0,0,0,0,1]])
    
    kondo_unitary_2 = Operator([[1,0,0,0,0,0,0,0],
                          [0,1,0,0,0,0,0,0],
                          [0,0,c1*l1,0,0,1j*s1*l1,0,0],
                          [0,0,0,l2,0,0,0,0],
                          [0,0,0,0,l2,0,0,0],
                          [0,0,1j*s1*l1,0,0,c1*l1,0,0],
                          [0,0,0,0,0,0,1,0],
                          [0,0,0,0,0,0,0,1]])
    
    return kondo_unitary_2

def circuit_3(N, pos_list, trotter_steps,angles = 0,theta_k = 0,theta_z = 0, num_cl_bits = 0, trotter_barriers = False, save = False):
    if num_cl_bits == 0:
        qc = fermi_state_circuit(N)
    else:
        qc = fermi_state_circuit(N,num_cl_bits)
    qc.x(N)
    qc.barrier()
    
    c = num_qubits//2
    for i in range(trotter_steps):
        add_fsim_half(qc,angles)
        qc.unitary(kondo_unitary(theta_k,theta_z),[c,c+1,c-1],label=r'$U_{k}(\theta_k,\theta_z)$')
        add_fsim_inv_half(qc,angles)
        if trotter_barriers:
            qc.barrier()
    if save == True:
        qc.save_statevector()
    #qc.save_statevector()  remove save for changing to operator
    return qc

def plot_mag_impurity(super_qc_list_20):
    for i in range(len(super_qc_list_20)):
        theta = super_qc_list_20[i][1]
        theta_k = super_qc_list_20[i][2]
        qc_list = super_qc_list_20[i][0]
        imp_observables = [SparsePauliOp('I'*N + 'Z' + 'I'*N)]*max_trotter_steps
        job_1 = estimator.run(qc_list,imp_observables,shots = None)
        expectation_list_1 = list(job_1.result().values)
        return expectation_list_1
    
H_t = 0
H_k = 0
for i in range(2*N):
    if i==N-1 or i==N:
        continue
    else:
        H_t += -theta*(SparsePauliOp('I'*(i) + 'XX' + 'I'*(2*N-i-1)) + SparsePauliOp('I'*(i) + 'YY' + 'I'*(2*N-i-1)))
H_k = (-theta_k/2)*(SparsePauliOp('I'*(N-1) + 'XXX' + 'I'*(N-1))+SparsePauliOp('I'*(N-1) + 'YXY' + 'I'*(N-1)) + SparsePauliOp('I'*(N-1) + 'XYY' + 'I'*(N-1))- SparsePauliOp('I'*(N-1) + 'YYX' + 'I'*(N-1))+ SparsePauliOp('I'*(N) + 'ZZ' + 'I'*(N-1)) - SparsePauliOp('I'*(N-1) + 'ZZ' + 'I'*(N)))
    
def plot_hexp(super_qc_list_50):


    for i in range(len(super_qc_list_50)):
        theta = super_qc_list_50[i][1]
        theta_k = super_qc_list_50[i][2]
        qc_list = super_qc_list_50[i][0]

        h_analytical = [H_t + H_k]*max_trotter_steps
        #print("Operator obtained from circuit")
        job_analytical = estimator.run(qc_list,h_analytical,shots = None)
        h_values1 = list(job_analytical.result().values)

        return h_values1
    
imp_op = SparsePauliOp('I'*N + 'Z' + 'I'*N)

def ferm_mag(pos):
    op1 = SparsePauliOp('I'*(N+pos) + 'Z' + 'I'*(N-pos))
    op2 = SparsePauliOp('I'*(N-pos) + 'Z' + 'I'*(N+pos))

    ferm_mag_op = 0.5*(op2 - op1)
    #print(ferm_mag_op)
    return ferm_mag_op

def correlator_expectation2(pos,qc_list):
    op1 = ferm_mag(pos)
    corr_op = op1 @ imp_op
    obs_list = [corr_op]*max_trotter_steps
    job = estimator.run(qc_list,obs_list,shots = None)
    exp_vals = list(job.result().values)
    return exp_vals

def reduced_corr(pos,qc_list):
    op1 = ferm_mag(pos)
    op2 = SparsePauliOp('I'*N + 'Z' + 'I'*N)
    obs_list1 = [op1]*max_trotter_steps
    obs_list2 = [op2]*max_trotter_steps
    job1 = estimator.run(qc_list,obs_list1,shots = None)
    job2 = estimator.run(qc_list,obs_list2,shots = None)
    exp_vals1 = list(job1.result().values)
    exp_vals2 = list(job2.result().values)
    exp_vals_red = [a*b for a,b in zip(exp_vals1,exp_vals2)]
    return exp_vals_red

def plot_correlator(qc_list,pos):
    exp_vals = correlator_expectation2(pos,qc_list)
    exp_vals_red = reduced_corr(pos,qc_list)
    final_vals = [a-b for a,b in zip(exp_vals,exp_vals_red)]
    return final_vals

def trace_norm(density_matrix):
    density_matrix = np.array(density_matrix)
    eigenvalues, eigenvectors = np.linalg.eig(density_matrix)
    sum = 0
    #print(eigenvalues)
    for i in eigenvalues:
        if i.real < 0:
            sum += i.real
    #print(sum)
    return abs(sum)

def calculate_entropy(rho):
    rho = np.array(rho)
    R = rho*(la.logm(rho))
    S = -np.matrix.trace(R)
    return S


def von_neumann_entropy(reduced_dm_list):
    entropy_list = []
    for dm in reduced_dm_list:
        entropy = calculate_entropy(dm)
        entropy_list.append(entropy)
    return entropy_list

def negativity(density_matrix_list):
    neg_list = []
    for density_matrix in density_matrix_list:
        dm_pt = density_matrix.partial_transpose([N])
        neg = trace_norm(dm_pt)
        #print(neg)
        neg_list.append(neg.real)
    return neg_list

def concurrence(reduced_dm_list):
    concurrence_list = []
    for dm in reduced_dm_list:
        c = purity(dm).real
        if c>1:
            c = 1
            print("Purity is greater than 1,by a value of:",c-1)
        concurrence_list.append(np.sqrt(2*(1-c)))
    return concurrence_list


###################    Step 4: The main code which generates <S^z-imp>, <H>(t) and correlator functions w.r.t time and space    ###########################

super_qc_list = []  #list to store circuits for each parameter combination
measured_bits =list(range(2*N + 1))  #list of qubits to measure
super_corr_list = []  #list to store correlator functions
pos_list = list(range(N) ) #list of positions to calculate correlator functions

estimator = Estimator(approximation=True) #estimator object to estimate the expectation values
sampler = Sampler()  #sampler object to sample the circuits

if theta_k > theta:
    print('Kondo interaction is greater than hopping parameter. Skipping over values')
else:
    print('Creating super list of circuits....')
    theta_z = -theta_k
    qc_list = []
    qc_list_2 = []
    for t in range(max_trotter_steps):
        qc = circuit_3(N,[0,1,2*N], t, [theta,0,0],theta_k,theta_z,num_cl_bits = len(measured_bits), trotter_barriers = True, save = True)
        qc.measure(measured_bits,list(range(len(measured_bits))))
        qc_list.append(qc)
    super_qc_list.append((qc_list,theta,theta_k))

    print("Super list genereted successfully!")
    print('Calculating avg impurity magnetization....')

    expectation_list_1 = plot_mag_impurity(super_qc_list)

    print('Calculating H(t)....')

    h_values1 = plot_hexp(super_qc_list)

    print("Calculating entanglement measures (Concurrence and Von-Neumann Entropy)")

    Z_imp_observables = [SparsePauliOp('I'*(N) + 'Z' + 'I'*(N))]*max_trotter_steps
    X_imp_observables = [SparsePauliOp('I'*(N) + 'X' + 'I'*(N))]*max_trotter_steps
    Y_imp_observables = [SparsePauliOp('I'*(N) + 'Y' + 'I'*(N))]*max_trotter_steps

    Z_avg = estimator.run(super_qc_list[0][0],Z_imp_observables,shots = None)
    X_avg = estimator.run(super_qc_list[0][0],X_imp_observables,shots = None)
    Y_avg = estimator.run(super_qc_list[0][0],Y_imp_observables,shots = None)

    Z_values = list(Z_avg.result().values)
    X_values = list(X_avg.result().values)
    Y_values = list(Y_avg.result().values)

    reduced_dm_tomo = []
    for i in range(max_trotter_steps):
        reduced_dm = DensityMatrix((1/2)*(np.eye(2) + X_values[i]*np.array([[0,1],[1,0]]) + Y_values[i]*np.array([[0,-1j],[1j,0]]) + Z_values[i]*np.array([[1,0],[0,-1]])))
        reduced_dm_tomo.append(reduced_dm)

    von_neumann_entropy_list = von_neumann_entropy(reduced_dm_tomo)
    von_neumann_entropy_list = [x.real for x in von_neumann_entropy_list]
    #negativity_list = negativity(density_matrix_list)    #Takes a tremendous amount of time.... best calculated in a separate script/cell
    concurrence_list = concurrence(reduced_dm_tomo)
    concurrence_list = [y.real for y in concurrence_list]

    



    """print('Calculating correlator functions as function of time and space....')
    for pos in pos_list:
        super_corr_list.append(plot_correlator(super_qc_list[0][0],pos))   #calculate correlator functions in a separate python file...."""

    print('All calculations done!')

###################    Step 5: Save the results in a file    ###########################

    time_list = list(range(max_trotter_steps))
    

    string = f"N = {N}, theta = {theta}, theta_k = {theta_k}, max_trotter_steps = {max_trotter_steps}"
    header_sz = string + "\n TIME || S_z impurity expectation value"
    header_h = string + "\n TIME || H(t) expectation value"
    header_ent = string + "\n TIME || CONCURRENCE || VON NEUMANN"

    data_sz = np.column_stack((time_list,expectation_list_1))
    np.savetxt(f"N = {N}, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}_sz.txt",data_sz,header = header_sz)
    data_h = np.column_stack((time_list,h_values1))
    np.savetxt(f"N = {N}, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}_h.txt",data_h,header = header_h)
    data_ent = np.column_stack((time_list,concurrence_list,von_neumann_entropy_list))
    np.savetxt(f"N = {N}, theta = {round(theta,2)}, theta_k = {round(theta_k,2)}_ent.txt",data_ent,header = header_ent)




    """header_corr = string + "\n TIME || Correlator function_vals"
    data_corr = np.column_stack((time_list,super_corr_list[0])) 
    for i in range(1,len(pos_list)):
        data_corr = np.column_stack((data_corr,super_corr_list[i]))

    np.savetxt(f"N = {N}, theta = {theta}, theta_k = {theta_k}_corr.txt",data_corr,header = header_corr)"""
    


       

  

    

