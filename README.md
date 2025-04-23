## Introduction
This repository contains code used for writing the paper "Design and Benchmarks for Emulating Kondo Dynamics on a Quantum Chip" [arxivlink](https://arxiv.org/abs/2301.02146). The naming convention in the code closely matches the notation used in the paper. For any questions / comments, please feel free to contact me. I will try my best to get back to you.

## Requirements
1. Python with Qiskit (https://docs.quantum.ibm.com/guides), Numpy, Scipy installed.

## Code
1. **test_sz.py,test_h.py,test_ent.py** : Present in "scaled_codes", used for plotting impurity magnetization, heating and entanglement measures for N = 6 and N = 10 with a Fermi Sea state or a translationally invariant TS state using Qiskit statevector simulation
2. **randomized_states.ipynb** : Present in "playground ipynb files", used for plots with RTP (randomized tensor product) initial state.
3. **bruteforce_ed.ipynb** : Present in "playground ipynb files", used for ED for N = 6 sites. Main data and ED codes present in [this link](https://indianinstituteofscience-my.sharepoint.com/:f:/g/personal/ssoumyadeep_iisc_ac_in/EhedPZhYjm9Hqn5V3cHudbYBHcF4tXZTnnC9-7X3J-D3pw?e=OyCC9E).
4. **Coding_circuits.ipynb** : Present in "playground ipynb files", last part of the file used for finding the dispersion relation plots in Appendix A.


The rest of the code is self-explanatory, and consists mainly of plotting functions and helper/test functions.

