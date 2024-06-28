import numpy as np
import sys

A = sys.argv[1]
B = 1./(1+A)

np.save("B.npy",B)
