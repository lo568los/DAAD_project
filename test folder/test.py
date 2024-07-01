import numpy as np
import sys

A = float(sys.argv[1])
B = 1./(1+A)

np.save("./Saved/A_{}.npy".format(A),B)
