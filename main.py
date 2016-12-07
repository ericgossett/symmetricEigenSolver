"""
A simple example calculating the eigenpairs for a random symmetric 10 x 10
matrix.
"""
import sys
import argparse
import symmetricEigenSolver as sym
from utils import matrixgenerators as matgen
import numpy as np

np.set_printoptions(linewidth=500)


A = matgen.randSymMat(10)

print 'inital matrix:'
print A

val, vec = sym.getEig(A)

print 'eigenvalues:'
print val

print 'eigenvectors:'
print vec
