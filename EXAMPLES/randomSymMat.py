"""
Calculates the eigenpairs for a randomly generated square symmetric
matrix A.

USEAGE:

    python randomSymMat.py

By default A is a 10 x 10 matrix. The size of A can be defined by including
the --dim flag:

    python randomSymMat.py --dim=100
"""
import sys
sys.path.append('../')
import argparse
import symmetricEigenSolver as sym
from utils import matrixgenerators as matgen
import numpy as np

np.set_printoptions(linewidth=500)
parser = argparse.ArgumentParser()
parser.add_argument(
    '--dim',
    type=int,
    help='Defines the dimension of the generated symmetric matrix.'
)

args = parser.parse_args()

dim = 10
if args.dim:
    dim = args.dim

A = matgen.randSymMat(dim)

print 'inital matrix:'
print A

val, vec = sym.getEig(A)

print 'eigenvalues:'
print val

print 'eigenvectors:'
print vec
