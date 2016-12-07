import numpy as np
from scipy import random


def randSymMat(size):
    """
    Generates a random symmetric matrix of size M x M.
    """
    A = random.rand(size, size)
    Q, R = np.linalg.qr(A)
    v = random.rand(size)
    D = np.diag(v)
    return np.dot(Q, np.dot(D, Q.T))


def isPosDef(x):
    """
    Check if matrix is postive definite.
    """
    return np.all(np.linalg.eigvals(x) > 0)
