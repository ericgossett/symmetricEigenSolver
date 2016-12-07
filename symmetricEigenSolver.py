import numpy as np
import time
from scipy import linalg, random


def eigenCheck(A, values, vectors, show_message=False):
    """
    Check if the eigenpairs are commensurate with the system.
    Determine if ||Ax -lx||_2 approaches zero.
    """
    A_orig = A.copy()
    all_small = True
    for i in range(0, len(values)):
        norm_diff = np.linalg.norm(
            A_orig.dot(vectors[0:len(values), i]) -
            values[i]*vectors[0:len(values), i]
        )
        if norm_diff > 1e-6:
            all_small = False
    if all_small is True:
        if show_message:
            print 'Eigenvectors consistent with original matrix.'
        return True
    else:
        if show_message:
            print 'WARNING: Eigenvectors NOT consistent with original matrix.'
        return False


def checkOrthonormality(A, show_message=False):
    """
    Check if the eigenvectors are mutually orthogonal and that
    each eigenvector is normalized to unity (orthonormality).
    """
    tol = 1e-4
    consistent = True
    for i in range(0, A.shape[1]):
        for j in range(0, A.shape[1]):
            inner_product = np.dot(
                np.conjugate(A[0:A.shape[0], i].T),
                A[0:A.shape[0], j]
            )
            if i > j:
                if inner_product.real > tol or inner_product.imag > tol:
                    if show_message:
                        msg = 'ERROR: Inner product is not zero;'
                        msg += ' not orthogonal'
                        print msg
                        print np.dot(
                            np.conjugate(A[0:A.shape[0], i].T),
                            A[0:A.shape[0], j]
                        )
                    consistent = False
            elif i == j:
                if inner_product.real-1.0 > tol or inner_product.imag > tol:
                    if show_message:
                        msg = 'ERROR: Inner product is not zero;'
                        msg += ' not normalized to unity'
                        print msg
                        print np.dot(
                            np.conjugate(A[0:A.shape[0], i].T),
                            A[0:A.shape[0], j]
                        )
                    consistent = False
    if consistent is True:
        if show_message:
            print 'Eigenvectors orthogonal!'
    if consistent is False:
        if show_message:
            print 'WARNING: Eigenvectors NOT orthogonal!'


def CholeskyDecomposition(A):
    """
    Computes the Cholesky factorization of a matrix A.
    Returns the upper triangular matrix R (alpha in this
    case).  Note: The matrix must be a Hermitian positive
    definite matrix.
    """
    D_new = A.copy()
    alpha = np.zeros(A.shape)
    for i in range(0, A.shape[1]):
        if A[i][i] > 0:
            alpha[i][i] = np.sqrt(D_new[i][i])
        else:
            print ':('
        for j in range(i+1, A.shape[1]):
            alpha[i][j] = D_new[i][j] / alpha[i][i]
            for k in range(i+1, j):
                D_new[k, j] = D_new[k, j] - (alpha[i][j] * alpha[i][k])
    return alpha


def householderGeneratorRankKVersion(A):
    """
    New blocked Householder method, discussed by Auckenthaler.
    Takes in a sub-matrix with k-columns and determines the
    corresponding factorization.  Returns Q and R of sub-matrix.
    """
    R = np.zeros(A.shape)
    Q = np.zeros(A.shape)
    I = np.identity(A.shape[0])
    A_prev = A.copy()
    D = np.dot(np.conjugate(A.T), A)
    alphas = np.linalg.cholesky(D).T
    A_new = np.zeros(A.shape)
    y = np.zeros((A.shape[0], 1))
    for i in range(A.shape[1]):
        sign = 1.0
        if A_prev[i][i].real != 0.0:
            sign = A_prev[i][i].real / abs(A_prev[i][i].real)
            beta = alphas[i][i] * sign
        else:
            beta = np.linalg.norm(A_prev[i:A.shape[0], i]) * sign
        R[i][i] = -beta
        tau = (A_prev[i][i] + beta) / beta
        for j in range(A.shape[1]):
            rho = (
                A_prev[i][j] + ((alphas[i][i] * alphas[i][j]) / beta)
            ) / (A_prev[i][i] + beta)
            for ii in range(i+1, A.shape[0]):
                A_new[ii][j] = A_prev[ii][j] - (rho * A_prev[ii][i])
            R[i][j] = -sign * alphas[i][j]
        for ii in range(A.shape[0]):
            if ii < i:
                y[ii] = 0
            elif ii == i:
                y[ii] = 1.0
            else:
                y[ii] = A_prev[ii][i] / (A_prev[i][i] + beta)
        A_prev = A_new.copy()
        if i == 0:
            Q = I - np.dot(y, (tau*np.conjugate(y.T)))
        else:
            Q_temp = Q.copy()
            Q = np.dot(Q_temp, (I - np.dot(y, (tau*np.conjugate(y.T)))))
    return Q, R


def symmetricBandReduction(A, b, store_Q=True):
    """
    Transforms a matrix into a symmetric banded form with
    width b.  Option: Store Q in order to know the trasformation
    between A and the symmetric banded matrix S.
    """
    n = A.shape[1]
    if store_Q is True:
        Q_total = np.identity(A.shape[1])
    col = 0
    for block_col in range(0, (n/b)-1):
        col = block_col * b
        Q, R = householderGeneratorRankKVersion(A[col+b:n, col:col+b])
        if np.allclose(A[col+b:n, col:col+b], np.dot(Q, R)) is False:
            print 'ERROR: A != QR'
            print 'A', A[col+b:n, col:col+b]
            print 'QR', np.dot(Q, R)
            exit()
        if store_Q is True:
            if block_col == 0:
                Q_total[col+b:n, col+b:n] = Q
            else:
                temp = np.identity(A.shape[1])
                temp[col+b:n, col+b:n] = Q
                Q_total = np.dot(temp, Q_total)
        A[col+b:n, col:col+b] = R
        A[col+b:n, col+b:n] = np.dot(
            np.dot(np.conjugate(Q).T, A[col+b:n, col+b:n]),
            Q
        )
        if store_Q is False:
            Q, R = householderGeneratorRankKVersion(A[col+b:n, col:col+b])
        A[col:col+b, col+b:n] = R.T
    if store_Q is True:
        Q_total = np.conjugate(Q_total.T)
        return Q_total, A
    return A


def symToTridiag(A, b):
    """
    Turns A into a tridiagonal matrix. First, transforms A into
    a symmetric banded form, then into its tridiagonal form.
    """
    A_new = A.copy()
    try:
        Q1, S1 = symmetricBandReduction(A, b)
        print 'symmetric banded:'
        print S1
        Q2, T = symmetricBandReduction(S1, 1)
        print 'tridiagonal:'
        print T
        return np.dot(Q1.T, Q2), T
    except np.linalg.LinAlgError:
        b = b - 1
        return symToTridiag(A, b)


def QRAlgorithm(A, convergence=True):
    """
    Pure QR algorithm to determine the corresponding eigenpairs
    of the system.  Uses the eigenCheck function to determine
    convergence of the eigenpairs, or until max iterations reached.
    """
    if convergence is True:
        max_iter = 100000
    elif convergence is False:
        max_iter = 25
    A_orig = A.copy()
    zero_tol = 1.0e-7
    eigen_vec = np.identity(A.shape[1])
    i = 0
    while eigenCheck(A_orig, A.diagonal(), eigen_vec) is False:
        if i > max_iter:
            break
        i += 1
        A[np.abs(A) < zero_tol] = 0
        Q, R = np.linalg.qr(A)
        A = np.dot(R, Q)
        eigen_vec = np.dot(eigen_vec, Q)
    return A.diagonal(), eigen_vec


def getEig(A):
    """
    Calls functions to obtain the eigenpairs of a matrix A.
    First, call routine to determine tridiagonal form of matrix,
    then determine the corresponding eigenpairs using the QR method.
    Returns the eigenvalues and the eigenvectors after a back-transformation
    corresponding to the original matrix A.
    """
    Q, S = symToTridiag(A, A.shape[0]/2)
    v, w = QRAlgorithm(S)
    return v, np.dot(Q, w)
