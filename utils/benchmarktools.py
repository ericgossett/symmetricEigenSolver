import os
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
from scipy import linalg, random
from scipy.linalg import lapack
import symmetricEigenSolver as sym
from utils import matrixgenerators as matgen
import time


def DetermineScaling():
    """
    This routine determines the scaling with respect to the matrix size
    of the different parts of the routine. It tracks the time for the
    tridiagonal process, the eigenpair determination, and the total time.
    These values are stored in the "BENCHMARK_DATA" directory.
    """
    n_size = [
        50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 220,
        240, 260, 280, 300, 350, 400, 450, 500, 600, 700, 800,
        900, 1000
    ]

    for n in n_size:
        print 'Benchmarking: Finding eigenvectors for matrix of size ', n, '..'
        A = matgen.randSymMat(n)

        Q, T = sym.symToTridiag(A, A.shape[0])

        # Time symmetric band reduction into tridiagonal
        sym_start = time.time()
        S = sym.symmetricBandReduction(A, 1, False)
        sym_end = time.time()

        # If directory does not exist, make it
        if not os.path.exists('BENCHMARK_DATA'):
            os.makedirs('BENCHMARK_DATA')

        # Store sbr time
        data = open('BENCHMARK_DATA/tridiagonalization_times.dat', 'a')
        data.write(str(n) + '    ' + str(sym_end-sym_start) + '\n')
        data.close()
        print 'Tridiagonalization time = ', str(sym_end - sym_start)

        # Time eigenvalue solver
        eig_start = time.time()
        v, w = sym.QRAlgorithm(S, False)
        eig_end = time.time()

        # Store eig time
        data = open('BENCHMARK_DATA/eigenpair_times.dat', 'a')
        data.write(str(n) + '    ' + str(eig_end - eig_start) + '\n')
        data.close()
        print 'Eigenpair time = ', str(eig_end - eig_start)

        # Store total time
        data = open('BENCHMARK_DATA/total_times.dat', 'a')
        data.write(str(n) + '    ' + str(eig_end - sym_start) + '\n')
        data.close()
        print 'Total time = ', str(eig_end - sym_start)


def logPlotTimes(n_size, times, name):
    """
    Plots the scaling data in a log-log plot. It determines the
    slope of the line, indicating the scaling (square, cubic, etc.).
    The plots are saved in "BENCHMARK_DATA".
    """
    plt.style.use('ggplot')

    figure()
    log_times = [math.log(time) for time in times]
    log_size = [math.log(n) for n in n_size]
    plot(log_size, log_times, 'ro')

    m, b = np.polyfit(log_size, log_times, 1)
    print 'Slope:', m
    fit = np.polyfit(log_size, log_times, 1)
    fit_fn = np.poly1d(fit)
    x = np.arange(0, 8, 0.1)
    plot(x, fit_fn(x), '-k', label='slope='+str(m))

    legend(loc=4)
    ylabel('Time [sec]')
    xlabel('Matrix size [n]')
    title(name + ' Scaling (Log-Log)')
    savefig('BENCHMARK_DATA/'+name.replace(" ", "") + '_plot.png')
