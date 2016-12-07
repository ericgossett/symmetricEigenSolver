"""
This is a routine used for benchmarking. This will store the times
of the tridiagonalization routine, determining the eigenpairs, and
the total time.

This will generate data in the DATA directory and plots in the PLOT
directory.

Note: Since our routine is not fully optimized, this benchmarking
procedure may take about an hour.

USAGE:
    python benchmark.py
"""
import sys
sys.path.append('../')
from utils import benchmarktools as bench
import numpy as np

bench.DetermineScaling()

# Plot Tridiagonalization times
n_size, times = np.loadtxt(
    'BENCHMARK_DATA/tridiagonalization_times.dat',
    unpack=True
)

bench.logPlotTimes(
    n_size,
    times,
    'Tridiagonalization Routine'
)

# Plot Eigenpair times
n_size, times = np.loadtxt(
    'BENCHMARK_DATA/eigenpair_times.dat',
    unpack=True
)

bench.logPlotTimes(
    n_size,
    times,
    'Eigenpair Routine'
)

# Plot Total times
n_size, times = np.loadtxt(
    'BENCHMARK_DATA/total_times.dat',
    unpack=True
)
bench.logPlotTimes(
    n_size,
    times,
    'Total Routine'
)
