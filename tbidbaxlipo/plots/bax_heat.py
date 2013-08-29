from tbidbaxlipo.util import fitting
import numpy as np
from matplotlib import pyplot as plt

def two_exp_func(t, fmax, k1, k2):
    return (fmax * (1 - np.exp(-k1 * (1 - np.exp(-k2*t)) * t)))

def fit_timecourse(time, y):
    fmax = fitting.Parameter(0.9)
    k1 = fitting.Parameter(0.1)
    k2 = fitting.Parameter(0.001)

    def two_exp_closure(t):
        return two_exp_func(t, fmax(), k1(), k2())

    fitting.fit(two_exp_closure, [fmax, k1, k2], y, time)

    return (fmax(), k1(), k2())

def fit_from_solver_sims(t, concs, simulations):
    """Fit a matrix of simulated dye release curves with a uniform
    time vector for all simulations, e.g., as produced by a series of
    deterministic simulations. The concentration index should be the
    first index into the simulations array, and the time index
    should be the second one.

    Returns
    -------
    list of numpy.array
        List with three elements (fmax_arr, k1_arr, k2_arr). The elements
        in the arrays represent the fmax, k1, and k2 values for each of
        the provided concentrations/timecourses.
    """
    num_concs = len(concs)
    fmax_arr = np.zeros(num_concs)
    k1_arr = np.zeros(num_concs)
    k2_arr = np.zeros(num_concs)
    for i in range(num_concs):
        (fmax_arr[i], k1_arr[i], k2_arr[i]) =  \
                        fit_timecourse(t, simulations[i, :])
    return (fmax_arr, k1_arr, k2_arr)

def fit_from_CptDataset(data):
    """Fit an HDF5 dataset of multi-compartment stochastic simulations."""
    time = data.sim_data[0,0,0,:]
    num_concs = data.sim_data.shape[0]
    fmax_arr = np.zeros(num_concs)
    k1_arr = np.zeros(num_concs)
    k2_arr = np.zeros(num_concs)
    for i in range(num_concs):
        (dr_mean, dr_sd) = data.get_mean_dye_release(i)
        (fmax_arr[i], k1_arr[i], k2_arr[i]) = fit_timecourse(time, dr_mean)
    return (fmax_arr, k1_arr, k2_arr)

def plot_fits_from_solver_sims(t, concs, simulations):
    plt.figure()
    (fmax_arr, k1_arr, k2_arr) = fit_from_solver_sims(t, concs, simulations)
    for i in range(len(concs)):
        # Plot data
        plt.plot(t, simulations[i,:], color='r')
        # Plot fit
        plt.plot(t, two_exp_func(t, fmax_arr[i], k1_arr[i], k2_arr[i]),
                 color='b')
    plt.show()
    return (fmax_arr, k1_arr, k2_arr)

def plot_fits_from_CptDataset(data):
    plt.figure()
    time = data.sim_data[0,0,0,:]
    num_concs = data.sim_data.shape[0]
    (fmax_arr, k1_arr, k2_arr) = fit_from_CptDataset(data)
    for i in range(num_concs):
        (dr_mean, dr_sd) = data.get_mean_dye_release(i)
        plt.plot(time, dr_mean, color='r')
        plt.plot(time, two_exp_func(time, fmax_arr[i], k1_arr[i], k2_arr[i]),
                 color='b')
    plt.show()
    return (fmax_arr, k1_arr, k2_arr)
