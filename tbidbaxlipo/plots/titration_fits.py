from tbidbaxlipo.util import fitting
import numpy as np
from matplotlib import pyplot as plt

class TitrationFit(object):
    """Superclass for fitting kinetic titrations using mathematical functions.

    Implementing subclasses must implement a ``@staticmethod`` called
    ``fit_func`` that takes the time vector and a list (or array) of parameter
    values.

    In addition, subclasses should implement an ``__init__`` method that
    calls the superclass ``__init__`` with a list of initial guesses
    for the parameters, in an order corresponding the parameters used in
    ``fit_func``.

    Parameters
    ----------
    initial_guesses : list of numbers
        Values to be used as the starting guesses for fitting. The length
        of this array is used to determine the number of parameters.

    Attributes
    ----------
    k_arr : numpy.array
        Array containing the fitted parameters. The array has shape
        (num_params, num_concentrations), i.e., there is a set of
        parameters from every fitted concentration timecourse.
    concs : numpy.array
        The list of concentrations used in the titration. Serves as
        the x-coordinate when plotting the parameter values as a function
        of concentration.
    num_params : int
        The number of parameters in the function used for fitting.
    initial_guesses : list of numbers
        Initial values used for fitting.
    """
    def __init__(self, initial_guesses):
        if initial_guesses is None or len(initial_guesses) == 0:
            raise ValueError('initial_guesses cannot be None or empty.')

        self.k_arr = None
        self.concs = None
        self.num_params = len(initial_guesses)
        self.initial_guesses = initial_guesses

    def fit_timecourse(self, time, y):
        """Fit a single timecourse with the desired function.

        Parameters
        ----------
        time : numpy.array
            Vector of time values.
        y : numpy.array
            Vector of y-axis values (e.g., dye release).

        Returns
        -------
        list of numbers
            The best-fit parameter values produced by the fitting procedure.
        """
        params = [fitting.Parameter(self.initial_guesses[i])
                  for i in range(self.num_params)]
        def fit_func_closure(t):
            return self.fit_func(t, [p() for p in params])
        fitting.fit(fit_func_closure, params, y, time)
        return [p() for p in params]

    def fit_from_dataframe(self, df):
        """Fit all of the timecourses from the data in the pandas DataFrame.

        In addition to returning the parameters, sets the attributes
        ``self.concs`` and ``self.k_arr`` to the concentrations and fitted
        values as a side effect.

        Parameters
        ----------
        df : pandas.DataFrame
            Data is expected to be in the form returned by
            :py:func:`tbidbaxlipo.util.plate_assay.to_dataframe`, with the
            concentrations given as the columns, and the time, mean, and SD
            values given in a MultiIndex on the rows (with indices 'TIME',
            'MEAN', and 'SD', respectively).

        Returns
        -------
        list of numbers
            The two-dimensional list of best-fit parameters, with
            shape (num_params, num_concs).
        """
        self.k_arr = np.zeros((self.num_params, len(df.columns)))
        self.concs = df.columns.values
        for i, bax_conc in enumerate(self.concs):
            conc_data = df[bax_conc]
            time = conc_data[:, 'TIME']
            y = conc_data[:, 'MEAN']
            self.k_arr[:,i] = self.fit_timecourse(time, y)
        return self.k_arr

    def plot_fits_from_dataframe(self, df):
        """Creates a figure showing the timecourses and the best fits.

        Parameters
        ----------
        df : pandas.DataFrame
            Data is expected to be in the form returned by
            :py:func:`tbidbaxlipo.util.plate_assay.to_dataframe`.
        """
        plt.figure()
        self.fit_from_dataframe(df)
        for i, bax_conc in enumerate(df.columns):
            conc_data = df[bax_conc]
            time = conc_data[:, 'TIME']
            y = conc_data[:, 'MEAN']
            plt.plot(time, y, color='r')
            plt.plot(time, self.fit_func(time, self.k_arr[:,i]), color='b')
        plt.xlabel('Time (sec)')
        plt.show()

##########################################
# Various Fitting functions (subclasses) #
##########################################

class OneExpNoFmax(TitrationFit):
    r"""Fit timecourses to a one-parameter exponential function.

    .. math::

        y(t) = 1 - e^{-k_1 t}
    """
    def __init__(self):
        super(OneExpNoFmax, self).__init__(initial_guesses=[1e-4])

    @staticmethod
    def fit_func(t, k_arr):
        return (1 - np.exp(-k_arr[0]*t))

class Linear(TitrationFit):
    r"""Fit timecourses to a straight line through the origin.

    .. math::

        y(t) = k_1 * t

    Useful for fitting dye release data transformed into average pore
    numbers via the Poisson assumption of Schwarz.
    """
    def __init__(self):
        super(Linear, self).__init__(initial_guesses=[1e-4])

    @staticmethod
    def fit_func(t, k_arr):
        return k_arr[0]*t

class OneExpFmax(TitrationFit):
    r"""Fit timecourses to a two-parameter exponential (rate and max).

    .. math::

        y(t) = k_2 (1 - e^{-k1 t})
    """

    def __init__(self):
        super(OneExpFmax, self).__init__(initial_guesses=[1e-4, 0.9])

    @staticmethod
    def fit_func(t, k_arr):
        return k_arr[1] * (1 - np.exp(-k_arr[0]*t))

class TwoExp(TitrationFit):
    r"""Fit timecourses to a three-parameter, two-phase exponential.

    .. math::

        y(t) = k_2 \left(1 - e^{-k1 * (1 - e^{-k3 t}) t \right)
    """
    def __init__(self):
        super(TwoExp, self).__init__(initial_guesses=[1e-4, 0.9, 1e-2])

    @staticmethod
    def fit_func(t, k_arr):
        return (k_arr[1]* (1 - np.exp(-k_arr[0] *
                                      (1 - np.exp(-k_arr[2]*t)) * t)))

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

if __name__ == '__main__':
    from tbidbaxlipo.plots.layout_130905 import *
    #fit = OneExpNoFmax()
    #reset_norm_means = truncate_timecourses(reset_norm_means, 0, 80)
    #pores = truncate_timecourses(pores, 0, 80)
    #pore_df = to_dataframe(pores)
    df = to_dataframe(reset_norm_means)
    ion()
    fit = TwoExp()
    fit.plot_fits_from_dataframe(df)
    figure()
    plot(fit.concs, fit.k_arr[0], marker='o')
