"""
This module contains code to fit dye release titration fits to a variety of
mathematical functions. The fitting routines can handle data produced from
multiple sources, including actual experimental data as well as synthetic data
from deterministic and stochastic simulation.

Most fitting code is implemented in the superclass :py:class:`TitrationFit`;
specific mathematical functions are implemented in the various subclasses.  The
subclasses implement the specific mathematical function as a static method, and
also specify the number of parameters in the model, the names of the
parameters, and their initial values.
"""

from tbidbaxlipo.util import fitting
import numpy as np
from matplotlib import pyplot as plt

class TitrationFit(object):
    """Superclass for fitting kinetic titrations using mathematical functions.

    Implementing subclasses must implement a ``@staticmethod`` called
    ``fit_func`` that takes the time vector and a list (or array) of parameter
    values.

    In addition, subclasses should implement an ``__init__`` method that calls
    the superclass ``__init__`` with a list of parameter names and initial
    guesses for the parameter values, in an order corresponding the parameters
    used in ``fit_func``.

    Parameters
    ----------
    param_names : list of strings
        List of strings giving human-readable names for the parameters used
        in the model. The order of the names should correspond to the order
        of the initial guesses.
    initial_guesses : list of numbers
        Values to be used as the starting guesses for fitting.

    Attributes
    ----------
    k_arr : numpy.array
        Array containing the fitted parameters. The array has shape
        (num_params, num_concentrations), i.e., there is a set of parameters
        from every fitted concentration timecourse.
    concs : numpy.array
        The list of concentrations used in the titration. Serves as the
        x-coordinate when plotting the parameter values as a function of
        concentration.
    num_params : int
        The number of parameters in the function used for fitting.
    initial_guesses : list
        Initial values used for fitting.
    """
    def __init__(self, param_names, initial_guesses):
        if initial_guesses is None or len(initial_guesses) == 0:
            raise ValueError('initial_guesses cannot be None or empty.')
        if len(param_names) != len(initial_guesses):
            raise ValueError('There must be the same number of parameter names '
                             'and initial guesses.')

        self.param_names = param_names
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
                  for i in range(len(self.initial_guesses))]
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
        self.k_arr = np.zeros((len(self.initial_guesses), len(df.columns)))
        self.concs = df.columns.values
        for i, bax_conc in enumerate(self.concs):
            conc_data = df[bax_conc]
            time = conc_data[:, 'TIME']
            y = conc_data[:, 'MEAN']
            self.k_arr[:,i] = self.fit_timecourse(time, y)
        return self.k_arr

    def fit_from_CptDataset(self, data):
        """Fit an HDF5 dataset of multi-compartment stochastic simulations.

        Parameters
        ----------
        data : tbidbaxlipo.models.simulation.CptDataset
            Wrapper around an HDF5 dataset containing stochastic simulation
            results.

        Returns
        -------
        list of numbers
            The two-dimensional list of best-fit parameters, with
            shape (num_params, num_concs).
        """
        self.k_arr = np.zeros((len(self.initial_guesses),
                               data.sim_data.shape[0]))
        time = data.sim_data[0,0,0,:]
        for i in range(data.sim_data.shape[0]):
            (dr_mean, dr_sd) = data.get_mean_dye_release(i)
            self.k_arr[:, i] = self.fit_timecourse(time, dr_mean)
        return self.k_arr


    def plot_fits_from_dataframe(self, df):
        """Creates figures showing the best fits and the parameters.

        Includes a figure including both the original timecourse data and the
        best fits of the given model to the data; also creates figures showing
        the concentration-dependence of each of the parameters in the model.

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

        for k_index in range(len(self.initial_guesses)):
            plt.figure()
            plt.plot(self.concs, self.k_arr[k_index], marker='o')
            plt.xlabel('Bax (nM)')
            plt.ylabel(self.param_names[k_index])
            plt.title('%s vs. Bax conc' % self.param_names[k_index])
            plt.show()

    def plot_fits_from_CptDataset(self, jobs, data):
        """Creates figures showing the best fits and the parameters.

        Includes a figure including both the original timecourse data and the
        best fits of the given model to the data; also creates figures showing
        the concentration-dependence of each of the parameters in the model.

        Parameters
        ----------
        jobs : list of tbidbaxlipo.models.simulation.Job objects
            The job objects are expected to contain entries for the
            ``Bax_0`` parameter setting the initial condition of Bax
            in their ``params_dict`` field.
        data : tbidbaxlipo.models.simulation.CptDataset
            Wrapper around an HDF5 dataset containing stochastic simulation
            results.

        Returns
        -------
        list of numbers
            The two-dimensional list of best-fit parameters, with
            shape (num_params, num_concs).
        """
        plt.figure()
        self.fit_from_CptDataset(data)
        time = data.sim_data[0,0,0,:]
        assert len(jobs) == data.sim_data.shape[0]
        for i, job in enumerate(jobs):
            (dr_mean, dr_sd) = data.get_mean_dye_release(i)
            plt.plot(time, dr_mean, color='r')
            plt.plot(time, self.fit_func(time, self.k_arr[:,i]), color='b')
        plt.xlabel('Time (sec)')
        plt.show()

        self.concs = [j.params_dict['Bax_0'] for j in jobs]
        for k_index in range(len(self.initial_guesses)):
            plt.figure()
            plt.plot(self.concs, self.k_arr[k_index], marker='o')
            plt.xlabel(self.param_names[k_index])
            plt.ylabel('Bax (nM)')
            plt.title('%s vs. Bax conc' % self.param_names[k_index])
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
        super(OneExpNoFmax, self).__init__(param_names=['$k_1$'],
                                           initial_guesses=[1e-4])
    def fit_func(self, t, k_arr):
        """Single parameter exponential fitting function."""
        return (1 - np.exp(-k_arr[0]*t))

class Linear(TitrationFit):
    r"""Fit timecourses to a straight line through the origin.

    .. math::

        y(t) = k_1 * t

    Useful for fitting dye release data transformed into average pore
    numbers via the Poisson assumption of Schwarz.
    """
    def __init__(self):
        super(Linear, self).__init__(param_names=['$k_1$'],
                                     initial_guesses=[1e-4])

    def fit_func(self, t, k_arr):
        """Linear fitting function."""
        return k_arr[0]*t

class OneExpFmax(TitrationFit):
    r"""Fit timecourses to a two-parameter exponential (rate and max).

    .. math::

        y(t) = F_{max} (1 - e^{-k_1 t})
    """

    def __init__(self):
        super(OneExpFmax, self).__init__(param_names=['$k_1$', '$F_{max}$'],
                                         initial_guesses=[1e-4, 0.9])

    def fit_func(self, t, k_arr):
        """Two-parameter exponential fitting function."""
        return k_arr[1] * (1 - np.exp(-k_arr[0]*t))

class Schwarz(TitrationFit):
    r"""Fit timecourses to a three-parameter, two-phase linear function.

    .. math::

        y(t) = k_2 t + (k_1 - k_2) \frac{1 - e^{-k t}}{k}
    """
    def __init__(self):
        super(Schwarz, self).__init__(
                    param_names=['$k_1$', '$k_2$', '$tau$'],
                    initial_guesses=[1e-4, 1e-5, 8e-5])

    def fit_func(self, t, k_arr):
        """Fitting function from Schwarz."""
        return (k_arr[1] * t) + \
               (k_arr[0] - k_arr[1]) * \
               ((1 - np.exp(-k_arr[2] * t)) / k_arr[2])

class TwoExp(TitrationFit):
    r"""Fit timecourses to a three-parameter, two-phase exponential.

    .. math::

        y(t) = F_{max} \left(1 - e^{-k_1 (1 - e^{-k_2 t}) t} \right)
    """
    def __init__(self):
        super(TwoExp, self).__init__(
                    param_names=['$k_1$', '$F_{max}$', '$k_2$'],
                    initial_guesses=[1e-5, 0.5, 1e-3])

    def fit_func(self, t, k_arr):
        """Two-exponential fitting function."""
        return (k_arr[1]* (1 - np.exp(-k_arr[0] *
                                      (1 - np.exp(-k_arr[2]*t)) * t)))

class TwoExpLinear(TitrationFit):
    r"""Fit timecourses to a three-parameter, two-phase exponential with
        a linear term.

    .. math::

        y(t) = F_{max} \left(1 - e^{-k_1 (1 - e^{-k_2 t}) t} \right) + k_3 t
    """
    def __init__(self):
        super(TwoExpLinear, self).__init__(
                    param_names=['$k_1$', '$F_{max}$', '$k_2$', '$k_3$'],
                    initial_guesses=[1e-5, 0.5, 1e-3, 1e-6])

    def fit_func(self, t, k_arr):
        """Two-exponential fitting function."""
        return (k_arr[1]* (1 - np.exp(-k_arr[0] *
                                      (1 - np.exp(-k_arr[2]*t)) * t)) +
                (k_arr[3] * t))

class TwoExpNoFmax(TitrationFit):
    def __init__(self):
        super(TwoExpNoFmax, self).__init__(
                    param_names=['$k_1$', '$k_2$'],
                    initial_guesses=[1e-5, 1e-3])

    def fit_func(self, t, k_arr):
        """Two-exponential fitting function."""
        return (1 - np.exp(- k_arr[0] * (1 - np.exp(-k_arr[1]*t)) * t))


class TwoExpWithBackground(TitrationFit):
    r"""Fit timecourses to a three-parameter, two-phase exponential.

    .. math::

        y(t) = F_{max} \left(1 - e^{-k_1 (1 - e^{-k_2 t}) t} \right)
    """
    def __init__(self, bg_rates):
        super(TwoExpWithBackground, self).__init__(
                    param_names=['$k_1$', '$F_{max}$', '$k_2$'],
                    initial_guesses=[1e-4, 0.9, 1e-2])
        self.bg_rates = bg_rates

    def fit_func(self, t, k_arr):
        """Two-exponential fitting function."""
        return ((k_arr[1]+self.bg_rates[1]) *
                (1 - np.exp(-(self.bg_rates[0] + k_arr[0]) *
                             (1 - np.exp(-(self.bg_rates[2] + k_arr[2])
                                 * t)) * t)))

