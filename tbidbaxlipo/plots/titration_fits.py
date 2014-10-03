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
from tbidbaxlipo.util.plate_assay import TIME, VALUE

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
        (num_params, num_concentrations), i.e., there is a
        set of parameters from every fitted concentration timecourse.
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

class TwoExpNoT(TitrationFit):
    def __init__(self):
        super(TwoExpNoT, self).__init__(
                    param_names=['$k_1$', '$k_2$'],
                    initial_guesses=[1e-4, 1])

    def fit_func(self, t, k_arr):
        """Two-exponential fitting function."""
        return (1 - np.exp(-k_arr[1] * (1 - np.exp(- 9e-5 * t))))
        #return (1 - np.exp(-k_arr[1] * (1 - np.exp(- k_arr[0] * t))))

class LinkedEq(TitrationFit):
    def __init__(self, bax_conc):
        super(LinkedEq, self).__init__(
                    param_names=['$k_1$', '$k_2$', '$F_{max}$'],
                    initial_guesses=[6e-3, 1.25e-4, 0.5])
        self.bax_conc = bax_conc

    def fit_func(self, t, k_arr):
        """Two-exponential fitting function."""
        #return ((k_arr[2] / float(k_arr[0] - k_arr[1])) *
        #        (k_arr[0] * (1 - np.exp(-k_arr[1] * t)) -
        #         k_arr[1] * (1 - np.exp(-k_arr[0] * t))))
        return ((k_arr[2] / float(6e-3 - 1.25e-4)) *
                (6e-3 * (1 - np.exp(-1.25e-4 * t)) -
                 1.25e-4 * (1 - np.exp(-6e-3 * t))))

        #return (1 - np.exp(-k_arr[1] * (1 - np.exp(- k_arr[0] * t))))

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
                    initial_guesses=[1e-5, 1e-6, 2e-4])

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
                    initial_guesses=[2e-4, 0.5, 1e-2])

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

##########################################
# Plotting and analysis functions        #
##########################################

def plot_two_exp_fits(data, layout, num_reps=3, conc_str_index=1, plot=True):
    fmax_arr = np.zeros((num_reps, len(layout.keys())))
    k1_arr = np.zeros((num_reps, len(layout.keys())))
    k2_arr = np.zeros((num_reps, len(layout.keys())))
    conc_list = np.zeros(len(layout.keys()))

    for i, conc_str in enumerate(layout.keys()):
        if plot:
            plt.figure()
        wells = layout[conc_str]
        well_conc = float(conc_str.split(' ')[conc_str_index])
        conc_list[i] = well_conc

        for j, well in enumerate(wells):
            well_data = data[well]
            time = np.array(well_data[TIME])
            y = np.array(well_data[VALUE])

            fit = TwoExp()
            (k1, fmax, k2) = fit.fit_timecourse(time, y)

            fmax_arr[j, i] = fmax
            k1_arr[j, i] = k1
            k2_arr[j, i] = k2

            if plot:
                plt.plot(time, y, 'b')
                plt.plot(time, fit.fit_func(time, (k1, fmax, k2)), 'r')
                plt.title(well_conc)
                plt.xticks([0, 2000, 4000, 6000, 8000])
                plt.ylabel('% Release')
                plt.xlabel('Time (sec)')

    if plot:
        plt.show()
    return (fmax_arr, k1_arr, k2_arr, conc_list)

def plot_fmax_curve(fmax_arr, conc_list):
    fmax_means = np.mean(fmax_arr, axis=0)
    fmax_sds = np.std(fmax_arr, axis=0)

    # Plot of Fmax
    plt.figure()
    plt.ylabel(r'$F_{max}$ value (% Release)')
    plt.xlabel('[Bax] (nM)')
    plt.title(r'$F_{max}$')

    # Try fitting the fmax values to a hill function
    hill_vmax = fitting.Parameter(1)
    kd = fitting.Parameter(100)
    def hill(s):
        return ((hill_vmax() * s) / (kd() + s))
    fitting.fit(hill, [kd], fmax_means[1:], conc_list[1:])
    plt.plot(conc_list, hill(conc_list), 'r', linewidth=2, label='Hill fit')

    plt.text(35, 0.93, '$K_D$ = %.2f' % kd(), fontsize=16)
    #plt.text(35, 0.99, '$V_{max}$ = %.2f' % hill_vmax(), fontsize=16)

    # Try fitting the fmax values to an exponential function
    exp_fmax = fitting.Parameter(1.0)
    exp_k = fitting.Parameter(0.01)
    def exp_func(s):
        return exp_fmax() * (1 - np.exp(-exp_k()*s))
    fitting.fit(exp_func, [exp_fmax, exp_k], fmax_means[1:], conc_list[1:])
    plt.plot(conc_list, exp_func(conc_list), 'g', linewidth=2,
             label='Exp fit')
    print exp_fmax()
    print exp_k()

    # Try fitting the fmax values to a double-exponential
    exp_fmax1 = fitting.Parameter(1.)
    exp_fmax2 = fitting.Parameter(1.)
    exp_k1 = fitting.Parameter(0.001)
    exp_k2 = fitting.Parameter(0.001)
    def two_exp(s):
        return exp_fmax1()*(1 - np.exp(-exp_k1()*s)) + \
               exp_fmax2()*(1 - np.exp(-exp_k2()*s))
    fitting.fit(two_exp, [exp_fmax1, exp_fmax2, exp_k1, exp_k2],
                fmax_means[1:], conc_list[1:])
    plt.plot(conc_list, two_exp(conc_list), 'm', linewidth=2,
            label='Two-exp fit')

    # Try fitting the fmax values to a double-exponential
    hill_exp_fmax = fitting.Parameter(1.)
    hill_exp_kd = fitting.Parameter(100.)
    hill_exp_k = fitting.Parameter(0.01)
    def hill_exp(s):
        return hill_exp_fmax()*(1 - np.exp(
                             ((-hill_exp_k() * s) / (hill_exp_kd() + s)) * s))
    fitting.fit(hill_exp, [hill_exp_kd, hill_exp_k],
                fmax_means[1:], conc_list[1:])
    plt.plot(conc_list, hill_exp(conc_list), color='k', linewidth=2,
            label='Hill-exp fit')

    # Plot the data
    plt.errorbar(conc_list[1:], fmax_means[1:],
                 yerr=fmax_sds[1:] / np.sqrt(3),
                 label='Data', linestyle='', linewidth=2)

    plt.legend(loc='lower right')
    plt.show()

    # Pores per liposome at endpoint, hill fit
    plt.figure()
    pores_per_lipo = -np.log(1 - fmax_means[1:])
    plt.plot(conc_list[1:], pores_per_lipo, marker='o')
    plt.title('Pores per liposome at endpoint')

    exp_fmax = fitting.Parameter(1.0)
    exp_k = fitting.Parameter(0.01)
    m = fitting.Parameter(0.01)
    def exp_func(s):
        return exp_fmax() * (1 - np.exp(-exp_k()*s)) + m()*s
    fitting.fit(exp_func, [exp_fmax, exp_k, m], pores_per_lipo,
                conc_list[1:])
    plt.plot(conc_list, exp_func(conc_list), 'g', linewidth=2,
             label='Exp fit')
    print exp_fmax()
    print exp_k()
    print m()

    #hill_vmax = fitting.Parameter(1)
    #kd = fitting.Parameter(100)
    #def hill(s):
    #    return ((hill_vmax() * s) / (kd() + s))
    #fitting.fit(hill, [kd, hill_vmax], pores_per_lipo, conc_list[1:])
    #plt.plot(conc_list[1:], hill(conc_list[1:]), color='r')

def plot_k1_curve(k1_arr, conc_list):
    k1_means = np.mean(k1_arr, axis=0)
    k1_sds = np.std(k1_arr, axis=0)

    plt.figure()
    plt.errorbar(conc_list, k1_means, yerr= k1_sds / np.sqrt(3), color='r',
                 linestyle='', linewidth=2)
    plt.title('$k_1$')
    plt.xlabel('[Bax] (nM)')
    plt.ylabel('$k_1\ (\mathrm{sec}^{-1})$')
    """
    # Fit with exponential-linear curve
    vi = fitting.Parameter(0.05)
    vf = fitting.Parameter(0.025)
    tau = fitting.Parameter(0.02)
    # Define fitting function
    def biphasic(t):
        return (vf()*t) + ( (vi() - vf()) *
                            ((1 - np.exp(-tau()*t))/tau()) )
    fitting.fit(biphasic, [vi, vf, tau], k1_means, conc_list)
    plt.plot(conc_list, biphasic(conc_list), 'r', linewidth=2)

    # Fit with "double-binding curve"
    ksat = fitting.Parameter(0.0005)
    knonsat = fitting.Parameter(20)
    R0 = fitting.Parameter(20)
    def double_binding(t):
        return (R0() * t)/(ksat() + t) + (knonsat()*t)
    fitting.fit(double_binding, [R0, ksat, knonsat], k1_means, conc_list)
    plt.plot(conc_list, double_binding(conc_list), 'g', linewidth=2)

    plt.text(400, 0.00025,
            r'$\frac{R_0 Bax_0}{K_{sat} + Bax_0} + K_{nonsat} Bax_0$')
    plt.text(400, 0.00020, r'$R_0$ = %.3g nM' % R0())
    plt.text(400, 0.00015, r'$K_{sat}$ = %.3g nM' % ksat())
    plt.text(400, 0.00010, r'$K_{nonsat}$ = %.3g nM' % knonsat())

    plt.text(35, 0.00054,
            r'$v_f + (v_i - v_f) \left(\frac{1 - e^{-\tau t}}{\tau}\right)$')
    plt.text(35, 0.00049, r'$v_i$ = %.3g sec$^{-1}$ nM$^{-1}$' % vi())
    plt.text(35, 0.00044, r'$v_f$ = %.3g sec$^{-1}$ nM$^{-1}$' % vf())
    plt.text(35, 0.00039, r'$\frac{1}{\tau}$ = %.2f nM' % (1/tau()))
    plt.title('Biphasic fit of $k_1$')
    plt.show()
    """

def plot_k2_curve(k2_arr, conc_list):
    k2_means = np.mean(k2_arr, axis=0)
    k2_sds = np.std(k2_arr, axis=0)
    # Plot k2 on a linear scale
    plt.figure()
    plt.errorbar(conc_list, k2_means, yerr=k2_sds / np.sqrt(3))
    plt.title('$k_2$')


