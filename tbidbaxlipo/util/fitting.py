r"""
A cookbook wrapper around the parameter optimization algorithm in
scipy.optimize. Drawn from
http://www.scipy.org/Cookbook/FittingData#head-5eba0779a34c07f5a596bbcf99dbc7886eac18e5.

This code is very useful for performing basic fitting of various mathematical
functions to data.

Example usage:

.. ipython:: python

    from tbidbaxlipo.util import fitting
    import numpy as np
    t = np.linspace(0, 100, 50)
    data = 3*t + 4
    m = fitting.Parameter(2)
    b = fitting.Parameter(5)
    def fit_func(x): return (m()*x + b())
    fitting.fit(fit_func, [m, b], data, t)
    print "m: %f\nb: %f" % (m(), b())
"""

import numpy as np
from scipy import optimize
from pysb.integrate import Solver
from matplotlib import pyplot as plt

class Parameter:
    """A simple object wrapper around parameter values.

    The parameter value is accessed using the overloaded call method, e.g.
    ``k()``.

    Parameters
    ----------
    value : number
        The initial guess for the parameter value.
    """
    def __init__(self, value):
            self.value = value

    def set(self, value):
        """Set the value."""
        self.value = value

    def __call__(self):
        """Get the value by calling the parameter, e.g. k()."""
        return self.value

def fit(function, parameters, y, x = None, maxfev=100000):
    """Fit the function to the data using the given parameters.

    Creates a wrapper around the fitting function and passes it to
    `scipy.optimize.leastsq`. Every evaluation of the wrapper function assigns
    new values (the values being tested at that evaluation) to the parameter
    objects in ``parameters`` as a side effect. This way, once
    `scipy.optimize.leastsq` is finished executing, the parameters have their
    new, optimized values.

    Parameters
    ----------
    function : function
        The fit function should itself be a closure, referencing the
        parameter values after they have been declared in the calling
        scope. The function should be able to operate on a input
        matching the argument ``x``, and produce an output that can
        be compared to the argument ``y``.
    parameters : list of :py:class:`Parameter`
        The parameters used in the fitting function.
    y : numpy.array
        The output of ``function`` evaluated on ``x`` will be compared to
        ``y``.
    x : numpy.array
        The input to ``function``, e.g., a time vector. If ``None`` (the
        default), an ordinal vector running from [0 ... y.shape[0]] is used.
    maxfev : int
        Maximum function evaluations. Passed to ``scipy.optimize.leastsq``
        as a keyword argument.

    Returns
    -------
    tuple: (residuals, result)
        The first entry in the tuple is a num_timepoints length array
        containing the residuals at the best fit parameter values. The second
        element is the result object returned by scipy.optimize.leastsq.
    """
    def f(params):
        i = 0
        for p in parameters:
            p.set(10 ** params[i])
            i += 1
        err = y - function(x)
        return err

    if x is None: x = np.arange(y.shape[0])
    p = [np.log10(param()) for param in parameters]
    result = optimize.leastsq(f, p, ftol=1e-12, xtol=1e-12, maxfev=maxfev,
                              full_output=True)
    # At this point the parameter instances should have their final fitted
    # values, NOT log transformed. To get the residuals at this final value,
    # we need to pass in log transformed values
    residuals = f([np.log10(param()) for param in parameters])

    return (residuals, result)

def fit_matrix(function, parameters, y, x = None, maxfev=100000):
    """Similar to fit, but takes fits a matrix of multiple replicates.

    The data array y is expected to be of the form
    (num_timepoints, num_replicates).

    Note that the residuals entry in the returned tuple has length
    num_timepoints * num_replicates.
    """

    (num_timepoints, num_replicates) = y.shape

    def f(params):
        err_matrix = np.zeros((num_timepoints, num_replicates))
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        # For every replicate, call the function
        for i in range(num_replicates):
            err_matrix[:, i] = y[:, i] - function(x)
        return err_matrix.reshape(num_timepoints * num_replicates)

    if x is None: x = np.arange(y.shape[0])
    p = [param() for param in parameters]
    result = optimize.leastsq(f, p, ftol=1e-12, xtol=1e-12, maxfev=maxfev,
                              full_output=True)
    p_final = [param() for param in parameters]
    residuals = f(p_final)

    return (residuals, result)


def mse(function, y, x):
    err = y - function(x)
    n = err.shape[0]
    sq_err = err**2
    return np.sum(sq_err) / n

def residuals(function, y, x):
    err = y - function(x)
    return err

def r_squared(y, yhat):
    ybar = np.mean(y)
    ssres = np.sum((y - yhat) ** 2)
    #ssreg = np.sum((yhat - ybar) ** 2)
    sstot = np.sum((y - ybar) ** 2)
    return 1 - (ssres / sstot)

class PySBFitResult(object):
    """Class to contain results of fitting PySB model by least squares."""
    def __init__(self, params, ypred, residuals, result):
        self.params = params
        self.ypred = ypred
        self.residuals = residuals
        self.result = result

def fit_pysb_builder(builder, obs_name, t, data,
                     obs_type='Expression'):
    """Use SciPy's leastsq algorithm to fit a PySB model to data.

    Parameters
    ----------
    builder : instance of pysb.builder.Builder
    obs_name : string
        Name of the observable or expression in the model that is to be fit
        to the data.
    t : numpy.array
        Vector containing the independent variable, e.g., time.
    data : numpy.array
        Vector of values to be fit.
    obs_type : string, "Expression" or "Observable"
        Indicates whether the named expression/observable specified by
        obs_name is to be found in the model's set of Expression objects
        or Observable objects.

    Returns
    -------
    Instance of PySBLeastSqFit containing results.
    """
    # Create lists of fitting parameters, one for all of the
    # models parameters, the other for just the ones we wish
    # to estimate
    all_params = []
    params_to_fit = []
    for p in builder.model.parameters:
        fit_p = Parameter(p.value)
        all_params.append(fit_p)
        if p in builder.estimate_params:
            params_to_fit.append(fit_p)

    # Function defining the model output
    s = Solver(builder.model, t)
    def fit_func(t):
        #import ipdb; ipdb.set_trace()
        param_values = [p() for p in all_params]
        s.run(param_values=param_values)
        if obs_type == 'Expression':
            return s.yexpr[obs_name]
        elif obs_type == 'Observable':
            return s.yobs[obs_name]
        else:
            raise ValueError('Unrecognized obs_type; must be '
                             '"Expression" or "Observable"')
    # The call to the fit function
    (residuals, result) = fit(fit_func, params_to_fit, data, t)
    pysb_fit = PySBFitResult([p() for p in all_params], fit_func(t),
                             residuals, result)
    return pysb_fit


class GlobalFit(object):
    """Least-squares fit of PySB model to a set of multiple timecourses, with a
    mix of globally and locally fit parameters.

    Parameters
    ----------
    builder : pysb.builder.Builder
        Builder containing the model to fit. Should contain an attribute
        builder.global_params for the parameters that are to be fit globally.
    time : np.array
        The time vector.
    data : list of np.array
        The experimental timecourses to fit.
    params : dict of lists
        The keys to the dict should be names of parameters in the PySB model
        (e.g., initial conditions); each value should be a list containing
        values for the parameter for each of the entries in the data list. The
        length of each value in the dict should match the length of the data
        list.
    obs_name : string
        The name of the model observable to compare against the data.
    """

    def __init__(self, builder, time, data, params, obs_name):
        # Check that the dimensions of everything that has been provided matches
        for tc in data:
            if not len(time) == len(tc):
                raise ValueError("Length of time vector must match the length "
                                 "of each data vector.")
        for p, vals in params.iteritems():
            if not len(vals) == len(data):
                raise ValueError("Each parameter in the params dict must have "
                                 "an entry for each entry in the data list.")

        self.builder = builder
        self.time = time
        self.data = data
        self.params = params
        self.obs_name = obs_name

    # A generic objective function
    def plot_func(self, x):
        s = Solver(self.builder.model, self.time)
        # Iterate over each entry in the data array
        for data_ix, data in enumerate(self.data):
            # Set the parameters appropriately for the simulation:
            # Iterate over the globally fit parameters
            for g_ix, p in enumerate(self.builder.global_params):
                p.value = x[g_ix]
            # Iterate over the locally fit parameters
            for l_ix, p in enumerate(self.builder.local_params):
                ix_offset = len(self.builder.global_params) + \
                            data_ix * len(self.builder.local_params)
                p.value = x[l_ix + ix_offset]
            # Now fill in the initial condition parameters
            for p_name, values in self.params.iteritems():
                p = self.builder.model.parameters[p_name]
                p.value = values[data_ix]
            # Now run the simulation
            s.run()
            # Plot the observable
            plt.plot(self.time, s.yobs[self.obs_name], color='g')

    def fit(self, maxfev=1e5):
        s = Solver(self.builder.model, self.time)

        # A generic objective function
        def obj_func(x):
            # The cumulative error over all timecourses
            err = 0
            # Iterate over each entry in the data array
            for data_ix, data in enumerate(self.data):
                # Set the parameters appropriately for the simulation:
                # Iterate over the globally fit parameters
                for g_ix, p in enumerate(self.builder.global_params):
                    p.value = x[g_ix]
                # Iterate over the locally fit parameters
                for l_ix, p in enumerate(self.builder.local_params):
                    ix_offset = len(self.builder.global_params) + \
                                data_ix * len(self.builder.local_params)
                    p.value = x[l_ix + ix_offset]
                # Now fill in the initial condition parameters
                for p_name, values in self.params.iteritems():
                    p = self.builder.model.parameters[p_name]
                    p.value = values[data_ix]
                # Now run the simulation
                s.run()
                # Calculate the squared error
                err += np.sum((data - s.yobs[self.obs_name]) ** 2)
            return err

        # Initialize the parameter array with initial values
        p = []
        for g_p in self.builder.global_params:
            p.append(g_p.value)
        for data_ix in range(len(self.data)):
            for l_p in self.builder.local_params:
                p.append(l_p.value)

        # Now, run the fit
        #result = optimize.leastsq(obj_func, p,
        #                    ftol=1e-12, xtol=1e-12, maxfev=maxfev,
        #                    full_output=True)
        result = optimize.minimize(obj_func, p, method='Nelder-Mead')

        return result

