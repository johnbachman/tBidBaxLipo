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
    """
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        err = y - function(x)
        return err

    if x is None: x = np.arange(y.shape[0])
    p = [param() for param in parameters]
    result = optimize.leastsq(f, p, ftol=1e-12, xtol=1e-12, maxfev=maxfev,
                              full_output=True)
    return result

def fit_matrix(function, parameters, y, x = None, maxfev=100000):
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
    residual_variance = np.var(f(p_final))

    return (residual_variance, result)


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

