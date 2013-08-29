from tbidbaxlipo.util import fitting
import numpy as np

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

