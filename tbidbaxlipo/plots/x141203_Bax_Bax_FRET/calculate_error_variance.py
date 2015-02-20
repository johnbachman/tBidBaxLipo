from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import fitting

def calc_err_var_linear(data, last_n_pts=50, plot=False):
    a1 = fitting.Parameter(0.005)
    b = fitting.Parameter(1.)
    def fit_func(t):
        return a1()*t + b()
    t = np.arange(len(data))
    data_subset = data[-last_n_pts:]
    t_subset = t[-last_n_pts:]
    result = fitting.fit(fit_func, [a1, b], data_subset, t_subset,
                         log_transform=False)
    residuals = result[0]
    if plot:
        plt.figure()
        plt.subplot(1, 2, 1)
        plt.plot(t, data)
        plt.plot(t_subset, fit_func(t_subset))
        plt.subplot(1, 2, 2)
        plt.hist(residuals)
    return np.std(residuals)

def calc_err_var_quadratic(data, last_n_pts=50, plot=False):
    a1 = fitting.Parameter(0.001)
    a2 = fitting.Parameter(0.005)
    b = fitting.Parameter(1.)
    def fit_func(t):
        return a1()*t**2 + a2()*t + b()

    t = np.arange(len(data))
    data_subset = data[-last_n_pts:]
    t_subset = t[-last_n_pts:]
    result = fitting.fit(fit_func, [a1, a2, b], data_subset, t_subset,
                         log_transform=False)
    residuals = result[0]
    if plot:
        plt.figure()
        plt.subplot(1, 2, 1)
        plt.plot(t, data)
        plt.plot(t_subset, fit_func(t_subset))
        plt.subplot(1, 2, 2)
        plt.hist(residuals)
    return np.std(residuals)

def calc_err_var_cubic(data, last_n_pts=50, plot=False):
    a1 = fitting.Parameter(0.0001)
    a2 = fitting.Parameter(0.001)
    a3 = fitting.Parameter(0.005)
    b = fitting.Parameter(1.)
    def fit_func(t):
        return a1()*t**3 + a2()*t**2 + a3()*t + b()

    t = np.arange(len(data))
    data_subset = data[-last_n_pts:]
    t_subset = t[-last_n_pts:]
    result = fitting.fit(fit_func, [a1, a2, a3, b], data_subset, t_subset,
                         log_transform=False)
    residuals = result[0]
    if plot:
        plt.figure()
        plt.subplot(1, 2, 1)
        plt.plot(t, data)
        plt.plot(t_subset, fit_func(t_subset))
        plt.subplot(1, 2, 2)
        plt.hist(residuals)
    return np.std(residuals)

if __name__ == '__main__':

    from preprocess_data import data_126, data_54

    plt.ion()
    plt.close('all')

    nbd_126_data = data_126[0, 0, :]
    nbd_126_std = calc_err_var_cubic(nbd_126_data, last_n_pts=80, plot=True)
    lbl = "126C NBD, sigma: %s" % nbd_126_std
    plt.title(lbl)
    print lbl

    nbd_54_data = data_54[0, 0, :]
    nbd_54_std = calc_err_var_cubic(nbd_54_data, last_n_pts=80, plot=True)
    lbl = "54C NBD, sigma: %s" % nbd_54_std
    plt.title(lbl)
    print lbl

    fret_126_data = data_126[0, 1, :]
    fret_126_std = calc_err_var_cubic(fret_126_data, last_n_pts=80, plot=True)
    lbl = "126C FRET, sigma: %s" % fret_126_std
    plt.title(lbl)
    print lbl

    fret_54_data = data_54[0, 1, :]
    fret_54_std = calc_err_var_cubic(fret_54_data, last_n_pts=80, plot=True)
    lbl = "54C FRET, sigma: %s" % fret_54_std
    plt.title(lbl)
    print lbl


