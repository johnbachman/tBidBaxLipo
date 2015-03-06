from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import fitting
import scipy.stats
import tbidbaxlipo.plots.titration_fits as tf

def calc_err_var(data, last_n_pts=50, fit_type='cubic', plot=False):
    # Prune the data
    t = np.arange(len(data))
    data_subset = data[-last_n_pts:]
    t_subset = t[-last_n_pts:]
    # Get the appropriate fit object
    if fit_type == 'linear':
        fit = tf.LinearIntercept(log_transform=False)
    elif fit_type == 'quadratic':
        fit = tf.Quadratic(log_transform=False)
    elif fit_type == 'cubic':
        fit = tf.Cubic(log_transform=False)
    else:
        raise ValueError("Unknown fit type! Must be 'linear', 'quadratic', or "
                         "'cubic'.")
    # Run the fit
    params = fit.fit_timecourse(t_subset, data_subset)
    print params
    # Get the best-fit curve
    ypred = fit.fit_func(t_subset, params)
    # Calculate residuals
    residuals = data_subset - ypred
    # Plot results
    fig = None
    if plot:
        # Plot of fit and distribution of residuals
        fig = plt.figure(figsize=(12, 5))
        plt.subplot(1, 3, 1)
        plt.plot(t, data)
        plt.plot(t_subset, ypred)
        plt.title('Fit of %s to last %s pts' % (fit_type, last_n_pts))
        plt.subplot(1, 3, 2)
        plt.hist(residuals)
        plt.title('Histogram of residuals')
        plt.subplot(1, 3, 3)
        scipy.stats.probplot(residuals, dist='norm', plot=plt)
        plt.title('Quantile-quantile plot vs. normal')
        plt.tight_layout()
    return (np.std(residuals), fig)

