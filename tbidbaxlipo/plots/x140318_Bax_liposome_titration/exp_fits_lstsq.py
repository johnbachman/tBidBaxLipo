"""
Fit the data from the liposome titration experiment from 140318 to a
simple two-parameter exponential equation by least squares minimization.
"""

import numpy as np
from matplotlib import pyplot as plt
from preprocess_data import data_to_fit, bg_time, lipo_concs_to_fit
from tbidbaxlipo.util import fitting, set_fig_params_for_publication, \
                             format_axis
from scipy.stats import linregress

def plot_exp_fits(time, data, concs, plot_fits=True, fmax_value=None):
    """Fit data to exponential, plot and return fitted values.

    Parameters
    ----------
    time : np.array
        Time vector for each timecourse.
    data : list of np.array
        List of arrays containing the fluorescence timecourses at each
        liposome concentration.
    concs : list
        List of concentrations of liposomes (not lipids) in nanomolar.
    plot_fits : boolean
        If True (default), the fits to the data are plotted, each in its own
        figure. If False, the fits are performed but the results are not
        plotted.

    Returns
    -------
    tuple : (fmax_arr, k_arr, sum_sq_err)
        fmax_arr and k_arr are numpy arrays containing the fitted values for
        the fmax and k parameters, respectively. sum_sq_err is the sum of the
        squares of the values of the residuals at the best fit parameters.
    """
    fmax_arr = np.zeros(len(concs))
    k_arr = np.zeros(len(concs))

    # Create the figure, if applicable
    if plot_fits:
        set_fig_params_for_publication()
        if fmax_value is None:
            figname = 'exp_fits_curves_fmax_var'
        else:
            figname = 'exp_fits_curves_fmax_fixed'
        plt.figure(figname, figsize=(1.5, 1.5), dpi=300)

    # Iterate over all of the liposome concentrations
    for i, conc in enumerate(concs):
        y = data[i]

        if fmax_value is None:
            fmax = fitting.Parameter(3.85)
        else:
            fmax = fitting.Parameter(fmax_value)

        k = fitting.Parameter(5e-4)
        def fit_func(t):
            return 1 + fmax() * (1 - np.exp(-k()*t))

        if fmax_value is None:
            fit_result = fitting.fit(fit_func, [k, fmax], y, time)
        else:
            fit_result = fitting.fit(fit_func, [k], y, time)

        fmax_arr[i] = fmax()
        k_arr[i] = k()

        if plot_fits:
            # Plot the data
            plt.plot(time, y, 'k')
            # Plot the fit to the data
            plt.plot(time, fit_func(time), 'r', zorder=3)

    # Format the finished plot
    if plot_fits:
        # Axis labels
        plt.ylabel('$F/F_0$')
        plt.xlabel(r'Time (sec $\times 10^3$)')
        # Axis bounds and ticks
        ax = plt.gca()
        ax.set_ylim([0.7, 5.2])
        ax.set_xlim([0, bg_time[-1] + 500])
        ax.set_xticks(np.linspace(0, 1e4, 6))
        ax.set_xticklabels([int(f) for f in np.linspace(0, 10, 6)])
        format_axis(ax)
        plt.subplots_adjust(bottom=0.24, left=0.21)
        plt.show()

    residuals = fit_result[0]
    sum_sq_err = np.sum(residuals ** 2)
    return (fmax_arr, k_arr, sum_sq_err)

def plot_k_fmax_varying(fmax_arr, k_arr, conc_arr):
    """Plot Fmax and k vs. liposome concentration for fits with Fmax
    allowed to vary.
    """
    set_fig_params_for_publication()
    plt.figure('exp_fits_fmax_var', figsize=(1.9, 1.5), dpi=300)
    # Plot Fmax vs. concentration on the left-hand axis
    plt.plot(conc_arr[:-1], fmax_arr[:-1], marker='o', markersize=3,
             color='b')
    plt.xlabel('[Liposomes] (nM)')
    ax1 = plt.gca()
    ax1.set_xscale('log')
    ax1.set_xlim([0.05, 30])
    ax1.set_ylabel('$F_{max}$', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    ax1.set_yscale('log')
    # Plot k vs. concentration on the right-hand axis
    ax2 = ax1.twinx()
    ax2.set_xlim([0.05, 30])
    ax2.plot(conc_arr[:-1], k_arr[:-1], marker='o', markersize=3, color='r')
    ax2.set_ylabel(r'k (sec $\times$ 10^{-3})', color='r')
    ax2.set_yscale('log')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    #ax2.set_yticks(np.linspace(6.6e-4, 7.8e-4, 7))
    #ax2.set_yticklabels(['%.1f' % f for f in np.linspace(6.6, 7.8, 7)])
    format_axis(ax1)
    format_axis(ax2, yticks_position='right')
    plt.subplots_adjust(left=0.18, bottom=0.19, right=0.75)

    plt.show()

def plot_k_fmax_fixed(k_arr, conc_arr):
    """Plot fits of k with Fmax fixed to a single value."""
    set_fig_params_for_publication()
    plt.figure('exp_fits_fmax_fixed', figsize=(1.5, 1.5), dpi=300)
    # Plot k vs. concentration on the left-hand axis
    plt.plot(conc_arr[:-1], k_arr[:-1], marker='o', markersize=3,
             color='r', linestyle='', zorder=3)
    plt.xlabel('[Liposomes] (nM)')
    ax1 = plt.gca()
    ax1.set_ylabel(r'k (sec$^{-1}$)', color='r')
    ax1.set_xscale('log')
    ax1.set_xlim([0.05, 30])
    ax1.set_yscale('log')
    ax1.set_ylim([1e-6, 1.2e-3])
    for tl in ax1.get_yticklabels():
        tl.set_color('r')
    format_axis(ax1)
    plt.subplots_adjust(left=0.30, bottom=0.19, right=0.93, top=0.93)

    # Fit to a line
    lin_fit = linregress(conc_arr[:-1], k_arr[:-1])
    plt.plot(conc_arr[:-1], lin_fit[0] * conc_arr[:-1] + lin_fit[1], color='b')
    print("---")
    print("Linear fit:")
    print(lin_fit)

    lslope = fitting.Parameter(1.0)
    lintercept = fitting.Parameter(np.exp(-9.6))
    def power_law(x):
        return lintercept() * x ** lslope()
    plaw_fit = fitting.fit(power_law, [lslope, lintercept], k_arr[:-1],
                           conc_arr[:-1])
    print("----")
    print("Power law fit (y = int * x ** slope):")
    print("intercept: %s" % lintercept())
    print("slope: %s" % lslope())
    plt.plot(conc_arr[:-1], power_law(conc_arr[:-1]), color='g')

    log_fit = linregress(np.log(conc_arr[:-1]), np.log(k_arr[:-1]))
    plt.plot(conc_arr[:-1],
             np.exp(log_fit[1]) * (conc_arr[:-1] ** log_fit[0]), color='k')
    print("----")
    print("Log-log linear fit:")
    print(log_fit)

    plt.show()

if __name__ == '__main__':
    plt.ion()

    # First run with fmax free to vary, show that this results in artifacts
    (fmax_arr, k_arr, sum_sq_err) = plot_exp_fits(bg_time, data_to_fit,
                                        lipo_concs_to_fit, plot_fits=True)
    plot_k_fmax_varying(fmax_arr, k_arr, lipo_concs_to_fit)
    plt.figure('exp_fits_fmax_var')
    plt.savefig('140318_exp_fits_lstsq_fmax_var.pdf')
    plt.figure('exp_fits_curves_fmax_var')
    plt.savefig('140318_exp_fits_lstsq_curves_fmax_var.pdf')
    print("Fmax variable, sum squared error: %s" % sum_sq_err)

    # Then re-run with Fmax-fixed
    (fmax_arr, k_arr, sum_sq_err) = plot_exp_fits(bg_time, data_to_fit,
                                        lipo_concs_to_fit, plot_fits=True,
                                        fmax_value=3.85)
    plot_k_fmax_fixed(k_arr, lipo_concs_to_fit)
    plt.figure('exp_fits_fmax_fixed')
    plt.savefig('140318_exp_fits_lstsq_fmax_fixed.pdf')
    plt.figure('exp_fits_curves_fmax_fixed')
    plt.savefig('140318_exp_fits_lstsq_curves_fmax_fixed.pdf')
    print("Fmax fixed, sum squared error: %s" % sum_sq_err)

