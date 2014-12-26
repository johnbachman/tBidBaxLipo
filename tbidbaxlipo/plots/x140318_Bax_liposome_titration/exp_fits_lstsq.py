import numpy as np
from matplotlib import pyplot as plt

from preprocess_data import data_to_fit, bg_time, lipo_concs_to_fit
from tbidbaxlipo.util import fitting, set_fig_params_for_publication, \
                             format_axis

def plot_exp_fits(time, data, concs, plot_fits=True):
    fmax_arr = np.zeros(len(concs))
    k_arr = np.zeros(len(concs))
    for i, conc in enumerate(concs):
        y = data[i]
        fmax = fitting.Parameter(3.85)
        k = fitting.Parameter(5e-4)
        def fit_func(t):
            return 1 + fmax() * (1 - np.exp(-k()*t))
        fitting.fit(fit_func, [k, fmax], y, time)

        fmax_arr[i] = fmax()
        k_arr[i] = k()

        if plot_fits:
            plt.figure()
            plt.plot(time, y, 'b')
            plt.plot(time, fit_func(time), 'r')
            plt.title('%f nM liposomes' % conc)
            plt.xticks([0, 2000, 4000, 6000, 8000])
            plt.ylabel('$F/F_0$')
            plt.xlabel('Time (sec)')
            plt.show()

    return (fmax_arr, k_arr)

def plot_fmax_k_curves(fmax_arr, k_arr, conc_arr):

    set_fig_params_for_publication()

    plt.figure('exp_fits', figsize=(1.7, 1.5), dpi=300)
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

    ax2 = ax1.twinx()
    ax2.set_xlim([0.05, 30])
    ax2.plot(conc_arr[:-1], k_arr[:-1], marker='o', markersize=3, color='r')
    ax2.set_ylabel(r'k (sec $\times$ 10^{-3})', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    #ax2.set_yticks(np.linspace(6.6e-4, 7.8e-4, 7))
    #ax2.set_yticklabels(['%.1f' % f for f in np.linspace(6.6, 7.8, 7)])

    format_axis(ax1)
    format_axis(ax2, yticks_position='right')
    plt.subplots_adjust(left=0.20, bottom=0.19, right=0.80)

    plt.show()
    ax2.set_yscale('log')

if __name__ == '__main__':
    plt.ion()

    (fmax_arr, k_arr) = plot_exp_fits(bg_time, data_to_fit, lipo_concs_to_fit,
                                      plot_fits=False)
    plot_fmax_k_curves(fmax_arr, k_arr, lipo_concs_to_fit)

