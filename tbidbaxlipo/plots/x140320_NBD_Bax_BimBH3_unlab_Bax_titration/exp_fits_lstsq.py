from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import fitting, format_axis, \
                             set_fig_params_for_publication
from tbidbaxlipo.util.plate_assay import plot_all, TIME, VALUE
from preprocess_data import timecourse_wells, nbd_wells, bax_bg_wells, \
                            timecourse_averages, bgsub_wells, \
                            bax_concs_to_fit, time, data_to_fit
from scipy import stats

def plot_data():
    """Plots the data and various transformations of it."""

    # Timecourse wells
    plt.figure()
    plot_all(timecourse_wells)
    plt.title("Raw timecourses")

    # NBD wells
    plt.figure()
    plot_all(nbd_wells)
    plt.title("NBD-Bax wells, Raw")

    # Lipo bg wells
    plt.figure()
    plot_all(bax_bg_wells)
    plt.title("Unlabeled Bax background wells, Raw")

    # Averages
    plt.figure()
    plot_all(timecourse_averages)
    plt.title("Average timecourses, all")

    # Background-subtracted wells
    plt.figure()
    plot_all(bgsub_wells)
    plt.title("Background-subtracted wells")

def plot_exp_fits(time, data_to_fit, bax_concs_to_fit):
    fmax_list = []
    k1_list = []

    # Create figure showing the fits to each timecourse
    set_fig_params_for_publication()
    plt.figure('exp_fit_curves', figsize=(1.5, 1.5), dpi=300)

    for i, bax_conc in enumerate(bax_concs_to_fit):
        # Try fitting the high conc trajectory
        y = data_to_fit[i, 0, :]
        y = y / np.mean(y[0:2])
        fmax = fitting.Parameter(25.)
        k1 = fitting.Parameter(np.log(2)/2000.)
        k2 = fitting.Parameter(1e-3)
        b = fitting.Parameter(np.mean(y[0:2]))
        k_bleach = fitting.Parameter(1.17e-5)
        #def exp_func(t):
        #    return (b() + fmax()*(1 - np.exp(-k1() *
        #            (1 - np.exp(-k2()*t))*t))) * np.exp(-k_bleach()*t)
        # One-parameter exp
        def exp_func(t):
            return (b() + fmax()*(1 - np.exp(-k1()*t))) * \
                    np.exp(-k_bleach()*t)

        #fitting.fit(exp_func, [fmax, k1, k2], y, t)
        fitting.fit(exp_func, [fmax, k1], y, time)

        # Add to list
        fmax_list.append(fmax())
        k1_list.append(k1())

        # Plot the data
        plt.plot(time, y, color='k', linewidth=1)
        # Plot fits to the data
        plt.plot(time, exp_func(time), color='r', zorder=3, linewidth=1)

    # Format the finished plot
    # Axis labels
    plt.xlabel(r'Time (sec $\times 10^3$)')
    plt.ylabel('$F/F_0$')
    # Axis bounds and ticks
    ax = plt.gca()
    ax.set_ylim([0.7, 4])
    ax.set_xlim([0, time[-1] + 500])
    ax.set_xticks(np.linspace(0, 1e4, 6))
    ax.set_xticklabels([int(f) for f in np.linspace(0, 10, 6)])
    format_axis(ax)
    plt.subplots_adjust(bottom=0.21, left=0.23)

    # Now, plot the scaling of the parameters with concentration
    plt.figure('exp_fits', figsize=(1.7, 1.5), dpi=300)
    log_concs = np.log10(np.array(bax_concs_to_fit[:-1]))
    slope = fitting.Parameter(1.)
    intercept = fitting.Parameter(1.)
    def log_line(x):
        log_concs = np.log10(x)
        return slope() * log_concs + intercept()
    fitting.fit(log_line, [slope, intercept], fmax_list[:-1],
                bax_concs_to_fit[:-1])
    plt.plot(bax_concs_to_fit[:-1], fmax_list[:-1], marker='o', markersize=3,
             color='b', linestyle='')
    plt.plot(bax_concs_to_fit[:-1], log_line(bax_concs_to_fit[:-1]), color='b')

    plt.xlabel('[Bax] (nM)')
    ax1 = plt.gca()
    ax1.set_xscale('log')
    ax1.set_xlim([8, 1000])
    ax1.set_ylabel('$F_{max}$', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')

    ax2 = ax1.twinx()
    ax2.set_xlim([5, 1000])
    ax2.plot(bax_concs_to_fit[:-1], k1_list[:-1], marker='o', markersize=3,
             color='r', linestyle='')
    slope = fitting.Parameter(5e-4)
    intercept = fitting.Parameter(6.6e-3)
    fitting.fit(log_line, [slope, intercept], k1_list[:-1],
                bax_concs_to_fit[:-1])
    ax2.plot(bax_concs_to_fit[:-1], log_line(bax_concs_to_fit[:-1]),
                                    color='r')

    ax2.set_ylabel(r'k (sec $\times\ 10^{-3}$)', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    ax2.set_yticks(np.linspace(6.6e-4, 7.8e-4, 7))
    ax2.set_yticklabels(['%.1f' % f for f in np.linspace(6.6, 7.8, 7)])

    format_axis(ax1)
    format_axis(ax2, yticks_position='right')
    plt.subplots_adjust(left=0.20, bottom=0.19, right=0.80)

if __name__ == '__main__':
    plot_exp_fits(time, data_to_fit, bax_concs_to_fit)

    plt.figure('exp_fit_curves')
    plt.savefig('140320_exp_fit_curves.pdf')

    plt.figure('exp_fits')
    plt.savefig('140320_exp_fits.pdf')
