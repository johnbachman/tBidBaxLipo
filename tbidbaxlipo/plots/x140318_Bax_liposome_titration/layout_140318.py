import sys
import pickle
import numpy as np
from matplotlib import pyplot as plt
from preprocess_data import timecourse_wells, lipo_bg_wells, bgsub_wells, \
                            layout, lipo_bg_layout, bax_lipo_layout, \
                            data_to_fit, bg_time, lipo_concs_to_fit
from tbidbaxlipo.util.plate_assay import plot_all, TIME, VALUE
from tbidbaxlipo.util import fitting, set_fig_params_for_publication, \
                             emcee_fit, format_axis
from tbidbaxlipo.models.nbd import multiconf
from tbidbaxlipo.models import one_cpt

def plot_emcee_fits(gf, sampler):
    set_fig_params_for_publication()
    fig = plt.figure(figsize=(1.5, 1.5), dpi=300)
    plt.ylabel('$F/F_0$')
    plt.xlabel(r'Time (sec $\times 10^3$)')
    plt.ylim([0.7, 5.2])
    plt.xlim([0, bg_time[-1] + 500])
    ax = plt.gca()
    ax.set_xticks(np.linspace(0, 1e4, 6))
    ax.set_xticklabels([int(f) for f in np.linspace(0, 10, 6)])
    plt.subplots_adjust(bottom=0.24, left=0.21)
    for data in data_to_fit:
        plt.plot(bg_time, data, 'k', linewidth=1)
    # Plot the final point
    gf.plot_func(sampler.flatchain[0,-1,:])
    format_axis(ax)

if __name__ == '__main__':
    import triangle
    plt.ion()

    if len(sys.argv) == 3:
        mcmc_path = sys.argv[1]
        output_base = sys.argv[2]
        # Load the sampler data
        (gf, sampler) = pickle.load(open(mcmc_path))
        # Plot the best fit vs. the data
        plot_emcee_fits(gf, sampler)
        plt.savefig('%s_fits.pdf' % output_base)
        # Triangle plot for parameters
        plt.figure()
        triangle.corner(sampler.flatchain[0])
        plt.savefig('%s_tri.pdf' % output_base)
        sys.exit()

