import sys
import pickle
import numpy as np
from matplotlib import pyplot as plt
from preprocess_data import data_to_fit, bg_time
from tbidbaxlipo.util import set_fig_params_for_publication, format_axis
import corner

def plot_emcee_fits(gf, sampler):
    """Plot fits from the MCMC chain vs. the data."""
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
    for data_ix in range(data_to_fit.shape[0]):
        plt.plot(bg_time, data_to_fit[data_ix, 0, :], 'k', linewidth=1)
    # Plot the final point (should probably plot max likelihood instead)
    gf.plot_func(sampler.flatchain[0,-1,:])
    format_axis(ax)

if __name__ == '__main__':
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
        corner.corner(sampler.flatchain[0])
        plt.savefig('%s_tri.pdf' % output_base)
        sys.exit()

