"""A script that shows the triangle plots of the parameters and the posterior
probability at different steps of the chain for inspection of convergence.
"""
import sys
import pickle
import triangle
import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.util import set_fig_params_for_publication, format_axis

def triangle_plots(gf, sampler):
    """Triangle plots of lowest and highest temperature chains."""
    chain = sampler.flatchain
    fig = triangle.corner(chain[0],
                          labels=[p.name for p in gf.builder.global_params])
    fig.suptitle('Triangle plot, lowest temp')
    #fig.savefig('triangle_low.png')
    fig = triangle.corner(chain[-1],
                          labels=[p.name for p in gf.builder.global_params])
    fig.suptitle('Triangle plot, highest temp')
    #fig.savefig('triangle_high.png')

def plot_chain_convergence(sampler):
    """Four-panel plot showing convergence of four lowest-temperature chains.

    Useful for determining if the temperature spacing at the lowest part of
    the chain is adequate.
    """
    plt.figure('Chain convergence')
    plt.subplot(2, 2, 1)
    plt.plot(sampler._lnprob[0,:,:].T, alpha=0.1)
    plt.title('0th chain')
    plt.subplot(2, 2, 2)
    plt.plot(sampler._lnprob[1,:,:].T, alpha=0.1)
    plt.title('1st chain')
    plt.subplot(2, 2, 3)
    plt.plot(sampler._lnprob[2,:,:].T, alpha=0.1)
    plt.title('2nd chain')
    #plt.subplot(2, 2, 4)
    #plt.plot(sampler._lnprob[3,:,:].T, alpha=0.1)
    #plt.title('3rd chain')
    #plt.savefig('chain_convergence.png')

def plot_emcee_fits(gf, sampler):
    """Plot fits from the MCMC chain vs. the data."""
    set_fig_params_for_publication()
    fig = plt.figure(figsize=(3, 3), dpi=300)
    plt.ylabel('$F/F_0$')
    plt.xlabel(r'Time (sec $\times 10^3$)')
    plt.ylim([0.7, 5.2])
    plt.xlim([0, gf.time[-1] + 500])
    ax = plt.gca()
    ax.set_xticks(np.linspace(0, 1e4, 6))
    ax.set_xticklabels([int(f) for f in np.linspace(0, 10, 6)])
    plt.subplots_adjust(bottom=0.24, left=0.21)
    for data in gf.data:
        plt.plot(gf.time, data, 'k', linewidth=1)
    # Plot the final point (should probably plot max likelihood instead)
    gf.plot_func(sampler.flatchain[0,-1,:])
    format_axis(ax)

if __name__ == '__main__':

    usage_msg =  "Usage:\n"
    usage_msg += " python show_chain.py chain_filename.mcmc\n"

    if len(sys.argv) < 2:
        print usage_msg
        sys.exit()

    chain_filename = sys.argv[1]

    # Unpickle the chain
    with open(chain_filename) as f:
        (gf, sampler) = pickle.load(f)

    # Show plots
    plt.ion()
    triangle_plots(gf, sampler)
    plot_chain_convergence(sampler)
    plot_emcee_fits(gf, sampler)

