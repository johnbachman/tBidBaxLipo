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

def plot_emcee_fits_subplots(gf, sampler):
    plt.figure()
    for data_ix, data in enumerate(gf.data):
        plt.subplot(3, 4, data_ix + 1)
        plt.plot(gf.time, data)
        gf.plot_func_single(sampler.flatchain[0,-1,:], data_ix)

def plot_emcee_fits(gf, sampler, sample=True, burn=None, nsamples=100):
    """Plot fits from the MCMC chain vs. the data."""
    set_fig_params_for_publication()

    # If we're plotting samples, get the indices now and use them for
    # all observables
    if sample:
        (nwalkers, nsteps) = sampler.chain.shape[1:3]
        if burn is None:
            burn = int(nsteps / 2)

        walker_indices = np.random.randint(0, nwalkers, size=nsamples)
        step_indices = np.random.randint(burn, nsteps, size=nsamples)

    for obs_ix in range(gf.data.shape[1]):
        fig = plt.figure(figsize=(3, 3), dpi=300)
        plt.ylabel('$F/F_0$')
        plt.xlabel(r'Time (sec $\times 10^3$)')
        #plt.ylim([0.7, 5.2])
        plt.xlim([0, gf.time[-1] + 500])
        ax = plt.gca()
        #ax.set_xticks(np.linspace(0, 1e4, 6))
        #ax.set_xticklabels([int(f) for f in np.linspace(0, 10, 6)])
        plt.subplots_adjust(bottom=0.24, left=0.21)
        # Plot the different observables
        for cond_ix in range(gf.data.shape[0]):
            data = gf.data[cond_ix, obs_ix, :]
            plt.plot(gf.time, data, 'k', linewidth=1)
        # Colors for different observables
        obs_colors = ['r', 'g', 'b', 'k']
        # If we're plotting samples:
        if sample:
            for i in xrange(nsamples):
                p = sampler.chain[0, walker_indices[i], step_indices[i], :]
                plot_args = {'color': obs_colors[obs_ix], 'alpha': 0.1}
                gf.plot_func(p, obs_ix=obs_ix, plot_args=plot_args)

        # Plot the maximum a posteriori fit
        maxp_flat_ix = np.argmax(sampler.lnprobability[0])
        maxp_ix = np.unravel_index(maxp_flat_ix,
                                   sampler.lnprobability[0].shape)
        maxp = sampler.lnprobability[0][maxp_ix]
        plot_args = {'color': 'm', 'alpha': 1}
        gf.plot_func(sampler.chain[0, maxp_ix[0], maxp_ix[1]], obs_ix=obs_ix,
                     plot_args=plot_args)

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
    plot_emcee_fits(gf, sampler, burn=None, sample=True)
    #plot_emcee_fits_subplots(gf, sampler)

