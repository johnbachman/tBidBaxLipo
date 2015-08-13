"""A script that shows the triangle plots of the parameters and the posterior
probability at different steps of the chain for inspection of convergence.
"""
import sys
import pickle
import triangle
import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.util import set_fig_params_for_publication, format_axis
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import os

def save_fig(fig, plot_filename, display):
    if plot_filename and not display:
        canvas = FigureCanvasAgg(fig)
        fig.set_canvas(canvas)
        fig.savefig('%s.pdf' % plot_filename)
        fig.savefig('%s.png' % plot_filename)
    elif plot_filename:
        plt.savefig('%s.pdf' % plot_filename)
        plt.savefig('%s.png' % plot_filename)
    else:
        raise Exception("No plot_filename given, doing nothing.")

def triangle_plots(gf, sampler, plot_filename=None):
    """Triangle plots of lowest and highest temperature chains."""
    chain = sampler.flatchain
    num_params = chain.shape[2]
    # Lowest temp
    if DISPLAY:
        low_fig = plt.figure(figsize=(5, 5))
        low_fig.subplots(num_params, num_params)
        med_fig = plt.figure(figsize=(5, 5))
        med_fig.subplots(num_params, num_params)
        hi_fig = plt.figure(figsize=(5, 5))
        hi_fig.subplots(num_params, num_params)
    else:
        low_fig = Figure(figsize=(5, 5))
        med_fig = Figure(figsize=(5, 5))
        hi_fig = Figure(figsize=(5, 5))
        ix = 1
        for i in range(num_params):
            for j in range(num_params):
                low_fig.add_subplot(num_params, num_params, ix)
                med_fig.add_subplot(num_params, num_params, ix)
                hi_fig.add_subplot(num_params, num_params, ix)
                ix += 1

    triangle.corner(chain[0], fig=low_fig,
                    labels=[p.name for p in gf.builder.global_params])
    low_fig.suptitle('Triangle plot, lowest temp')
    # Intermediate temp
    temp_ix = 4
    triangle.corner(chain[temp_ix], fig=med_fig,
                    labels=[p.name for p in gf.builder.global_params])
    med_fig.suptitle('Triangle plot, temp %s' % temp_ix)
    # Highest temp
    triangle.corner(chain[-1], fig=hi_fig,
                    labels=[p.name for p in gf.builder.global_params])
    hi_fig.suptitle('Triangle plot, highest temp')

    if plot_filename:
        save_fig(low_fig, '%s.low' % plot_filename, DISPLAY)
        save_fig(med_fig, '%s.med' % plot_filename, DISPLAY)
        save_fig(hi_fig, '%s.hi' % plot_filename, DISPLAY)

def plot_chain_convergence(sampler, plot_filename):
    """Four-panel plot showing convergence of four lowest-temperature chains.

    Useful for determining if the temperature spacing at the lowest part of
    the chain is adequate.
    """
    if not DISPLAY and plot_filename is None:
        raise ValueError("DISPLAY is set to False but plot_filename is None, "
                         "so there will be no output.")

    ntemps = sampler.chain.shape[0]
    max_plots = 12
    interval = int(ntemps / (max_plots - 1))
    ncols = 4
    nrows = 3
    if DISPLAY:
        fig = plt.figure('Chain convergence', figsize=(16, 9))
    else:
        fig = Figure(figsize=(16, 9))
    temp_indices = np.round(np.linspace(0, ntemps - 1, max_plots))
    for fig_ix, temp_ix in enumerate(temp_indices):
        ax = fig.add_subplot(nrows, ncols, fig_ix + 1)
        ax.plot(sampler._lnprob[temp_ix,:,:].T, alpha=0.1)
        ax.set_title('Chain %d' % temp_ix)
    if plot_filename:
        save_fig(fig, plot_filename, DISPLAY)
    #plt.tight_layout()
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

def plot_emcee_fits(gf, sampler, sample=True, burn=None, nsamples=100,
                    plot_filename=None):
    """Plot fits from the MCMC chain vs. the data."""
    if not DISPLAY and plot_filename is None:
        raise ValueError("DISPLAY is set to False but plot_filename is None, "
                         "so there will be no output.")
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
        if DISPLAY:
            fig = plt.figure(figsize=(3, 3), dpi=300)
        else:
            fig = Figure(figsize=(3, 3), dpi=300)

        ax = fig.gca()
        ax.set_ylabel('$F/F_0$')
        #ax.set_xlabel(r'Time (sec $\times 10^3$)')
        ax.set_xlabel(r'Time')
        #plt.ylim([0.7, 5.2])
        ax.set_xlim([0, gf.time[-1] + 500])
        #ax.set_xticks(np.linspace(0, 1e4, 6))
        #ax.set_xticklabels([int(f) for f in np.linspace(0, 10, 6)])
        fig.subplots_adjust(bottom=0.21, left=0.20)
        # Plot the different observables
        for cond_ix in range(gf.data.shape[0]):
            data = gf.data[cond_ix, obs_ix, :]
            ax.plot(gf.time, data, 'k', linewidth=1)
        # Colors for different observables
        obs_colors = ['r', 'g', 'b', 'k']
        # If we're plotting samples:
        if sample:
            for i in xrange(nsamples):
                p = sampler.chain[0, walker_indices[i], step_indices[i], :]
                plot_args = {'color': obs_colors[obs_ix], 'alpha': 0.1}
                gf.plot_func(p, obs_ix=obs_ix, ax=ax, plot_args=plot_args)

        # Plot the maximum a posteriori fit
        maxp_flat_ix = np.argmax(sampler.lnprobability[0])
        maxp_ix = np.unravel_index(maxp_flat_ix,
                                   sampler.lnprobability[0].shape)
        maxp = sampler.lnprobability[0][maxp_ix]
        plot_args = {'color': 'm', 'alpha': 1}
        gf.plot_func(sampler.chain[0, maxp_ix[0], maxp_ix[1]], ax=ax,
                     obs_ix=obs_ix, plot_args=plot_args)

        format_axis(ax)
        if plot_filename:
            obs_plot_filename = '%s.obs%d' % (plot_filename, obs_ix)
            save_fig(fig, obs_plot_filename, DISPLAY)

def plot_conformations(gf, sampler, sample=True, burn=None, nsamples=100,
                       plot_filename=None):
    """Plot fluorescence conformations from the MCMC chain vs. the data."""
    if not DISPLAY and plot_filename is None:
        raise ValueError("DISPLAY is set to False but plot_filename is None, "
                         "so there will be no output.")
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
        if DISPLAY:
            fig = plt.figure(figsize=(3, 3), dpi=300)
        else:
            fig = Figure(figsize=(3, 3), dpi=300)
        ax = fig.gca()
        ax.set_ylabel('$F/F_0$')
        #ax.set_xlabel(r'Time (sec $\times 10^3$)')
        ax.set_xlabel(r'Time')
        #plt.ylim([0.7, 5.2])
        ax.set_xlim([0, gf.time[-1] + 500])
        #ax.set_xticks(np.linspace(0, 1e4, 6))
        #ax.set_xticklabels([int(f) for f in np.linspace(0, 10, 6)])
        fig.subplots_adjust(bottom=0.24, left=0.21)
        # Plot the different observables
        for cond_ix in range(gf.data.shape[0]):
            data = gf.data[cond_ix, obs_ix, :]
            ax.plot(gf.time, data, 'k', linewidth=1)
        # Colors for different observables
        obs_colors = ['r', 'g', 'b', 'm', 'k']

        def plot_confs(x, plot_args):
            for cond_ix in range(gf.data.shape[0]):
                gf.set_parameters(x, obs_ix=obs_ix, cond_ix=cond_ix)
                gf.solver.run()
                # Iterate over the number of conformations
                for conf_ix in range(gf.builder.num_confs):
                    conf_y = gf.solver.yobs['Bax_c%d' % conf_ix]
                    conf_scaling = \
                       gf.builder.model.parameters['c%d_scaling' %
                                                     conf_ix].value
                    y = conf_y * conf_scaling
                    if 'color' not in plot_args:
                        ax.plot(gf.solver.tspan, y, label='c%d' % conf_ix,
                                 color=obs_colors[conf_ix], **plot_args)
                    else:
                        ax.plot(gf.solver.tspan, y, label='c%d' % conf_ix,
                                 **plot_args)

        # If we're plotting samples:
        if sample:
            for i in xrange(nsamples):
                p = sampler.chain[0, walker_indices[i], step_indices[i], :]
                plot_args = {'alpha': 0.1}
                plot_confs(p, plot_args=plot_args)

        # Plot the maximum a posteriori fit
        maxp_flat_ix = np.argmax(sampler.lnprobability[0])
        maxp_ix = np.unravel_index(maxp_flat_ix,
                                   sampler.lnprobability[0].shape)
        maxp = sampler.lnprobability[0][maxp_ix]
        plot_args = {'color': 'k', 'alpha': 1, 'linewidth':0.5}
        plot_confs(sampler.chain[0, maxp_ix[0], maxp_ix[1]],
                   plot_args=plot_args)

        format_axis(ax)
        if plot_filename:
            save_fig(fig, plot_filename, DISPLAY)

if __name__ == '__main__':

    usage_msg =  "Usage:\n"
    usage_msg += " python show_chain.py chain_filename.mcmc output_dir\n"

    if len(sys.argv) < 3:
        print usage_msg
        sys.exit()

    chain_filename = sys.argv[1]
    output_dir = sys.argv[2]

    # Unpickle the chain
    print("Loading %s" % chain_filename)
    with open(chain_filename) as f:
        (gf, sampler) = pickle.load(f)

    DISPLAY = False

    # Show plots
    #plt.ion()
    #print("Plotting triangle plots")
    output_base = os.path.join(output_dir, os.path.basename(chain_filename))
    #triangle_plots(gf, sampler, plot_filename=chain_filename + '.tri')
    print("Plotting convergence")
    plot_chain_convergence(sampler, output_base + '.conv')
    print("Plotting sample fits")
    plot_emcee_fits(gf, sampler, burn=None, sample=True,
                    plot_filename=output_base + '.fits')
    print("Plotting conformation timecourses")
    #plot_conformations(gf, sampler, burn=None, sample=True,
    #                   plot_filename=output_base + '.confs')
    #plot_emcee_fits_subplots(gf, sampler)

