from bayessb.report import reporter, Result, ThumbnailResult
import numpy as np
from matplotlib import pyplot as plt

reporter_group_name = 'Experiments'
num_samples = 100

@reporter('tBid 10x')
def increase_tBid_10x(mcmc_set):
    """ .. todo:: document the basis for this"""
    tspan = mcmc_set.chains[0].options.tspan
    plt.figure()
    plot_filename = '%s_tBid_10x.png' % mcmc_set.name
    thumbnail_filename = '%s_tBid_10x_th.png' % mcmc_set.name

    # Make sure we can call the method 'get_observable_timecourse'
    if not hasattr(mcmc_set.chains[0], 'get_observable_timecourse') or \
       not hasattr(mcmc_set.chains[0], 'plot_data'):
        return Result('None', None)

    # Plot the original data
    mcmc_set.chains[0].plot_data()

    # Plot a sampling of trajectories from the original parameter set
    for i in range(num_samples):
        position = mcmc_set.get_sample_position()
        x = mcmc_set.chains[0].get_observable_timecourse(position=position)
        plt.plot(tspan, x, color='g', alpha=0.5)

    # Now increase tBid 10-fold...
    model = mcmc_set.chains[0].options.model
    old_tBid_0 = model.parameters['tBid_0'].value
    model.parameters['tBid_0'].value = old_tBid_0 * 10

    # ...and do the sampling again
    for i in range(num_samples):
        position = mcmc_set.get_sample_position()
        x = mcmc_set.chains[0].get_observable_timecourse(position=position)
        plt.plot(tspan, x, color='r', alpha=0.5)

    plt.savefig(plot_filename)
    plt.savefig(thumbnail_filename, dpi=10)

    # Make sure to reset the tBid initial condition to its original value!
    model.parameters['tBid_0'].value = old_tBid_0

    return ThumbnailResult(thumbnail_filename, plot_filename)

@reporter('Bax 5x')
def increase_Bax_5x(mcmc_set):
    """ .. todo:: document the basis for this"""
    tspan = mcmc_set.chains[0].options.tspan
    plt.figure()
    plot_filename = '%s_Bax_5x.png' % mcmc_set.name
    thumbnail_filename = '%s_Bax_5x_th.png' % mcmc_set.name

    # Make sure we can call the method 'get_observable_timecourse'
    if not hasattr(mcmc_set.chains[0], 'get_observable_timecourse') or \
       not hasattr(mcmc_set.chains[0], 'plot_data'):
        return Result('None', None)

    # Plot the original data
    mcmc_set.chains[0].plot_data()

    # Plot a sampling of trajectories from the original parameter set
    for i in range(num_samples):
        position = mcmc_set.get_sample_position()
        x = mcmc_set.chains[0].get_observable_timecourse(position=position)
        plt.plot(tspan, x, color='g', alpha=0.5)

    # Now increase Bax 5-fold...
    model = mcmc_set.chains[0].options.model
    old_Bax_0 = model.parameters['Bax_0'].value
    model.parameters['Bax_0'].value = old_Bax_0 * 5

    # ...and do the sampling again
    for i in range(num_samples):
        position = mcmc_set.get_sample_position()
        x = mcmc_set.chains[0].get_observable_timecourse(position=position)
        plt.plot(tspan, x, color='r', alpha=0.5)

    plt.savefig(plot_filename)
    plt.savefig(thumbnail_filename, dpi=10)

    # Make sure to reset the Bax initial condition to its original value!
    model.parameters['Bax_0'].value = old_Bax_0

    return ThumbnailResult(thumbnail_filename, plot_filename)

