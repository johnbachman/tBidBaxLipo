from bayessb.report import reporter, Result, ThumbnailResult
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

reporter_group_name = 'Experiments'
num_samples = 100

@reporter('tBid 0.1x')
def increase_tBid_01x(mcmc_set):
    """ .. todo:: document the basis for this"""
    return change_tBid_concentration(mcmc_set, 0.1)

@reporter('tBid 0.5x')
def increase_tBid_05x(mcmc_set):
    """ .. todo:: document the basis for this"""
    return change_tBid_concentration(mcmc_set, 0.5)

@reporter('tBid 10x')
def increase_tBid_10x(mcmc_set):
    """ .. todo:: document the basis for this"""
    return change_tBid_concentration(mcmc_set, 10)

def change_tBid_concentration(mcmc_set, fold_change):
    tspan = mcmc_set.chains[0].options.tspan
    fig = Figure()
    ax = fig.gca()
    plot_filename = '%s_tBid_%.1fx.png' % (mcmc_set.name, fold_change)
    thumbnail_filename = '%s_tBid_%.1fx_th.png' % (mcmc_set.name, fold_change)

    # Make sure we can call the method 'get_observable_timecourses'
    if not hasattr(mcmc_set.chains[0], 'get_observable_timecourses') or \
       not hasattr(mcmc_set.chains[0], 'plot_data'):
        return Result('None', None)

    # Plot the original data
    mcmc_set.chains[0].plot_data(ax)

    # Plot a sampling of trajectories from the original parameter set
    for i in range(num_samples):
        position = mcmc_set.get_sample_position()
        timecourses = mcmc_set.chains[0].get_observable_timecourses(position=position)
        for obs_name, timecourse in timecourses.iteritems():
            ax.plot(tspan, timecourse, color='g', alpha=0.5, label=obs_name)

    # Now change tBid x-fold...
    model = mcmc_set.chains[0].options.model
    old_tBid_0 = model.parameters['tBid_0'].value
    model.parameters['tBid_0'].value = old_tBid_0 * fold_change

    # ...and do the sampling again
    for i in range(num_samples):
        position = mcmc_set.get_sample_position()
        timecourses = mcmc_set.chains[0].get_observable_timecourses(position=position)
        for obs_name, timecourse in timecourses.iteritems():
            ax.plot(tspan, timecourse, color='r', alpha=0.5, label=obs_name)

    canvas = FigureCanvasAgg(fig)
    fig.set_canvas(canvas)
    fig.savefig(plot_filename)
    fig.savefig(thumbnail_filename, dpi=10)

    # Make sure to reset the tBid initial condition to its original value!
    model.parameters['tBid_0'].value = old_tBid_0

    return ThumbnailResult(thumbnail_filename, plot_filename)

@reporter('Bax 5x')
def increase_Bax_5x(mcmc_set):
    """ .. todo:: document the basis for this"""
    tspan = mcmc_set.chains[0].options.tspan
    fig = Figure()
    ax = fig.gca()
    plot_filename = '%s_Bax_5x.png' % mcmc_set.name
    thumbnail_filename = '%s_Bax_5x_th.png' % mcmc_set.name

    # Make sure we can call the method 'get_observable_timecourses'
    if not hasattr(mcmc_set.chains[0], 'get_observable_timecourses') or \
       not hasattr(mcmc_set.chains[0], 'plot_data'):
        return Result('None', None)

    # Plot the original data
    mcmc_set.chains[0].plot_data(ax)

    # Plot a sampling of trajectories from the original parameter set
    for i in range(num_samples):
        position = mcmc_set.get_sample_position()
        timecourses = mcmc_set.chains[0].get_observable_timecourses(position=position)
        for obs_name, timecourse in timecourses.iteritems():
            ax.plot(tspan, timecourse, color='g', alpha=0.5, label=obs_name)

    # Now increase Bax 5-fold...
    model = mcmc_set.chains[0].options.model
    old_Bax_0 = model.parameters['Bax_0'].value
    model.parameters['Bax_0'].value = old_Bax_0 * 5

    # ...and do the sampling again
    for i in range(num_samples):
        position = mcmc_set.get_sample_position()
        timecourses = mcmc_set.chains[0].get_observable_timecourses(position=position)
        for obs_name, timecourse in timecourses.iteritems():
            ax.plot(tspan, timecourse, color='r', alpha=0.5, label=obs_name)

    canvas = FigureCanvasAgg(fig)
    fig.set_canvas(canvas)
    fig.savefig(plot_filename)
    fig.savefig(thumbnail_filename, dpi=10)

    # Make sure to reset the Bax initial condition to its original value!
    model.parameters['Bax_0'].value = old_Bax_0

    return ThumbnailResult(thumbnail_filename, plot_filename)

