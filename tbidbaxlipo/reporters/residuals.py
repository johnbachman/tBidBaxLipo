from bayessb.report import reporter, Result, ThumbnailResult
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from bayessb.multichain import NoPositionsException
import numpy as np
import math

reporter_group_name = "Residuals"

@reporter('Residuals at Maximum Likelihood')
def residuals_at_max_likelihood(mcmc_set):
    # Get the maximum likelihood parameters
    try:
        (max_likelihood, max_likelihood_position) = mcmc_set.maximum_likelihood()
    except NoPositionsException as npe:
        return Result(None, None)

    # Get the residuals
    residuals = mcmc_set.chains[0].get_residuals(max_likelihood_position)

    # Make the residuals plot
    fig = Figure()
    ax = fig.gca()
    plot_filename = '%s_max_likelihood_residuals.png' % mcmc_set.name
    thumbnail_filename = '%s_max_likelihood_residuals_th.png' % mcmc_set.name
    ax.plot(residuals[0], residuals[1])
    ax.set_title('Residuals at Maximum Likelihood')
    #ax.xlabel('Time')
    #ax.ylabel('Residual')
    canvas = FigureCanvasAgg(fig)
    fig.set_canvas(canvas)
    fig.savefig(plot_filename)
    fig.savefig(thumbnail_filename, dpi=10)

    return ThumbnailResult(thumbnail_filename, plot_filename)

@reporter('ACF of ML Residuals')
def acf_of_ml_residuals(mcmc_set):
    # Get the maximum likelihood parameters
    try:
        (max_likelihood, max_likelihood_position) = mcmc_set.maximum_likelihood()
    except NoPositionsException as npe:
        return Result(None, None)

    # Get the residuals
    residuals = mcmc_set.chains[0].get_residuals(max_likelihood_position)

    # Plot the autocorrelation function
    acf = np.correlate(residuals[1], residuals[1], mode='full')

    plot_filename = '%s_acf_of_ml_residuals.png' % mcmc_set.name
    thumbnail_filename = '%s_acf_of_ml_residuals_th.png' % mcmc_set.name
    fig = Figure()
    ax = fig.gca()
    ax.plot(acf)
    ax.set_title('Autocorrelation of Maximum Likelihood Residuals')

    canvas = FigureCanvasAgg(fig)
    fig.set_canvas(canvas)
    fig.savefig(plot_filename)
    fig.savefig(thumbnail_filename, dpi=10)

    return ThumbnailResult(thumbnail_filename, plot_filename)

@reporter('Durbin-Watson Statistic of ML Residuals')
def durbin_watson_of_ml_residuals(mcmc_set):
    # Get the maximum likelihood parameters
    try:
        (max_likelihood, max_likelihood_position) = mcmc_set.maximum_likelihood()
    except NoPositionsException as npe:
        return Result(None, None)

    # Get the residuals
    residuals = mcmc_set.chains[0].get_residuals(max_likelihood_position)
    e = residuals[1]

    # Calculate Durbin-Watson
    d = float(np.sum([math.pow(e[i+1] - e[i], 2) for i in range(len(e)-1)])) / \
        float(np.sum([math.pow(i, 2) for i in e]))
    return Result(d, None, expectation=2.0)


