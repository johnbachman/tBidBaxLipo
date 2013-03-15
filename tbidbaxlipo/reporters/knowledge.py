from bayessb.report import reporter, Result
import numpy as np
from matplotlib import pyplot as plt

reporter_group_name = 'Prior knowledge'
num_samples = 100

@reporter('tBid/Bax increases monotonically')
def tBid_Bax_monotonically_increasing(mcmc_set):
    """ .. todo:: document the basis for this"""
    tspan = mcmc_set.chains[0].options.tspan
    plt.figure()
    plot_filename = '%s_tBidBax_monotonic_increasing.png' % mcmc_set.name

    num_true = 0
    for i in range(num_samples): 
        x  = mcmc_set.get_sample_simulation()
        plt.plot(tspan, x['tBidBax'], color='r', alpha=0.5)
        if monotonic_increasing(x['tBidBax']):
            num_true += 1
    plt.savefig(plot_filename)
    return Result(float(num_true) / float(num_samples), plot_filename)

@reporter('iBax increases monotonically')
def iBax_monotonically_increasing(mcmc_set):
    """ .. todo:: document the basis for this"""
    tspan = mcmc_set.chains[0].options.tspan
    plt.figure()
    plot_filename = '%s_iBax_monotonic_increasing.png' % mcmc_set.name

    num_true = 0
    for i in range(num_samples): 
        x  = mcmc_set.get_sample_simulation()
        plt.plot(tspan, x['iBax'], color='r', alpha=0.5)
        if monotonic_increasing(x['iBax']):
            num_true += 1
    plt.savefig(plot_filename)
    return Result(float(num_true) / float(num_samples), plot_filename)

# Helper functions
# ================
def monotonic_increasing(x):
    # TODO rewrite so doesn't allow fixed, unchanging values
    dx = np.diff(x)
    return bool(np.all(dx >= 0))

def monotonic_decreasing(x):
    # TODO rewrite so doesn't allow fixed, unchanging values
    dx = np.diff(x)
    return bool(np.all(dx <= 0))


