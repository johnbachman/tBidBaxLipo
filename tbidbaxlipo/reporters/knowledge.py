from bayessb.report import reporter, Result, MeanSdResult, FuzzyBooleanResult
from bayessb.report.evidence import Evidence, Citation
import numpy as np
from matplotlib import pyplot as plt

reporter_group_name = 'Prior knowledge'
num_samples = 100

e = Evidence("""
When DAC-tBid (20 nM) and NBD-126-Bax (100 nM) were incubated together, energy
transfer (detected as a decrease in donor fluorescence) was observed only if
membranes (~2.5 nM liposomes) were also present (Figure 2A), suggesting that
the two proteins interact only when bound to membranes.

The extent of energy transfer between DAC-tBid and NBD-126-Bax did not decrease
over time (Figure 2A), suggesting that at steady state some fraction of tBid
remained bound to Bax in the lipid membrane. This result was unexpected given
that tBid and Bax do not cofractionate when liposomes or mitochondria are
solubilized with the detergent CHAPS (Billen et al., 2008).
""",
image='http://ars.els-cdn.com/content/image/1-s2.0-S0092867408014396-gr2.jpg',
citation=Citation(
"""Lovell, J. F., Billen, L. P., Bindner, S., Shamas-Din, A., Fradin, C.,
Leber, B., & Andrews, D. W. (2008). Membrane binding by tBid initiates an
ordered series of events culminating in membrane permeabilization by Bax. Cell,
135(6), 1074-1084.""",
pmid='19062087',
doi='doi:10.1016/j.cell.2008.11.010'))

@reporter('tBid/Bax increases monotonically', evidence=e)
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
    return FuzzyBooleanResult(float(num_true) / float(num_samples),
                              plot_filename, expectation=1.0)

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
    return FuzzyBooleanResult(float(num_true) / float(num_samples),
                              plot_filename, expectation=1.0)

@reporter('tBid/Bax K_D')
def tBidBax_kd(mcmc_set):
    """ .. todo:: document the basis for this"""
    num_kds = 10000
    # Get indices for the tBid/Bax binding constants
    estimate_params = mcmc = mcmc_set.chains[0].options.estimate_params
    tBid_iBax_kf_index = None
    tBid_iBax_kr_index = None
    for i, p in enumerate(estimate_params):
        if p.name == 'tBid_iBax_kf':
            tBid_iBax_kf_index = i
        elif p.name == 'tBid_iBax_kr':
            tBid_iBax_kr_index = i
    # If we couldn't find the parameters, return None for the result
    if tBid_iBax_kf_index is None or tBid_iBax_kr_index is None:
        return Result(None, None)
    # Sample the kr/kf ratio across the pooled chains
    kd_dist = np.zeros(num_kds)
    for i in range(num_kds):
        position = mcmc_set.get_sample_position()
        kd_dist[i] = ((10 ** position[tBid_iBax_kr_index]) /
                      (10 ** position[tBid_iBax_kf_index]))
    # Calculate the mean and variance
    mean = kd_dist.mean()
    sd = kd_dist.std()
    # Plot the Kd distribution
    plot_filename = '%s_tBidiBax_kd_dist.png' % mcmc_set.name
    plt.figure()
    plt.hist(kd_dist)
    plt.savefig(plot_filename)
    return MeanSdResult(mean, sd, plot_filename)

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


