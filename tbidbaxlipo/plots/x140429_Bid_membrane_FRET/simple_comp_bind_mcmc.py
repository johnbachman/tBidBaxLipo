import numpy as np
from matplotlib import pyplot as plt
import emcee
import calculate_fret as cf
import corner
from scipy.stats import distributions as dist

(bid_concs, fret_means, fret_ses) = cf.get_fret_from_endpoints()
bid568_conc = 10.

nwalkers = 500
burn_steps = 200
sample_steps = 100
ndim = 5

# Parameters are in order: log10(ic50), f, nonspec, total, sigma

def model_func(position, conc):
    (log_ic50, f, nonspec, total, sigma) = position
    frac_bound = nonspec + ((total - nonspec) /
                            (1 + 10 ** (np.log10(conc) - log_ic50)))
    return f * frac_bound

def prior(position):
    (log_ic50, f, nonspec, total, sigma) = position
    # log10(ic50) should be between -3 and 3, say
    if log_ic50 < -3 or log_ic50 > 3:
        return -np.inf
    # F, nonspec, and total are percentages, so they should be between 0 and 1
    if f < 0 or f > 1 or nonspec < 0 or nonspec > 1 or total < 0 or total > 1:
        return -np.inf
    # We don't expect the SD of the error to be greater than 20% FRET
    # when the range of the assay is between 0 and 40% FRET
    if sigma < 0 or sigma > 0.2:
        return -np.inf
    return 0

def likelihood(position):
    sigma = position[-1]
    # Get the model prediction for this parameter set
    ypred = model_func(position, bid_concs)
    errs = ypred - fret_means
    loglkl = np.sum(dist.norm.logpdf(errs, loc=0, scale=sigma))
    if np.isnan(loglkl):
        return -np.inf
    else:
        return loglkl

def posterior(position):
    return prior(position) + likelihood(position)

def run_mcmc(p0=None):
    if p0 is None:
        p0 = np.zeros((nwalkers, ndim))
        for walk_ix in range(nwalkers):
            p0[walk_ix, 0] = np.random.uniform(-3, 3) # log10(ic50)
            p0[walk_ix, 1] = np.random.uniform(0, 1) # f
            p0[walk_ix, 2] = np.random.uniform(0, 1) # nonspec
            p0[walk_ix, 3] = np.random.uniform(0, 1) # total
            p0[walk_ix, 4] = np.random.uniform(0, 0.20) # sigma

    sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior)

    # Burn-in
    print("Burn-in sampling...")
    pos, prob, state = sampler.run_mcmc(p0, burn_steps, storechain=False)
    sampler.reset()
    # Main sampling
    print("Main sampling...")
    sampler.run_mcmc(pos, sample_steps)

    return sampler

def plot_chain(sampler):
    # Check convergence
    plt.figure()
    plt.plot(sampler.lnprobability.T)

    # Plot maximum likelihood
    ml_ix = np.unravel_index(np.argmax(sampler.lnprobability),
                             sampler.lnprobability.shape)
    ml_pos = sampler.chain[ml_ix]
    bid_pred = np.logspace(-1, 3, 100)
    plt.figure()
    plt.errorbar(np.log10(bid_concs), fret_means, yerr=fret_ses, color='k')
    plt.plot(np.log10(bid_pred), model_func(ml_pos, bid_pred), color='r')

    # Plot sampling of trajectories parameters
    plt.figure()
    num_plot_samples = 2500
    num_tot_steps = sampler.flatchain.shape[0]
    for s_ix in range(num_plot_samples):
        p_ix = np.random.randint(num_tot_steps)
        p_samp = sampler.flatchain[p_ix]
        plt.plot(np.log10(bid_pred), model_func(p_samp, bid_pred),
                 alpha=0.01, color='r')
    plt.errorbar(np.log10(bid_concs), fret_means, yerr=fret_ses, color='k',
                 linewidth=2)

    # Triangle plots
    corner.corner(sampler.flatchain)



if __name__ == '__main__':
    plt.ion()
    sampler = run_mcmc()
    plot_chain(sampler)
