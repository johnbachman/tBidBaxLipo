"""
Fitting equation taken from:

Wang, Zhi-Xin, "An exact mathematical expression for describing competitive
binding of two different ligands to a protein molecule." FEBS Letters 360
(1995) 111-114.
"""

import numpy as np
from matplotlib import pyplot as plt
import emcee
import calculate_fret as cf
import triangle

(bid_concs, fret_means, fret_ses) = cf.get_fret_from_endpoints()
bid568_conc = 10.
bid_concs = bid_concs - bid568_conc
# A represents the labeled ligand, B the unlabeled ligand.
A0 = bid568_conc

nwalkers = 500
burn_steps = 1500
sample_steps = 200
ndim = 4

# Parameters are in order:
# log(K), P0, N, F

def model_func(position, B0):
    # Transform all params to linear scale
    log_params = 10 ** position[:2]
    lin_params = position[2:]
    (K, P0) = log_params
    (N, F) = lin_params
    # From paper
    a = (K + K + A0 + B0 - P0)
    b = (K * (A0 - P0)) + (K * (B0 - P0)) + (K * K)
    c = -K * K * P0
    theta = np.arccos((-2 * a**3 + 9*a*b - 27*c) / (2*np.sqrt((a**2 - 3*b)**3)))
    expr = (2 * np.sqrt(a**2 - 3*b) * np.cos(theta / 3.) - a)
    frac_A_bound = expr / (3*K + expr)
    return F * frac_A_bound + N

def prior(position):
    (K, P0, N, F) = position
    # Concentration parameters should be between 1pm and 1 uM
    if K < -4 or K > 4:
        return -np.inf
    # There are roughly 1 nM liposomes, so we'll assume that there must be at
    # least 1 nM liposome binding sites (lower bound of log10(P0) = -0.5)
    if P0 < -0.5 or P0 > 3:
        return -np.inf
    # F is a percentage, so should be between 0 and 1
    if F < 0 or F > 1:
        return -np.inf
    # N (nonspecific FRET) should also be between 0 and 1 between 0 and 1:
    # However, after running with bounds between 0 and 1, there was a local
    # of values around 0.2 (with terrible fits), but with no values greater
    # than 0.2. Restricting to a max of 0.165 prevents getting stuck in these
    # bad fits.
    if N < 0 or N > 1.0:
        return -np.inf
    # We don't expect the SD of the error to be greater than 20% FRET
    # when the range of the assay is between 0 and 40% FRET
    #if sigma < 0 or sigma > 0.2:
    #    return -np.inf
    return 0

def likelihood(position):
    # Get the model prediction for this parameter set
    ypred = model_func(position, bid_concs)
    loglkl = -np.sum(((ypred - fret_means) ** 2) / (2 * fret_ses ** 2))
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
            p0[walk_ix, 0] = np.random.uniform(-4, 4) # log(K)
            p0[walk_ix, 1] = np.random.uniform(-0.5, 3) # log(P0)
            p0[walk_ix, 2] = np.random.uniform(0, 1) # N (see prior func)
            p0[walk_ix, 3] = np.random.uniform(0, 1) # F

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
    bid_pred = np.logspace(-1, 4, 100)
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
    triangle.corner(sampler.flatchain)

if __name__ == '__main__':
    plt.ion()
    sampler = run_mcmc()
    plot_chain(sampler)
    """
    import sys
    sys.exit()
    plt.figure()
    plt.errorbar(np.log10(bid_concs), fret_means, yerr=fret_ses, color='k',
                 linewidth=2)
    position = np.array([np.log10(1.0), np.log10(1.0), np.log10(1.5),
                         np.log10(1.5), 0.10, 1.])
    plt.plot(np.log10(bid_concs), model_func(position, bid_concs), marker='o')
    insulin = np.logspace(-2, 5, 100)
    position = np.array([np.log10(1.0), np.log10(1.0), np.log10(1.5),
                         np.log10(1.5), 0.10, 1.])
    ypred = model_func(position, insulin)
    plt.figure()
    plt.errorbar(np.log10(bid_concs), fret_means, yerr=fret_ses, color='k',
                 linewidth=2)
    plt.plot(np.log10(insulin), ypred, marker='o')
    posterior(position)
    """
