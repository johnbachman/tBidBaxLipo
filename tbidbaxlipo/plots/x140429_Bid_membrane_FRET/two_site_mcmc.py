"""
Fitting equation taken from:

Wang, Zhi-Xin, Jiang, Ruo-Fan, "A novel two-site binding equation
presented in terms of the total ligand concentration." FEBS Letters 392
(1996) 245-249.
"""

import numpy as np
from matplotlib import pyplot as plt
import emcee
import calculate_fret as cf
import triangle

(bid_concs, fret_means, fret_ses) = cf.get_fret_from_endpoints()
bid568_conc = 10.
Lstar = bid568_conc

nwalkers = 500
burn_steps = 1000
sample_steps = 200
ndim = 6

# Parameters are in order:
# log(K1), log(K2), log(R1), log(R2), N, F



def model_func(position, L):
    # Transform all params to linear scale
    log_params = 10 ** position[:4]
    lin_params = position[4:]
    (K1, K2, R1, R2) = log_params
    (N, F) = lin_params
    Lt = L + Lstar
    # From paper
    a = ((((K1 + K2) * (1 + N)) + R1 + R2 - Lt) / (1 + N))
    b = (((K1 * K2 * (1 + N)) + (K2 * R1) + (K1 * R2) - ((K1 + K2) * Lt)) /
         (1 + N))
    c = (-K1 * K2 * Lt) / (1 + N)
    theta = np.arccos((-2 * a**3 + 9*a*b - 27*c) / (2*np.sqrt((a**2 - 3*b)**3)))
    Lstar_b = Lstar - ((Lstar / (L + Lstar)) *
                          ((2/3.) * np.sqrt(a**2 - 3*b) * np.cos(theta/3.) - (a/3.)))
    frac_bound = Lstar_b / Lstar
    return F * frac_bound

def prior(position):
    (K1, K2, R1, R2, N, F) = position
    # Concentration parameters should be between 1pm and 1 uM
    if K1 < -4 or K1 > 4 or K2 < -4 or K2 > 4 or \
       R1 < -4 or R1 > 4 or R2 < -4 or R2 > 4:
        return -np.inf
    # F is a percentage, so should be between 0 and 1
    if F < 0 or F > 1:
        return -np.inf
    # N (nonspecific binding partition coefficient) should also be
    # between 0 and 1:
    if N < 0 or N > 1:
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
            p0[walk_ix, 0] = np.random.uniform(-4, 3) # log(K1)
            p0[walk_ix, 1] = np.random.uniform(-4, 4) # log(K2)
            p0[walk_ix, 2] = np.random.uniform(-4, 3) # log(R1)
            p0[walk_ix, 3] = np.random.uniform(-4, 5) # log(R2)
            p0[walk_ix, 4] = np.random.uniform(0, 1) # N
            p0[walk_ix, 5] = np.random.uniform(0, 1) # F

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
