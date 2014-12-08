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
import pickle
from os.path import exists
from tbidbaxlipo.util import set_fig_params_for_publication, format_axis

(bid_concs, fret_means, fret_ses) = cf.get_fret_from_endpoints()
bid568_conc = 10.
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

def run_mcmc(p0=None, filename=None):
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

    if filename is not None:
        with open(filename, 'w') as f:
            pickle.dump(sampler, f)

    return sampler

def plot_chain(sampler):
    # Check convergence
    """
    plt.figure('Chains')
    plt.plot(sampler.lnprobability.T)
    plt.xlabel('MCMC Step')
    plt.ylabel('log(Posterior)')
    plt.title('Chain convergence')
    """

    # Get sampling of trajectories parameters
    set_fig_params_for_publication()
    plt.figure('Fits', figsize=(1.5, 1.5), dpi=300)
    num_tot_steps = sampler.flatchain.shape[0]
    num_plot_samples = 2500
    n_pts = 100
    (x_lb, x_ub) = (0.1, 1500)
    ypred_samples = np.zeros((num_plot_samples, n_pts))
    bid_pred = np.logspace(np.log10(x_lb), np.log10(x_ub), n_pts)
    for s_ix in range(num_plot_samples):
        p_ix = np.random.randint(num_tot_steps)
        p_samp = sampler.flatchain[p_ix]
        ypred = model_func(p_samp, bid_pred)
        ypred_samples[s_ix] = ypred
        #plt.plot(bid_pred, ypred, alpha=0.01, color='r')
    samp_lb = np.zeros(n_pts)
    samp_ub = np.zeros(n_pts)
    # Plot 95% confidence interval
    ax = plt.gca()
    ypred_lb = np.percentile(ypred_samples, 2.5, axis=0)
    ypred_ub = np.percentile(ypred_samples, 97.5, axis=0)
    ax.fill_between(bid_pred, ypred_lb, ypred_ub, color='lightgray')
    # Plot maximum likelihood
    ml_ix = np.unravel_index(np.argmax(sampler.lnprobability),
                             sampler.lnprobability.shape)
    ml_pos = sampler.chain[ml_ix]
    plt.plot(bid_pred, model_func(ml_pos, bid_pred), color='r',
             linewidth=0.5)
    # Plot no competitor line, with stderr
    plt.hlines(fret_means[-1], x_lb, x_ub, linestyle='solid', color='k',
               linewidth=0.5)
    plt.hlines(fret_means[-1] + fret_ses[-1], x_lb, x_ub,
               linestyle='dashed', color='k', linewidth=0.5)
    plt.hlines(fret_means[-1] - fret_ses[-1], x_lb, x_ub,
               linestyle='dashed', color='k', linewidth=0.5)
    # Plot data
    plt.errorbar(bid_concs[:-1], fret_means[:-1], yerr=fret_ses[:-1], color='k',
                 linewidth=1, capsize=1.5, zorder=3)
    # Label axes
    plt.xlabel('[cBid] (nM)')
    plt.ylabel(r'FRET (\%)')
    # Format axes
    format_axis(ax)
    ax.set_xscale('log')
    plt.xlim([x_lb, x_ub])
    #ax.set_xticks(np.logspace(0.1, 1000, 5))
    ax.set_yticks([0.05, 0.15, 0.25, 0.35])
    ax.set_ylim([0.05, 0.37])
    #ax.set_xticklabels([int(f) for f in np.linspace(-1, 3, 5)])
    plt.subplots_adjust(bottom=0.18, left=0.25, right=0.94, top=0.94)

    # Triangle plots
    #triangle.corner(sampler.flatchain)
    import ipdb; ipdb.set_trace()

if __name__ == '__main__':
    plt.ion()

    pck_filename = 'exact_comp_bind_mcmc.pck'
    if exists(pck_filename):
        with open(pck_filename) as f:
            sampler = pickle.load(f)
    else:
        sampler = run_mcmc(filename=pck_filename)

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
