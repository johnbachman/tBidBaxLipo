"""
Fitting equation taken from:

Wang, Zhi-Xin, "An exact mathematical expression for describing competitive
binding of two different ligands to a protein molecule." FEBS Letters 360
(1995) 111-114.
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import emcee
import calculate_fret as cf
import triangle
import pickle
from os.path import exists
from tbidbaxlipo.util import set_fig_params_for_publication, format_axis, \
                             fontsize, capsize
import os.path
import sys

(bid_concs, fret_means, fret_ses) = cf.get_fret_from_endpoints()
bid568_conc = 10.
# A represents the labeled ligand, B the unlabeled ligand.
A0 = bid568_conc

nwalkers = 500
burn_steps = 100
ntemps = 10
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
    conc_A_spec = ((A0 * expr) / (3*K + expr))
    conc_A_nonspec = N * A0
    frac_A_bound = (conc_A_spec + conc_A_nonspec) / float(A0)
    return F * frac_A_bound

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
    # N (nonspecific FRET) should also be between 0 and 1:
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

def run_pt_mcmc(p0=None, filename=None):
    if p0 is None:
        p0 = np.zeros((ntemps, nwalkers, ndim))
        for temp_ix in range(ntemps):
            for walk_ix in range(nwalkers):
                p0[temp_ix, walk_ix, 0] = np.random.uniform(-4, 4) # log(K)
                p0[temp_ix, walk_ix, 1] = np.random.uniform(-0.5, 3) # log(P0)
                p0[temp_ix, walk_ix, 2] = np.random.uniform(0, 1) # N
                p0[temp_ix, walk_ix, 3] = np.random.uniform(0, 1) # F

    sampler = emcee.PTSampler(ntemps, nwalkers, ndim, likelihood, prior,
                              pool=None)

    print "Burn in sampling..."
    nstep = 0
    for p, lnprob, lnlike in sampler.sample(p0, iterations=burn_steps,
                            storechain=False):
        if nstep % 10 == 0:
            print "nstep %d of %d, MAP: %f" % (nstep, burn_steps,
                                               np.max(lnprob[0]))
        nstep +=1
    sampler.reset()
    print "Main sampling..."
    nstep = 0
    for p, lnprob, lnlike in sampler.sample(p, lnprob0=lnprob, lnlike0=lnlike,
                        iterations=sample_steps):
        if nstep % 10 == 0:
            print "nstep %d of %d, MAP: %f" % (nstep, sample_steps,
                                               np.max(lnprob[0]))
        nstep +=1

    if filename is not None:
        with open(filename, 'w') as f:
            pickle.dump(sampler, f)

    return sampler

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

def plot_chain(flatchain, lnprob):
    # Check convergence
    plt.figure('Chains')
    plt.plot(lnprob.T)
    plt.xlabel('MCMC Step')
    plt.ylabel('log(Posterior)')
    plt.title('Chain convergence')

    # Get sampling of trajectories parameters
    set_fig_params_for_publication()
    plt.figure('Fits', figsize=(1.5, 1.5), dpi=300)
    num_tot_steps = flatchain.shape[0]
    num_plot_samples = 2500
    n_pts = 100
    (x_lb, x_ub) = (0.1, 1500)
    ypred_samples = np.zeros((num_plot_samples, n_pts))
    bid_pred = np.logspace(np.log10(x_lb), np.log10(x_ub), n_pts)
    for s_ix in range(num_plot_samples):
        p_ix = np.random.randint(num_tot_steps)
        p_samp = flatchain[p_ix]
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
    ml_ix = np.argmax(lnprob)
    ml_pos = flatchain[ml_ix]
    plt.plot(bid_pred, model_func(ml_pos, bid_pred), color='r',
             linewidth=0.5)
    # Plot no competitor line, with stderr
    plt.hlines(fret_means[-1], x_lb, x_ub, linestyle='solid', color='k',
               linewidth=0.5)
    # Draw dashed line
    dash_line_pts = np.logspace(np.log10(x_lb), np.log10(x_ub), 50)
    for pt_ix in range(0, len(dash_line_pts) - 1, 2):
        plt.hlines([fret_means[-1] + fret_ses[-1], fret_means[-1] - fret_ses[-1]],
                   dash_line_pts[pt_ix], dash_line_pts[pt_ix + 1], color='k',
                   linewidth=0.5, zorder=3)
    # Plot data
    plt.errorbar(bid_concs[:-1], fret_means[:-1], yerr=fret_ses[:-1], color='k',
                 linewidth=1, capsize=capsize, zorder=3, linestyle='', marker='o',
                 markersize=2)
    # Label axes
    plt.xlabel('[cBid] (nM)')
    plt.ylabel(r'FRET (\%)')
    # Format axes
    ax.set_xscale('log')
    plt.xlim([x_lb, x_ub])
    ax.set_xticks([0.1, 1, 10, 100, 1000])
    ax.set_xticklabels([0.1, 1, 10, 100, 1000])
    ax.set_yticks([0.05, 0.15, 0.25, 0.35])
    ax.set_ylim([0.05, 0.37])
    format_axis(ax)
    plt.subplots_adjust(bottom=0.18, left=0.25, right=0.94, top=0.94)

    # Triangle plots
    #(tri_fig, axes) = plt.subplots(4, 4, figsize=(6, 6))
    triangle.corner(flatchain,
                    labels=['$log_{10}(K_D)$ (nM)', '$log_{10}(P)$ (nM)',
                            'Nonspec', 'FRET'])
    plt.subplots_adjust(right=0.96, top=0.96)

def plot_saturation_binding_predictions(flatchain):
    def binding_func(position, bid_tot, lipo_scale=1):
        # Transform all params to linear scale
        log_params = 10 ** position[:2]
        lin_params = position[2:]
        (kd, l_tot) = log_params
        (f0, fmax) = lin_params
        # 1.55 was concentration of lipos used in the experiment
        l_tot = l_tot / 1.55
        l_tot *= lipo_scale
        frac_bid_bound = ((1.0 / (2 * bid_tot)) *
                          (l_tot + bid_tot + kd -
                           np.sqrt((l_tot + bid_tot + kd)**2 -
                                   (4 * l_tot * bid_tot))))
        return frac_bid_bound

    plt.figure('Satbinding', figsize=(1.5, 1.5), dpi=300)
    num_tot_steps = flatchain.shape[0]
    num_plot_samples = 2500
    n_pts = 100
    (x_lb, x_ub) = (0.1, 1500)
    ypred_samples = np.zeros((num_plot_samples, n_pts))
    bid_pred = np.logspace(np.log10(x_lb), np.log10(x_ub), n_pts)
    for s_ix in range(num_plot_samples):
        p_ix = np.random.randint(num_tot_steps)
        p_samp = flatchain[p_ix]
        ypred = binding_func(p_samp, bid_pred, lipo_scale=1.55)
        ypred_samples[s_ix] = ypred
        #plt.plot(bid_pred, ypred, alpha=0.01, color='r')
    samp_lb = np.zeros(n_pts)
    samp_ub = np.zeros(n_pts)
    # Plot 95% confidence interval
    ax = plt.gca()
    # Label axes
    plt.xlabel('[cBid] (nM)')
    plt.ylabel(r'\% [cBid] bound')
    # Format axes
    ax.set_xscale('log')
    plt.xlim([x_lb, x_ub])
    ax.set_xticks([0.1, 1, 10, 100, 1000])
    ax.set_xticklabels([0.1, 1, 10, 100, 1000])
    plt.subplots_adjust(bottom=0.18, left=0.22, right=0.94, top=0.94)
    ypred_lb = np.percentile(ypred_samples, 2.5, axis=0)
    ypred_ub = np.percentile(ypred_samples, 97.5, axis=0)
    ax.fill_between(bid_pred, ypred_lb, ypred_ub, color='lightgray')
    format_axis(ax)
    # Plot mean of predictions
    plt.plot(bid_pred, np.mean(ypred_samples, axis=0), color='r')

    # Plot binding curve for different liposome concentrations
    plt.figure('Lipobinding', figsize=(1.5, 1.5), dpi=300)
    num_tot_steps = flatchain.shape[0]
    num_plot_samples = 2500
    n_pts = 100
    (x_lb, x_ub) = (0.01, 50)
    ypred_samples = np.zeros((num_plot_samples, n_pts))
    lipo_pred = np.logspace(np.log10(x_lb), np.log10(x_ub), n_pts)
    bid_pred = 20.
    for s_ix in range(num_plot_samples):
        p_ix = np.random.randint(num_tot_steps)
        p_samp = flatchain[p_ix]
        ypred = binding_func(p_samp, bid_pred, lipo_scale=lipo_pred)
        ypred_samples[s_ix] = ypred
        #plt.plot(lipo_pred, ypred, alpha=0.01, color='r')
    samp_lb = np.zeros(n_pts)
    samp_ub = np.zeros(n_pts)
    # Plot 95% confidence interval
    ax = plt.gca()
    # Label axes
    plt.xlabel('[Lipos] (nM)')
    plt.ylabel(r'\% cBid bound')
    # Format axes
    ax.set_xscale('log')
    plt.xlim([x_lb, x_ub])
    ax.set_xticks([0.01, 0.1, 1, 10,])
    ax.set_xticklabels([0.01, 0.1, 1, 10])
    plt.subplots_adjust(bottom=0.18, left=0.22, right=0.94, top=0.94)
    ypred_lb = np.percentile(ypred_samples, 2.5, axis=0)
    ypred_ub = np.percentile(ypred_samples, 97.5, axis=0)
    ax.fill_between(lipo_pred, ypred_lb, ypred_ub, color='lightgray')
    format_axis(ax)
    # Plot mean of predictions
    plt.plot(lipo_pred, np.mean(ypred_samples, axis=0), color='r')

if __name__ == '__main__':
    plt.ion()

    usage_msg =  "\nUsage:\n"
    usage_msg += "To run the fits and save chain to a pickle file:\n"
    usage_msg += "     python %s sample output_filename\n" % sys.argv[0]
    usage_msg += "To plot results from a pickled chain:\n"
    usage_msg += "     python %s plot input_filename\n" % sys.argv[0]

    # Check command-line args
    if len(sys.argv) < 3:
        print usage_msg
        sys.exit()

    # Sample
    pck_filename = sys.argv[2]
    if sys.argv[1] == 'sample':
        sampler = run_pt_mcmc(filename=pck_filename)
    # Plot
    elif sys.argv[1] == 'plot':
        with open(pck_filename) as f:
            sampler = pickle.load(f)
        plot_chain(sampler.flatchain[0], sampler.lnprobability[0])
        plot_saturation_binding_predictions(sampler.flatchain[0])
        plt.figure('Fits')
        plt.savefig('140429_exact_comp_bind_fit.pdf')
        plt.figure(3)
        plt.savefig('140429_exact_comp_bind_marginals.pdf')
    else:
        print usage_msg
        sys.exit()


