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
import corner
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
burn_steps = 200
ntemps = 10
sample_steps = 200
ndim = 5

# Model taken from:
# Schwarz, G., & Beschiaschvili, G. (1989). Thermodynamic and kinetic studies
# on the association of melittin with a phospholipid bilayer. Biochimica Et
# Biophysica Acta, 979(1), 82-90.

def calc_b(temp=37):
    # See equation A5c
    # e0 = "the elementary charge"
    # Value from Wikipedia: https://en.wikipedia.org/wiki/Elementary_charge
    e0 = 1.6021766208e-19 # coulombs
    # A_L = The area per lipid on the outer surface
    # A_L = 0.75 # nm^2 / lipid
    A_L = 0.75 * 10**-18 # M^2 / lipid
    # Beta = the fraction of lipid making up the outer leaflet
    beta = 0.5
    # eps_v = the permittivity of vacuum
    # Value from Wikipedia: https://en.wikipedia.org/wiki/Vacuum_permittivity
    eps_v = 8.854187817e-12 # (farads per meter)
    # eps = the dielectric constant of the aqueous phase
    # Wikipedia: 80.1 at 20C
    eps = 80.1
    R = 8.3144598 # J / (K*mol)
    # Temperature in Kelvin (37C = 310K)
    T = 273. + temp
    # Concentration of salt
    cs = 0.135 # Molar (135 mM KCl)

    # Units of numerator:
    # e0 / (A_L * B)
    # coulombs * m^-2

    # Units of denominator (inside square root):
    # 8 ev*e * RT * cs
    #  = farads m^-1 J K^-1 mol^-1 K mol
    #  = farads * m^-1 * J
    # In SI base units:
    # (coul^2 s^-2 s^4 kg^-1 m^-2) * m^-1 * (kg m^2 s^-2)
    # (all seconds cancel, kg cancels)
    # coul^2 m^-1
    # After square root:
    # sqrt(coul^2 m^-1)
    # coul m^-1/2
    #
    # numerator / denominator:
    # coul * m^-2 * coul-1 m^1/2
    # m^-1.5 ???

# Parameters are in order:
# v, b, pc, fret, f0
def binding_func(position, bid_concs, liposome_conc=1.55):
    lipid_conc_per_nm = 83.787 # 129.87 uM / 1.55 nM
    lipid_conc = lipid_conc_per_nm * liposome_conc

    # v: charges per peptide chain
    # b: dimensionless
    # pc: partition coefficient
    # f: FRET efficiency
    # f0: baseline FRET
    (v, b, pc, fret, f0) = position
    # r: cb/L, ratio of bound peptide to lipid
    r = np.logspace(-6, 0, 1000)

    # cf: concentration of free peptide
    # Equation (6), Schwarz et al. 1989
    def cf(r):
        return (r * alpha(r)) / pc
    # alpha: activity coefficient, depends on r
    # Equation (7) in Schwarz et al. 1989
    def alpha(r):
        log_alpha = 2 * v * np.arcsinh(v * b * r)
        return np.exp(log_alpha)

    # Bound peptide is (cb/L) * L
    cb = r * lipid_conc
    # cf, calculated from r
    cfr = cf(r)
    # ctot: total peptide, calculated from the r values
    ctot = cb + cfr
    fracs_bound = cb / ctot
    # Interpolate to find frac bound for 
    interp_fracs = np.interp(bid_concs, ctot, fracs_bound)
    return interp_fracs

def fret_func(position, bid_concs, liposome_conc=1.55):
    (v, b, pc, fret, f0) = position
    interp_fracs = binding_func(position, bid_concs, liposome_conc)
    return (fret * interp_fracs) + f0

def prior(position):
    (v, b, pc, fret, f0) = position
    # v: charges per peptide chain
    if v < 0 or v > 10:
        return -np.inf
    # b: dimensionless
    if b < 0 or b > 100:
        return -np.inf
    # pc: partition cofficient is a percentage, so should be between 0 and 1
    if pc < 0 or pc > 1:
        return -np.inf
    # fret: between 0 and 1:
    if fret < 0 or fret > 1:
        return -np.inf
    # F0 (nonspecific FRET) should also be between 0 and 1:
    if f0 < 0 or f0 > 1.0:
        return -np.inf
    return 0

def likelihood(position):
    # Get the model prediction for this parameter set
    ypred = fret_func(position, bid_concs)
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
                p0[temp_ix, walk_ix, 0] = np.random.uniform(0, 10) # v
                p0[temp_ix, walk_ix, 1] = np.random.uniform(0, 50) # b
                p0[temp_ix, walk_ix, 2] = np.random.uniform(0, 1) # pc
                p0[temp_ix, walk_ix, 3] = np.random.uniform(0, 1) # fret
                p0[temp_ix, walk_ix, 4] = np.random.uniform(0, 1) # f0

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
        ypred = fret_func(p_samp, bid_pred)
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
    plt.plot(bid_pred, fret_func(ml_pos, bid_pred), color='r',
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
    # Triangle plots
    #(tri_fig, tri_axes) = plt.subplots(5, 5, figsize=(5, 5), dpi=150)
    (tri_fig, tri_axes) = plt.subplots(2, 2, figsize=(2, 2), dpi=300)
    corner.corner(flatchain[:, 0:2], fig=tri_fig,
                  #labels=['V', 'B', 'PC', 'FRET', 'F0'],
                  labels=['V (cBid eff. charge)', 'B'],
                  plot_datapoints=False, no_fill_contours=True)
    for row_ix, row in enumerate(tri_axes):
        for col_ix, ax in enumerate(row):
            #format_axis(ax, label_padding=20)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_tick_params(labelsize=fontsize, pad=1, length=1.5,
                                     width=0.5)
            ax.yaxis.set_tick_params(labelsize=fontsize, pad=1, length=1.5,
                                     width=0.5)
            ax.xaxis.label.set_size(fontsize)
            ax.yaxis.label.set_size(fontsize)
            ax.xaxis.set_label_coords(0.5, -0.2)
            ax.yaxis.set_label_coords(-0.2, 0.5)
            if col_ix == 0:
                ax.axvline(2, color='r', alpha=0.5)
            if col_ix == 0 and row_ix == 1:
                ax.axhline(11.5, color='r', alpha=0.5)
            if row_ix == 1:
                ax.axvline(11.5, color='r', alpha=0.5)

    tri_fig.subplots_adjust(right=0.96, top=0.96, bottom=0.15, left=0.15,
                            hspace=0.15, wspace=0.15)

def plot_saturation_binding_predictions(flatchain, plot_filename=None):
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
        ypred = np.array([binding_func(p_samp, bid_pred, liposome_conc=l)
                          for l in lipo_pred])
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

    # Plot previously published value
    plt.plot(0.9, 0.5, marker='o', markersize=4, color='b')

    if plot_filename:
        plt.savefig('%s.pdf' % plot_filename)
        plt.savefig('%s.png' % plot_filename, dpi=300)

if __name__ == '__main__':
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
        plot_saturation_binding_predictions(sampler.flatchain[0],
                                            '140429_gouy_chap_bind_pred')
        print "MAP:"
        print np.max(sampler.lnprobability[0])
        print "Max lkl:"
        print np.max(sampler.lnlikelihood[0])
        plt.figure('Fits')
        plt.savefig('140429_gouy_chap_fit.pdf')
        plt.figure(3)
        plt.savefig('140429_gouy_chap_marginals.pdf')
    else:
        print usage_msg
        sys.exit()



