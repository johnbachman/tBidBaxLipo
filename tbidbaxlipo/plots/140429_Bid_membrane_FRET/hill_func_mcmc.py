import numpy as np
from matplotlib import pyplot as plt
import emcee
import calculate_fret as cf
import triangle

(bid_concs, fret_means, fret_ses) = cf.get_fret_from_endpoints()

nwalkers = 500
burn_steps = 100
sample_steps = 100
ndim = 4

# Parameters are in order: log10(kd), f, f0, sigma

def model_func(position):
    pass

def prior(position):
    (logkd, f, f0, sigma) = position
    # log10(kd) should be between -3 and 3, say
    if logkd < -3 or logkd > 3:
        return -np.inf
    # F and F0 are percentages, so they should be between 0 and 1
    if f < 0 or f > 1 or f0 < 0 or f0 > 1:
        return -np.inf
    # We don't expect the SD of the error to be greater than 20% FRET
    # when the range of the assay is between 0 and 40% FRET
    if sigma < 0 or sigma > 0.2:
        return -np.inf
    return 0

def likelihood(position):
    loglkl = 0
    return loglkl

def posterior(position):
    return prior(position) + likelihood(position)

def run_mcmc(p0=None):
    if p0 is None:
        p0 = np.zeros((nwalkers, ndim))
        for walk_ix in range(nwalkers):
            p0[walk_ix, 0] = np.random.uniform(-3, 3) # log10(kd)
            p0[walk_ix, 1] = np.random.uniform(0, 1) # F
            p0[walk_ix, 2] = np.random.uniform(0, 1) # F0
            p0[walk_ix, 3] = np.random.uniform(0, 0.20) # sigma

    sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior)
    pos, prob, state = sampler.run_mcmc(p0, burn_steps, storechain=True)

    return sampler

if __name__ == '__main__':
    sampler = run_mcmc()
    triangle.corner(sampler.flatchain)
