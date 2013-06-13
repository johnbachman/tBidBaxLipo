"""
A script to fit Bid binding data from the Octet Red to a mechanistic model
of Bid binding.
"""

import numpy as np
from tbidbaxlipo.util import fitting
from matplotlib import pyplot as plt
from pysb import *
from pysb.integrate import Solver, odesolve
import bayessb

def likelihood(mcmc, position):
    """The likelihood function. Incorporates scaling parameters for the
    signal magnitudes of the different Bid conformational states."""
    yout = mcmc.simulate(position, observables=True)
    params = mcmc.cur_params(position)
    obs = (yout['Bid_m'] * params[0]) + (yout['Bid_i'] * params[1])
    # Assume a SD of 0.05
    return np.sum((y - obs)**2 / (2 * 0.05**2))

def step(mcmc):
    if mcmc.iter % 20 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  ' \
              'glob_acc=%-.3f  lkl=%g  prior=%g  post=%g' % \
              (mcmc.iter, mcmc.sig_value, mcmc.T,
               mcmc.acceptance/(mcmc.iter+1.), mcmc.accept_likelihood,
               mcmc.accept_prior, mcmc.accept_posterior)

if __name__ == '__main__':

    # Define a simple model -------------------------
    Model()
    Monomer('Bid', ['loc'], {'loc': ['c', 'm', 'i']})

    Parameter('scaling_Bid_m', 0.17)
    Parameter('scaling_Bid_i', 0.40)

    Initial(Bid(loc='c'), Parameter('Bid_0', 1))

    Rule('Bid_binds_liposomes', Bid(loc='c') >> Bid(loc='m'),
            Parameter('Bid_c_to_m_k', 0.215))
    #Rule('Bid_conf_change', Bid(loc='m') >> Bid(loc='i'),
    #        Parameter('Bid_m_to_i_k', 0.02))
    Rule('Bid_dimerize',
         Bid(loc='m') + Bid(loc='m') <> Bid(loc='i'),
         Parameter('Bid_m_to_i_k', 0.02), Parameter('Bid_i_to_m_k', 0.02))

    Observable('Bid_m', Bid(loc='m'))
    Observable('Bid_i', Bid(loc='i'))
    Observable('Bid_c', Bid(loc='c'))
    # ------------------------------------------------

    # Load the data
    data = np.loadtxt('cBid1uM.txt')
    time = data[0:2500:5,0]
    y = data[0:2500:5,1]

    # Initialize the MCMC
    opts = bayessb.MCMCOpts()
    opts.model = model
    opts.tspan = time
    opts.estimate_params = [p for p in model.parameters
                            if not p.name.endswith('_0')]

    opts.initial_values = [p.value for p in opts.estimate_params]
    opts.nsteps = 4000
    opts.likelihood_fn = likelihood
    opts.T_init = 1
    opts.anneal_length = 0
    opts.step_fn = step
    opts.use_hessian = True
    opts.hessian_scale = 1
    opts.hessian_period = opts.nsteps / 10

    opts.sigma_step = 0
    opts.norm_step_size = 0.001
    opts.seed = 1
    mcmc = bayessb.MCMC(opts)

    # Plot before
    x = odesolve(model, time)
    x_obs = (scaling_Bid_m.value * x['Bid_m']) + \
            (scaling_Bid_i.value * x['Bid_i'])
    plt.ion()
    plt.figure()
    plt.plot(time, y)
    plt.plot(time, x_obs)

    # Do estimation
    mcmc.initialize()
    mcmc.estimate()

    # Plot after
    best_fit_position = mcmc.positions[np.argmin(mcmc.posteriors)]
    best_fit_params = mcmc.cur_params(position=best_fit_position)
    x = mcmc.simulate(position=best_fit_position, observables=True)
    obs = (x['Bid_m'] * best_fit_params[0]) + \
          (x['Bid_i'] * best_fit_params[1])
    plt.plot(time, obs)

    plt.figure()
    plt.plot(time, x['Bid_c'])
    plt.plot(time, x['Bid_m'])
    plt.plot(time, x['Bid_i'])
    plt.plot(time, y)

