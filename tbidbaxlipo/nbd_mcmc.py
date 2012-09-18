"""Functions for fitting the tBid/Bax/Lipo models to the NBD-Bax mutant
fluorescence data using MCMC.
"""

import mcmc_hessian
from pysb.integrate import odesolve
import numpy
import matplotlib.pyplot as plt
from tBid_Bax_1c import tBid_Bax_1c
import nbd_analysis
import pickle

# Prepare the data
# ================

tspan = nbd_analysis.time_c62
nbd_avgs, nbd_stds = nbd_analysis.calc_norm_avg_std()
ydata_norm = nbd_avgs[1] # NBD_62c 

# Some useful parameter sets
# ==========================

# Some parameters that looked reasonable after manual fitting
params_dict_before = {
    'tBid_transloc_kf':1e-1,
    'tBid_transloc_kr':0,
    'Bax_transloc_kf':1e-1,
    'Bax_transloc_kr':100,
    'tBid_mBax_kf':1,
    'tBid_mBax_kr':0.01,
    'tBid_iBax_kc':10,
    'iBax_reverse_k':1e-1,
    'Bax_dimerization_kf':1e-3,
    'Bax_dimerization_kr':1e-2,
    'Bax_tetramerization_kf':1e-3,
    'Bax_tetramerization_kr':1e-4
}

# Some nice parameters from fitting starting from the above
params_dict_after = {'Bax_dimerization_kf': 0.00039997477906346336,
    'Bax_dimerization_kr': 0.013749131382056693,
    'Bax_tetramerization_kf': 0.00093066359281583318,
    'Bax_tetramerization_kr': 0.0002291780988903229,
    'Bax_transloc_kf': 5.007883980982208,
    'Bax_transloc_kr': 171.1669217159442,
    'iBax_reverse_k': 0.2713605496448912,
    'tBid_iBax_kc': 11.621431266275815,
    'tBid_mBax_kf': 0.062823108082142434,
    'tBid_mBax_kr': 0.036115738840543629,
    'tBid_transloc_kf': 0.046433226280901206,
    'tBid_transloc_kr': 0.0
}

# Build the model
# ===============

m1c = tBid_Bax_1c(params_dict=params_dict_before)
m1c.build_model2()
model = m1c.model

# MCMC Functions
# ==============

def do_fit():
    """Runs MCMC on the globally defined model."""

    # Define the likelihood function
    def likelihood(mcmc, position):
        """The likelihood function. Calculates the error between model and data
        to give a measure of the likelihood of observing the data given the
        current parameter values.
        """
        yout = mcmc.simulate(position)

        # TODO need a way to get the desired observable(s) from the array
        obs = 2*yout[:,6]+4*yout[:,7]
        obs_max = max(obs)
        yout_norm = obs / obs_max    

        return numpy.sum((ydata_norm - yout_norm) ** 2 / (2 * sigmas ** 2))

    # Set the random number generator seed
    seed = 2
    random = numpy.random.RandomState(seed)

    # Set the standard deviation
    sigmas = nbd_stds[1]

    # Initialize the MCMC arguments
    opts = mcmc_hessian.MCMCOpts()
    opts.model = model
    opts.tspan = tspan

    print model.parameters
    # estimate rates only (not initial conditions) from wild guesses
    opts.estimate_params = [p for p in model.parameters if not p.name.endswith('_0') ]
    opts.initial_values = [p.value for p in opts.estimate_params]
    opts.nsteps = 10000
    opts.likelihood_fn = likelihood
    opts.step_fn = step
    opts.use_hessian = True
    opts.hessian_period = opts.nsteps / 10 # Calculate the Hessian 10 times
    opts.seed = seed
    mcmc = mcmc_hessian.MCMC(opts)

    # Run it!
    mcmc.run()

    # Pickle it!
    output_file = open('nbd62c_m1c_10k.pck', 'w')
    pickle.dump(mcmc, output_file)
    output_file.close()

    # Plot "After" curves
    if True:
        plt.ion()
        plt.figure()
        plt.plot(tspan, ydata_norm)
        after = mcmc.simulate()
        #after_array = before.view().reshape(len(tspan), len(before.dtype))
        iBax_after = 2*after[:,6] + 4*after[:,7]
        iBax_after_norm = iBax_after / max(iBax_after)
        plt.plot(tspan, iBax_after_norm[:])
        plt.show()

    
    return mcmc

def prior(mcmc, position):
    # TODO Need to put some decent priors on the on and off rates
    mean = math.log10([1e-2, 1e7])
    var = [100, 100]
    return numpy.sum((position - means) ** 2 / ( 2 * var))

def step(mcmc):
    """The function to call at every iteration. Currently just prints
    out a few progress indicators.
    """
    if mcmc.iter % 20 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  acc=%-.3f, lkl=%g  prior=%g  post=%g' % \
            (mcmc.iter, mcmc.sig_value, mcmc.T, mcmc.acceptance/(mcmc.iter+1),
             mcmc.accept_likelihood, mcmc.accept_prior, mcmc.accept_posterior)


# Plotting
# ========

def plot_from_params(params_dict):
    """Plot the model output using the given parameter value."""

    # The specific model may need to be changed here
    m1c = tBid_Bax_1c(params_dict=params_dict)
    m1c.build_model2()
    model = m1c.model

    plt.ion()

    plt.plot(tspan, ydata_norm)
    output = odesolve(model, tspan)
    #output_array = output.view().reshape(len(tspan), len(output.dtype))
    iBax = 2*output['Bax2'] + 4*output['Bax4']
    iBax_norm = iBax / max(iBax)
    plt.plot(tspan, iBax_norm[:])

    plt.show()

"""
Notes

Basic version:

We fit the kinetic parameters in the model and make assumptions about the observables
and their relation to the NBD signals (e.g., we can first do Bax2 for the Bax62C data,
and compare it to a case where Bax62C is tBid/Bax driven).

So--need to load the data (the three curves, and then normalize it to 0->1)
Then run the model, fitting only the kinetic parameters (not the initial conditions),
evaluate the objective function over the timecourse. Use a figure for the error based on the
variance in the data in a relatively straight-line area, and/or over a short distance.

Build out the ODE/core model with Bax2 binding and Bax4. So will have many parameters...

Advanced version

Ideally, would have a way of pre-equilibrating the system for the just-Bax condition,
and then perturb it with the addition of tBid.

Could develop more comprehensive/enumerated set of hypotheses where the Bax binding was due
to other states
"""
