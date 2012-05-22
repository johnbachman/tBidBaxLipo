import mcmc_hessian
from pysb.integrate import odesolve
import numpy
import matplotlib.pyplot as plt

from tBid_Bax_1c import tBid_Bax_1c

import nbd_analysis


m1c = tBid_Bax_1c()
m1c.build_model0()
model = m1c.model

seed = 2
random = numpy.random.RandomState(seed)

#ntimes = 20; # TODO Number of timepoints?
#tspan = numpy.linspace(0, 40, ntimes); # TODO define tspan
tspan = nbd_analysis.time_other
nbd_avgs, nbd_stds = nbd_analysis.calc_norm_avg_std()
ydata_norm = nbd_avgs[2] # 120C

#ysim = odesolve(model, tspan)
#ysim_array = ysim.view().reshape(len(tspan), len(ysim.dtype))
#yspecies = ysim_array[:, :len(model.species)]

sigma = 0.05; # TODO Set sigma appropriately
#ydata = yspecies * (random.randn(*yspecies.shape) * sigma + 1);
#ysim_max = yspecies.max(0)
#ydata_norm = ydata / ysim_max

def likelihood(mcmc, position):
    yout = mcmc.simulate(position)

    # TODO need a way to get the desired observable(s) from the array
    iBax = yout[:,5]
    iBax_max = max(iBax)
    yout_norm = iBax / iBax_max    

    return numpy.sum((ydata_norm - yout_norm) ** 2 / (2 * sigma ** 2))

#def prior(mcmc, position):
#    mean = math.log10([1e-2, 1e7])
#    var = [100, 100]
#    return numpy.sum((position - means) ** 2 / ( 2 * var))

def step(mcmc):
    if mcmc.iter % 20 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  acc=%-.3f, lkl=%g  prior=%g  post=%g' % \
            (mcmc.iter, mcmc.sig_value, mcmc.T, mcmc.acceptance/(mcmc.iter+1),
             mcmc.accept_likelihood, mcmc.accept_prior, mcmc.accept_posterior)

opts = mcmc_hessian.MCMCOpts()
opts.model = model
opts.tspan = tspan

print model.parameters
# estimate rates only (not initial conditions) from wild guesses
opts.estimate_params = [p for p in model.parameters if not p.name.endswith('_0') ]
opts.initial_values = [p.value for p in opts.estimate_params]
#opts.initial_values = [1e-4, 1e3, 1e6] # TODO set to nominal values?

opts.nsteps = 1000
opts.likelihood_fn = likelihood
opts.step_fn = step
opts.use_hessian = True
opts.hessian_period = opts.nsteps / 10
opts.seed = seed
mcmc = mcmc_hessian.MCMC(opts)


# Plot "Before" curves

mcmc.run()



"""
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
