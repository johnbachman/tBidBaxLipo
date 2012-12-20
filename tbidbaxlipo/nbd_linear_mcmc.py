"""Functions for fitting the linear NBD insertion models to the NBD-Bax mutant
fluorescence data using MCMC.
"""

import bayessb
from pysb.integrate import odesolve
import numpy
import matplotlib.pyplot as plt
import nbd_analysis as nbd
import pickle
from tbidbaxlipo.util.report import Report
from nbd_model import model
from scipy.interpolate import interp1d

# Prepare the data
# ================

tspan = nbd.time_other
nbd_avgs, nbd_stds = nbd.calc_norm_avg_std()

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
        yout = mcmc.simulate(position, observables=True)

        params = mcmc.cur_params(position)
        
        c3_scaling   = params[1]
        c62_scaling  = params[2]
        c120_scaling = params[3]
        c122_scaling = params[4]
        c126_scaling = params[5]

        c3_model = yout['Baxc3'] * c3_scaling

        c62_f = numpy.vectorize(interp1d(nbd.time_other, yout['Baxc62'],
                                         bounds_error=False))
        # Only use c62 up until 3589 seconds so that we don't get a bounds
        # error (nbd.time_other only goes up to 3590, whereas nbd.time_c2
        # goes up to 3594.7 seconds)
        c62_interp = c62_f(nbd.time_c62[0:1794])
        c62_model = c62_interp * c62_scaling

        c120_model = yout['Baxc120'] * c120_scaling
        c122_model = yout['Baxc122'] * c122_scaling
        c126_model = yout['Baxc126'] * c126_scaling

        err  = numpy.sum((nbd_avgs[0] - c3_model)**2 / (2 * nbd_stds[0]**2))
        err += numpy.sum((nbd_avgs[1][0:1794] - c62_model)**2 /
                (2 * nbd_stds[1][0:1794]**2))
        err += numpy.sum((nbd_avgs[2] - c120_model)**2 / (2 * nbd_stds[2]**2))
        err += numpy.sum((nbd_avgs[3] - c122_model)**2 / (2 * nbd_stds[3]**2))
        err += numpy.sum((nbd_avgs[4] - c126_model)**2 / (2 * nbd_stds[4]**2))

        return err

    # Set the random number generator seed
    seed = 4
    random = numpy.random.RandomState(seed)

    # Initialize the MCMC arguments
    opts = bayessb.MCMCOpts()
    opts.model = model
    opts.tspan = nbd.time_other

    print model.parameters
    # estimate rates only (not initial conditions) from wild guesses
    opts.estimate_params = [p for p in model.parameters
                              if not p.name.endswith('_0')]
                                 
    opts.initial_values = [p.value for p in opts.estimate_params]
    opts.nsteps = 
    opts.likelihood_fn = likelihood
    opts.step_fn = step
    opts.use_hessian = False
    opts.hessian_period = opts.nsteps / 10 # Calculate the Hessian 10 times
    opts.seed = seed
    mcmc = bayessb.MCMC(opts)

    # Run it!
    mcmc.run()

    # Pickle it!
    #mcmc.solver = None # FIXME This is a hack to make the MCMC pickleable
    #output_file = open('nbd62c_m1c_10k.pck', 'w')
    #pickle.dump(mcmc, output_file)
    #output_file.close()

    # Set to best fit position
    mcmc.position = mcmc.positions[numpy.argmin(mcmc.likelihoods)]

    # Plot "After" curves
    plt.ion()
    plt.figure()
    plot_data()

    best_fit_params = mcmc.cur_params(position=mcmc.position)

    x = odesolve(model, tspan)
    plt.plot(tspan, x['Baxc3'] * best_fit_params[1], 'r', label='c3 model')
    plt.plot(tspan, x['Baxc62'] * best_fit_params[2], 'g', label='c62 model')
    plt.plot(tspan, x['Baxc120'] * best_fit_params[3], 'b', label='c120 model')
    plt.plot(tspan, x['Baxc122'] * best_fit_params[4], 'm', label='c122 model')
    plt.plot(tspan, x['Baxc126'] * best_fit_params[5], 'k', label='c126 model')
    plt.show()
    plt.legend(loc='lower right')

    plt.show()

    rep = Report()
    rep.addCurrentFigure()
    rep.writeReport('nbd_mcmc_fit')

    p_name_vals = zip([p.name for p in model.parameters], best_fit_params)
    print('\n'.join(['%s: %g' % (p_name_vals[i][0], p_name_vals[i][1])
                     for i in range(0, len(p_name_vals))]))

    return mcmc

def plot_data():
    alpha = 0.5
    plt.plot(nbd.time_other, nbd_avgs[0], 'r.', label='c3 data', alpha=alpha)
    plt.plot(nbd.time_c62, nbd_avgs[1], 'g.', label='c62 data', alpha=alpha)
    plt.plot(nbd.time_other, nbd_avgs[2], 'b.', label='c120 data', alpha=alpha)
    plt.plot(nbd.time_other, nbd_avgs[3], 'm.', label='c122 data', alpha=alpha)
    plt.plot(nbd.time_other, nbd_avgs[4], 'k.', label='c126 data', alpha=alpha)

def prior(mcmc, position):
    # TODO Need to put some decent priors on the on and off rates
    #mean = math.log10([1e-2, 1e7])
    #var = [100, 100]
    #return numpy.sum((position - means) ** 2 / ( 2 * var))
    pass

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
    m1c = one_cpt.Builder(params_dict=params_dict)
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


# Main function
# =============

if __name__ == '__main__':
    mcmc = do_fit()



"""
Notes

Basic version:

We fit the kinetic parameters in the model and make assumptions about the
observables and their relation to the NBD signals (e.g., we can first do Bax2
for the Bax62C data, and compare it to a case where Bax62C is tBid/Bax driven).

So--need to load the data (the three curves, and then normalize it to 0->1)
Then run the model, fitting only the kinetic parameters (not the initial
conditions), evaluate the objective function over the timecourse. Use a figure
for the error based on the variance in the data in a relatively straight-line
area, and/or over a short distance.

Build out the ODE/core model with Bax2 binding and Bax4. So will have many
parameters...

Advanced version

Ideally, would have a way of pre-equilibrating the system for the just-Bax
condition, and then perturb it with the addition of tBid.

Could develop more comprehensive/enumerated set of hypotheses where the Bax
binding was due to other states """
