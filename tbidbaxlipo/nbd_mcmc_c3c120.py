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
#from nbd_model import model
#import nbd_parallel_model
from nbd_parallel_model import model, prior, random_initial_values
from scipy.interpolate import interp1d
import sys

# Prepare the data
# ================

tspan = nbd.time_other
nbd_avgs, nbd_stds = nbd.calc_norm_avg_std()

# -- c62 calculation --
#nbd62_f = numpy.vectorize(interp1d(nbd.time_c62, nbd_avgs[1],
#                                 bounds_error=False))
#nbd62_interp = nbd62_f(nbd.time_other)
#nbd_avgs[1] = nbd62_interp
nbd_avgs[1] = nbd_avgs[1][0:-2]
nbd_stds[1] = nbd_stds[1][0:-2]

# FIXME FIXME FIXME 
# A temporary hack to see if downsampling the timecourse improves the
# objective function landscape
#for i in range(0, len(nbd_avgs)):
#    nbd_avgs[i] = nbd_avgs[i][::100]
#    nbd_stds[i] = nbd_stds[i][::100]
likelihood_matrix = []
#nbd.time_other = nbd.time_other[::100]
#tspan = nbd.time_other
# FIXME FIXME FIXME

# MCMC Functions
# ==============

def do_fit(initial_values=None, basename='nbd_mcmc', random_seed=1):
    """Runs MCMC on the globally defined model."""

    # Initialize the MCMC arguments
    opts = bayessb.MCMCOpts()
    opts.model = model
    opts.tspan = nbd.time_other

    print model.parameters
    # estimate rates only (not initial conditions) from wild guesses
    #opts.estimate_params = [p for p in model.parameters
    #                          if not p.name.endswith('_0')]
    # TODO Set the params to fit here TODO TODO TODO
    opts.estimate_params = [model.parameters['c3_scaling'],
                            model.parameters['c120_scaling'],
                            model.parameters['c3_insertion_rate'],
                            model.parameters['c120_insertion_rate']]

    if initial_values is not None:
        opts.initial_values = initial_values
    else:
        opts.initial_values = [p.value for p in opts.estimate_params]

    opts.nsteps = 1000 #2000
    opts.likelihood_fn = likelihood
    #opts.prior_fn = prior
    opts.step_fn = step
    opts.use_hessian = True #True
    opts.hessian_period = opts.nsteps / 10 #10 # Calculate the Hessian 10 times
    opts.seed = random_seed
    mcmc = bayessb.MCMC(opts)

    mcmc.initialize()
    
    # Plot "Before" curves -------
    plt.ion()
    plt.figure()
    plot_data()

    initial_params = mcmc.cur_params(position=mcmc.initial_position)

    # TODO TODO TODO Set the curves to plot here
    x = mcmc.simulate(position=mcmc.initial_position, observables=True)
    plt.plot(tspan, x['Baxc3'] * initial_params[1], 'r', label='c3 model')
    #plt.plot(tspan, x['Baxc62'] * initial_params[2], 'g', label='c62 model')
    plt.plot(tspan, x['Baxc120'] * initial_params[3], 'b', label='c120 model')
    #plt.plot(tspan, x['Baxc122'] * initial_params[4], 'm', label='c122 model')
    #plt.plot(tspan, x['Baxc126'] * initial_params[5], 'k', label='c126 model')
    plt.legend(loc='lower right')
    plt.title('Before')
    plt.show()

    p_name_vals = zip([p.name for p in model.parameters], initial_params)
    print('\n'.join(['%s: %g' % (p_name_vals[i][0], p_name_vals[i][1])
                     for i in range(0, len(p_name_vals))]))
    # ---------------------------

    # Run it!
    mcmc.run()

    # Pickle it!
    output_basename = '%s_%d_steps_seed_%d' % \
                      (basename, opts.nsteps, random_seed)
    mcmc.options.likelihood_fn = None
    mcmc.options.step_fn = None
    output_file = open('%s.pck' % output_basename, 'w')
    pickle.dump(mcmc, output_file)
    output_file.close()


    # Plot "After" curves ------------
    # Set to best fit position
    best_fit_position = mcmc.positions[numpy.argmin(mcmc.posteriors)]

    plt.figure()
    plot_data()

    best_fit_params = mcmc.cur_params(position=best_fit_position)

    x = mcmc.simulate(position=best_fit_position, observables=True)

    # TODO TODO TODO Set the curves to plot here
    plt.plot(tspan, x['Baxc3'] * best_fit_params[1], 'r', label='c3 model')
    #plt.plot(tspan, x['Baxc62'] * best_fit_params[2], 'g', label='c62 model')
    plt.plot(tspan, x['Baxc120'] * best_fit_params[3], 'b', label='c120 model')
    #plt.plot(tspan, x['Baxc122'] * best_fit_params[4], 'm', label='c122 model')
    #plt.plot(tspan, x['Baxc126'] * best_fit_params[5], 'k', label='c126 model')
    plt.legend(loc='lower right')
    plt.title('After')

    plt.show()

    rep = Report()
    rep.addCurrentFigure()
    rep.writeReport(output_basename)

    p_name_vals = zip([p.name for p in model.parameters], best_fit_params)
    print('\n'.join(['%s: %g' % (p_name_vals[i][0], p_name_vals[i][1])
                     for i in range(0, len(p_name_vals))]))

    return mcmc

def plot_data():
    alpha = 0.5
    plt.plot(nbd.time_other, nbd_avgs[0], 'r.', label='c3 data', alpha=alpha)
    plt.plot(nbd.time_other, nbd_avgs[1], 'g.', label='c62 data', alpha=alpha)
    plt.plot(nbd.time_other, nbd_avgs[2], 'b.', label='c120 data', alpha=alpha)
    plt.plot(nbd.time_other, nbd_avgs[3], 'm.', label='c122 data', alpha=alpha)
    plt.plot(nbd.time_other, nbd_avgs[4], 'k.', label='c126 data', alpha=alpha)


# Define the likelihood function
def likelihood(mcmc, position):
    """The likelihood function. Calculates the error between model and data
    to give a measure of the likelihood of observing the data given the
    current parameter values.
    """
    yout = mcmc.simulate(position, observables=True)

    params = mcmc.cur_params(position)
    #print params

    c3_scaling   = params[1]
    c62_scaling  = params[2]
    c120_scaling = params[3]
    c122_scaling = params[4]
    c126_scaling = params[5]

    c3_model = yout['Baxc3'] * c3_scaling
    c62_model = yout['Baxc62'] * c62_scaling
    c120_model = yout['Baxc120'] * c120_scaling
    c122_model = yout['Baxc122'] * c122_scaling
    c126_model = yout['Baxc126'] * c126_scaling

    # TODO TODO TODO Set the components of the obj func here
    err  = numpy.sum((nbd_avgs[0] - c3_model)**2 / (2 * nbd_stds[0]**2))
    #err += numpy.sum((nbd_avgs[1] - c62_model)**2 / (2 * nbd_stds[1]**2))
    err += numpy.sum((nbd_avgs[2] - c120_model)**2 / (2 * nbd_stds[2]**2))
    #err += numpy.sum((nbd_avgs[3] - c122_model)**2 / (2 * nbd_stds[3]**2))
    #err += numpy.sum((nbd_avgs[4] - c126_model)**2 / (2 * nbd_stds[4]**2))

    """
    #err  = numpy.sum((nbd_avgs[0][::100] - c3_model[::100])**2 /
    #        (2 * nbd_stds[0][::100] **2))
    #err  += numpy.sum((nbd_avgs[1][::100] - c62_model[::100])**2 /
    #                 (2 * nbd_stds[1][::100]**2))
    #err  += numpy.sum((nbd_avgs[2][::100] - c120_model[::100])**2 /
    #                 (2 * nbd_stds[2][::100]**2))
    #err  += numpy.sum((nbd_avgs[3][::100] - c122_model[::100])**2 /
    #                 (2 * nbd_stds[3][::100]**2))
    #err  += numpy.sum((nbd_avgs[4][::100] - c126_model[::100])**2 /
    #                 (2 * nbd_stds[4][::100]**2))

    err  = numpy.sum((nbd_avgs[0][::100] - c3_model[::100])**2 /
            (2 * 0.1 **2))
    err  += numpy.sum((nbd_avgs[1][::100] - c62_model[::100])**2 /
            (2 * 0.1 **2))
    err  += numpy.sum((nbd_avgs[2][::100] - c120_model[::100])**2 /
            (2 * 0.1 **2))
    err  += numpy.sum((nbd_avgs[3][::100] - c122_model[::100])**2 /
            (2 * 0.1 **2))
    err  += numpy.sum((nbd_avgs[4][::100] - c126_model[::100])**2 /
            (2 * 0.1 **2))

    c3err  = numpy.sum((nbd_avgs[0][::100] - c3_model[::100])**2 /
            (2 * nbd_stds[0][::100] **2))
    c62err  = numpy.sum((nbd_avgs[1][::100] - c62_model[::100])**2 /
                     (2 * nbd_stds[1][::100]**2))
    c120err  = numpy.sum((nbd_avgs[2][::100] - c120_model[::100])**2 /
                     (2 * nbd_stds[2][::100]**2))
    c122err  = numpy.sum((nbd_avgs[3][::100] - c122_model[::100])**2 /
                     (2 * nbd_stds[3][::100]**2))
    c126err  = numpy.sum((nbd_avgs[4][::100] - c126_model[::100])**2 /
                     (2 * nbd_stds[4][::100]**2))
    likelihood_row = [c3err, c62err, c120err, c122err, c126err]

    likelihood_matrix.append(likelihood_row)
    return numpy.sum(numpy.array(likelihood_row))
    """
    return err

def step(mcmc):
    """The function to call at every iteration. Currently just prints
    out a few progress indicators.
    """
    if mcmc.iter % 20 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  acc=%-.3f, lkl=%g  prior=%g  post=%g' % \
            (mcmc.iter, mcmc.sig_value, mcmc.T, mcmc.acceptance/(mcmc.iter+1.),
             mcmc.accept_likelihood, mcmc.accept_prior, mcmc.accept_posterior)

# Main function
# =============

if __name__ == '__main__':
    if (len(sys.argv) == 4):
        # The first command line argument should contain the name of the
        # file to unpickle to get the list of starting values
        input_file = open(sys.argv[1], 'r')
        initial_value_list = pickle.load(input_file)

        # The second command line argument contains the index into the
        # matrix of starting values identifying the vector of values to use.
        # This value is also used here as the random seed.
        index = int(sys.argv[2])
        mcmc = do_fit(initial_values=initial_value_list[index],
                      basename=sys.argv[3], random_seed=index)
    else:
        print("Running do_fit() with the default arguments...")
        initial_values = random_initial_values(num_sets=3, num_residues=2)
        #mcmc = do_fit(initial_values=initial_values[0]) # Run with the defaults
        #initial_values = [0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1]
        mcmc = do_fit(initial_values=initial_values[0], random_seed=1,
                basename='nbd_mcmc_c3n120')
        mcmc = do_fit(initial_values=initial_values[1], random_seed=2,
                basename='nbd_mcmc_c3n120')
        mcmc = do_fit(initial_values=initial_values[2], random_seed=3,
                basename='nbd_mcmc_c3n120')
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
