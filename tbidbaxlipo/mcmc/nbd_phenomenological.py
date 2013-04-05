"""
Functions for fitting the linear or parallel NBD insertion models to the
NBD-Bax mutant fluorescence data using MCMC.
"""

import bayessb
from pysb.integrate import odesolve
import numpy
import matplotlib.pyplot as plt
import tbidbaxlipo.nbd_analysis as nbd
import pickle
from tbidbaxlipo.util.report import Report
#from tbidbaxlipo.models.nbd_parallel_model 
#     import model, prior, random_initial_values
from tbidbaxlipo.models.nbd_linear_model \
     import model, random_initial_values, prior
from scipy.interpolate import interp1d
import sys

# Prepare the data
# ================

tspan = nbd.time_other
nbd_avgs, nbd_stds = nbd.calc_avg_std()

# Cut off the last two timepoints from C62 so that it's the same length
nbd_avgs[1] = nbd_avgs[1][0:-2]
nbd_stds[1] = nbd_stds[1][0:-2]

# MCMC Functions
# ==============

def do_fit(initial_values=None, basename='nbd_mcmc', random_seed=1,
           show_plot=False):
    """Runs MCMC on the globally defined model."""

    # Initialize the MCMC arguments
    opts = bayessb.MCMCOpts()
    opts.model = model
    opts.tspan = nbd.time_other

    # estimate rates only (not initial conditions) from wild guesses
    opts.estimate_params = [p for p in model.parameters
                              if not p.name.endswith('_0')]

    if initial_values is not None:
        opts.initial_values = initial_values
    else:
        opts.initial_values = [p.value for p in opts.estimate_params]

    opts.nsteps = 4000 #2000
    opts.likelihood_fn = likelihood
    opts.prior_fn = prior
    opts.step_fn = step

    opts.use_hessian = True
    #opts.anneal_length = 4000
    opts.hessian_scale = 1
    opts.hessian_period = opts.nsteps / 10 #10 # Calculate the Hessian 10 times

    opts.sigma_adj_interval = 20
    opts.sigma_step = 0.9
    opts.sigma_max = 100
    opts.sigma_min = 0

    #opts.norm_step_size = 1.
    opts.seed = random_seed
    mcmc = bayessb.MCMC(opts)

    mcmc.initialize()

    # Plot "Before" curves -------
    plt.ion()
    plt.figure()
    plot_data()

    initial_params = mcmc.cur_params(position=mcmc.initial_position)

    # Set the curves to plot here
    x = mcmc.simulate(position=mcmc.initial_position, observables=True)
    plt.plot(tspan, x['Baxc3'] * initial_params[1], 'r', label='c3 model')
    plt.plot(tspan, x['Baxc62'] * initial_params[2], 'g', label='c62 model')
    plt.plot(tspan, x['Baxc120'] * initial_params[3], 'b', label='c120 model')
    plt.plot(tspan, x['Baxc122'] * initial_params[4], 'm', label='c122 model')
    plt.plot(tspan, x['Baxc126'] * initial_params[5], 'k', label='c126 model')
    plt.legend(loc='lower right')
    plt.title('Before')
    if show_plot:
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

    # Set the curves to plot here
    plt.plot(tspan, x['Baxc3'] * best_fit_params[1], 'r', label='c3 model')
    plt.plot(tspan, x['Baxc62'] * best_fit_params[2], 'g', label='c62 model')
    plt.plot(tspan, x['Baxc120'] * best_fit_params[3], 'b', label='c120 model')
    plt.plot(tspan, x['Baxc122'] * best_fit_params[4], 'm', label='c122 model')
    plt.plot(tspan, x['Baxc126'] * best_fit_params[5], 'k', label='c126 model')
    plt.legend(loc='lower right')
    plt.title('After')
    if show_plot:
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

    err  = numpy.sum((nbd_avgs[0] - c3_model)**2 / (2 * nbd_stds[0]**2))
    err += numpy.sum((nbd_avgs[1] - c62_model)**2 / (2 * nbd_stds[1]**2))
    err += numpy.sum((nbd_avgs[2] - c120_model)**2 / (2 * nbd_stds[2]**2))
    err += numpy.sum((nbd_avgs[3] - c122_model)**2 / (2 * nbd_stds[3]**2))
    err += numpy.sum((nbd_avgs[4] - c126_model)**2 / (2 * nbd_stds[4]**2))

    return err

def step(mcmc):
    """The function to call at every iteration. Currently just prints
    out a few progress indicators.
    """
    window = 200.
    if mcmc.iter % 20 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  loc_acc=%-.3f  ' \
              'glob_acc=%-.3f  lkl=%g  prior=%g  post=%g' % \
              (mcmc.iter, mcmc.sig_value, mcmc.T,
               numpy.sum(mcmc.accepts[mcmc.iter - window:mcmc.iter+1]) /
               float(window),
               mcmc.acceptance/(mcmc.iter+1.), mcmc.accept_likelihood,
               mcmc.accept_prior, mcmc.accept_posterior)

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
                      basename=sys.argv[3], random_seed=index, show_plot=False)
    else:
        print("Running do_fit() with the default arguments...")

        # Run the pathological case
        #input_file = open('nbd_mcmc_parallel_random_initial_values.pck')
        #initial_value_list = pickle.load(input_file)
        #mcmc11 = do_fit(initial_values=initial_value_list[11], random_seed=11,
        #                basename='nbd_mcmc_parallel11_gausspr')

        numpy.random.seed(1)
        initial_values = random_initial_values(num_sets=3)
        mcmc = do_fit(initial_values=initial_values[0], show_plot=True,
                      random_seed=1) # Run with the defaults
