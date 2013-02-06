"""Functions for fitting the linear NBD insertion models to the NBD-Bax mutant
fluorescence data using MCMC.

.. todo:: Likelihood should be scaled by number of timepoints, otherwise it
overwhelms the prior for "no reason".
"""

import bayessb
from pysb.integrate import odesolve
import numpy
import matplotlib.pyplot as plt
import nbd_analysis as nbd
import pickle
from tbidbaxlipo.util.report import Report
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

# MCMC Functions
# ==============

def do_fit(model, prior_func=None, estimate_params=None,
           initial_values=None, basename='nbd_mcmc',
           random_seed=1, nsteps=2000):
    """Runs MCMC on the globally defined model."""

    # Initialize the MCMC arguments
    opts = bayessb.MCMCOpts()
    opts.model = model
    opts.tspan = nbd.time_other

    # Set the params to estimate
    if estimate_params is None:
        opts.estimate_params = [p for p in model.parameters
                                  if not p.name.endswith('_0')]
    else:
        opts.estimate_params = estimate_params

    # Set the initial values for the parameters
    if initial_values is not None:
        opts.initial_values = initial_values
    else:
        opts.initial_values = [p.value for p in opts.estimate_params]

    opts.nsteps = nsteps
    opts.likelihood_fn = likelihood
    opts.step_fn = step
    opts.prior_fn = prior_func

    opts.sigma_adj_interval = 20
    opts.sigma_step = 0.9
    opts.sigma_max = 100
    opts.sigma_min = 0.005

    opts.use_hessian = True
    opts.hessian_scale = 1
    opts.hessian_period = opts.nsteps / 5 #10

    opts.seed = random_seed
    mcmc = bayessb.MCMC(opts)

    mcmc.initialize()

    # Run it!
    mcmc.run()

    # Pickle it!
    mcmc.options.likelihood_fn = None
    mcmc.options.prior_fn = None
    mcmc.options.step_fn = None
    output_file = open('%s.pck' % basename, 'w')
    pickle.dump(mcmc, output_file)
    output_file.close()
    mcmc.options.likelihood_fn = likelihood
    mcmc.options.prior_fn = prior_func

    return mcmc

def generate_figures(mcmc, do_report=True, basename='nbd_mcmc'):   
    """Takes an MCMC chain and plots a series of useful visualizations of the
    walk, the quality of the fit, etc.
    """

    if do_report:
        rep = Report()

    # Plot "Before" curves -------
    plt.ion()
    plt.figure()
    plot_data()

    initial_params = mcmc.cur_params(position=mcmc.initial_position)

    # TODO TODO TODO Set the curves to plot here
    x = mcmc.simulate(position=mcmc.initial_position, observables=True)
    plt.plot(mcmc.options.tspan,
             x['Baxc3'] * initial_params[3], 'r', label='c3 model')
    #plt.plot(tspan, x['Baxc62'] * initial_params[2], 'g', label='c62 model')
    #plt.plot(tspan, x['Baxc120'] * initial_params[3], 'b', label='c120 model')
    #plt.plot(tspan, x['Baxc122'] * initial_params[4], 'm', label='c122 model')
    #plt.plot(tspan, x['Baxc126'] * initial_params[5], 'k', label='c126 model')
    plt.legend(loc='lower right')
    plt.title('Before')
    plt.show()

    if do_report:
        rep.add_current_figure()

    p_name_vals = zip([p.name for p in mcmc.options.model.parameters],
                      initial_params)
    print "Initial values:"
    print('\n'.join(['%s: %g' % (p_name_vals[i][0], p_name_vals[i][1])
                     for i in range(0, len(p_name_vals))]))

    # Plot "After" curves ------------

    # Set to best fit position
    #best_fit_position = mcmc.positions[numpy.argmin(mcmc.posteriors)]

    # Set to last fit position
    mixed_start = mcmc.options.nsteps / 2
    #mixed_start = 30000
    mixed_positions = mcmc.positions[mixed_start:,:]
    mixed_accepted_positions = mixed_positions[mcmc.accepts[mixed_start:]]
    best_fit_position = mixed_accepted_positions[-1,:]

    plt.figure()
    plot_data()

    best_fit_params = mcmc.cur_params(position=best_fit_position)

    x = mcmc.simulate(position=best_fit_position, observables=True)

    # TODO TODO TODO Set the curves to plot here
    plt.plot(mcmc.options.tspan,
             x['Baxc3'] * best_fit_params[3], 'r', label='c3 model')
    #plt.plot(tspan, x['Baxc62'] * best_fit_params[2], 'g', label='c62 model')
    #plt.plot(tspan, x['Baxc120'] * best_fit_params[3], 'b', label='c120 model')
    #plt.plot(tspan, x['Baxc122'] * best_fit_params[4], 'm', label='c122 model')
    #plt.plot(tspan, x['Baxc126'] * best_fit_params[5], 'k', label='c126 model')
    plt.legend(loc='lower right')
    plt.title('Final accepted position')

    plt.show()

    if do_report:
        rep.add_current_figure()

    # -- Print final step parameter values --------------------
    p_name_vals = zip([p.name for p in mcmc.options.model.parameters],
                       best_fit_params)
    p_name_vals_string = '\n'.join(['%s: %g' % (p_name_vals[i][0],
                                                p_name_vals[i][1])
                                    for i in range(0, len(p_name_vals))])
    print p_name_vals_string

    if do_report:
        rep.add_text(p_name_vals_string)
   
    # Plot sampling of fits ----------------------------------
    plt.figure()
    plot_data()

    max_position_index = len(mixed_accepted_positions) - 1
    num_to_plot = min(500, max_position_index)
    for i in range(num_to_plot):
        cur_position = mixed_accepted_positions[max_position_index - i,:]
        x = mcmc.simulate(position=cur_position, observables=True)
        cur_position_linear = mcmc.cur_params(cur_position)
        # FIXME magic number 3!!!)
        plt.plot(mcmc.options.tspan,
                 x['Baxc3'] * cur_position_linear[3], alpha=0.02, color='r')

    plt.title("Sampling of %d accepted positions" % num_to_plot)
    plt.show()

    if do_report:
        rep.add_current_figure()

    # Plot marginal distributions of parameters ---------------
    for i, cur_param in enumerate(mcmc.options.estimate_params):
        plt.figure()
        plt.hist(mixed_accepted_positions[:,i])
        plt.title("%s, last %d accepts: initval %f" %
                (cur_param.name, len(mixed_accepted_positions[:,i]),
                 cur_param.value))
        plt.show()

        if do_report:
            rep.add_current_figure()

    # Plot convergence traces of all parameters
    plt.figure()
    plt.plot(mcmc.positions)
    plt.title("Parameter traces")
    plt.legend([p.name for p in mcmc.options.estimate_params], loc='lower left',
               prop={'size':6})
    plt.show()

    if do_report:
        rep.add_current_figure()

    if do_report:
        # Add the code for the fitting (this file)
        rep.add_python_code('nbd_mcmc_c3.py')
       
        # Add the code for the model
        rep.add_python_code('models/core.py')
        rep.add_python_code('models/one_cpt.py')

        # Write the report 
        rep.write_report(basename)

def plot_data():
    alpha = 0.5
    plt.plot(nbd.time_other, nbd_avgs[0], 'r.', label='c3 data', alpha=alpha)
    """
    plt.plot(nbd.time_other, nbd_avgs[1], 'g.', label='c62 data', alpha=alpha)
    plt.plot(nbd.time_other, nbd_avgs[2], 'b.', label='c120 data', alpha=alpha)
    plt.plot(nbd.time_other, nbd_avgs[3], 'm.', label='c122 data', alpha=alpha)
    plt.plot(nbd.time_other, nbd_avgs[4], 'k.', label='c126 data', alpha=alpha)
    """

# Define the likelihood function
def likelihood(mcmc, position):
    """The likelihood function. Calculates the error between model and data
    to give a measure of the likelihood of observing the data given the
    current parameter values.
    """
    yout = mcmc.simulate(position, observables=True)

    params = mcmc.cur_params(position)
    #print params

    # TODO Need to be able to get the indices from the model so that
    # they're not hardcoded
    c3_scaling   = params[3]
    #c62_scaling  = params[2]
    #c120_scaling = params[3]
    #c122_scaling = params[4]
    #c126_scaling = params[5]

    c3_model = yout['Baxc3'] * c3_scaling
    #c62_model = yout['Baxc62'] * c62_scaling
    #c120_model = yout['Baxc120'] * c120_scaling
    #c122_model = yout['Baxc122'] * c122_scaling
    #c126_model = yout['Baxc126'] * c126_scaling

    # TODO TODO TODO Set the components of the obj func here
    err  = numpy.sum((nbd_avgs[0] - c3_model)**2 / (2 * nbd_stds[0]**2))
    #err += numpy.sum((nbd_avgs[1] - c62_model)**2 / (2 * nbd_stds[1]**2))
    #err += numpy.sum((nbd_avgs[2] - c120_model)**2 / (2 * nbd_stds[2]**2))
    #err += numpy.sum((nbd_avgs[3] - c122_model)**2 / (2 * nbd_stds[3]**2))
    #err += numpy.sum((nbd_avgs[4] - c126_model)**2 / (2 * nbd_stds[4]**2))

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
                      basename=sys.argv[3], random_seed=index)
    else:
        # TODO TODO TODO Set num_residues and basename here
        print("Running do_fit() with the default arguments...")

        from tbidbaxlipo.models import one_cpt
        b = one_cpt.Builder()

        # -- select model ---
        #b.build_model_ta()

        #b.build_model_tar()

        b.build_model_tai()
        # -----

        random_seed = 1
        nsteps = 40000
        basename = '%s_%d_steps_seed_%d' % ('nbd_mcmc', nsteps, random_seed)

        mcmc = do_fit(b.model, b.prior, estimate_params=b.estimate_params,
                      basename=basename, random_seed=random_seed, nsteps=nsteps)

        generate_figures(mcmc, basename=basename)

        #initial_values = random_initial_values(num_sets=3, num_residues=1)
        #mcmc = do_fit(initial_values=initial_values[0]) # Run with the defaults
        #initial_values = [0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1]
        #mcmc = do_fit(initial_values=initial_values[0], random_seed=1,
        #        basename='nbd_mcmc_c3')
        #mcmc = do_fit(initial_values=initial_values[1], random_seed=2,
        #        basename='nbd_mcmc_c3')
        # mcmc = do_fit(initial_values=initial_values[2], random_seed=3,
        #       basename='nbd_mcmc_c3')

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
