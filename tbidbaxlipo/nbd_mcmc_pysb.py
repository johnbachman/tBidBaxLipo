"""Functions for fitting the linear NBD insertion models to the NBD-Bax mutant
fluorescence data using MCMC.

.. todo:: Likelihood should be scaled by number of timepoints, otherwise it
overwhelms the prior for "no reason".

.. todo:: Ideally, would have a way of pre-equilibrating the system for the
just-Bax condition, and then perturb it with the addition of tBid.
"""

import bayessb
from pysb.integrate import odesolve
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbd_analysis as nbd
import pickle
from tbidbaxlipo.util.report import Report
from scipy.interpolate import interp1d
import sys

model_names = ['ta', 'tai', 'taid', 'taidt', 'tair', 'taird', 'tairdt',
               'tad', 'tadt', 'tar', 'tard', 'tardt']

nbd_site_names = ['c3', 'c62']

# Prepare the data
# ================

tspan = nbd.time_other
nbd_avgs, nbd_stds = nbd.calc_norm_avg_std()

# The c62 data is on a slightly different timescale and has two extra
# datapoints. For simplicity, we simply remove them.
nbd_avgs[1] = nbd_avgs[1][0:-2]
nbd_stds[1] = nbd_stds[1][0:-2]

# MCMC Functions
# ==============

def do_fit(model, likelihood, prior=None, estimate_params=None,
           initial_values=None, basename='nbd_mcmc',
           random_seed=1, nsteps=2000):
    """Runs MCMC on the given model."""

    # Initialize the MCMC arguments
    opts = bayessb.MCMCOpts()
    opts.model = model
    opts.tspan = nbd.time_other
    opts.estimate_params = estimate_params

    # Set the initial values for the parameters
    if initial_values is not None:
        opts.initial_values = initial_values
    else:
        opts.initial_values = [p.value for p in opts.estimate_params]

    opts.nsteps = nsteps
    opts.likelihood_fn = likelihood
    opts.step_fn = step
    if prior is not None:
        opts.prior_fn = prior

    opts.sigma_step = 0.9
    opts.sigma_max = 50
    opts.sigma_min = 0.01
    opts.accept_rate_target = 0.23
    opts.accept_window = 300
    opts.sigma_adj_interval = 300
    opts.anneal_length = nsteps / 2
    opts.use_hessian = True
    opts.hessian_scale = 1
    opts.hessian_period = opts.nsteps / 5 #10

    opts.seed = random_seed
    mcmc = bayessb.MCMC(opts)

    mcmc.initialize()

    # Print initial parameter values
    init_vals = zip([p.name for p in mcmc.options.model.parameters],
                      mcmc.cur_params(position=mcmc.initial_position))
    init_vals_str = 'Initial values:\n'
    init_vals_str += '\n'.join(['%s: %g' % (init_vals[i][0],
                                             init_vals[i][1])
                                 for i in range(0, len(init_vals))])
    print "------------------------"
    print init_vals_str
    print "------------------------"

    # Run it!
    mcmc.run()

    # Pickle it, setting functions to None:
    mcmc.options.likelihood_fn = None
    mcmc.options.prior_fn = None
    mcmc.options.step_fn = None
    output_file = open('%s.mcmc' % basename, 'w')
    pickle.dump(mcmc, output_file)
    output_file.close()
    # Restore the functions for interactive use
    mcmc.options.likelihood_fn = likelihood
    mcmc.options.prior_fn = prior
    mcmc.options.step_fn = step

    return mcmc

def nbd_timecourse(mcmc, position, nbd_observable):
    """Simulates the model at the given parameter position and returns
    the appropriately scaled timecourse for the given NBD site."""

    x = mcmc.simulate(position=position, observables=True)

    total_Bax = mcmc.options.model.parameters['Bax_0'].value
    cur_params = mcmc.cur_params(position=position)
    scaling_factor = cur_params[3]
    return (x[nbd_observable] / total_Bax) * scaling_factor 

def generate_figures(mcmc, nbd_site, nbd_observable, do_report=True,
                     mixed_start=None, basename='nbd_mcmc', num_samples=500):   
    """Takes an MCMC chain and plots a series of useful visualizations of the
    walk, the quality of the fit, etc.
    """
    plt.ion()

    if do_report:
        rep = Report()

    # Plot "Before" curves -------
    x = nbd_timecourse(mcmc, mcmc.initial_position, nbd_observable)
    plt.figure()
    plot_data(nbd_site)
    plt.plot(mcmc.options.tspan, x, 'b', label=nbd_site)
    plt.legend(loc='lower right')
    plt.title('Before')
    plt.show()
    if do_report:
        rep.add_current_figure()

    # Print initial parameter values
    init_vals = zip([p.name for p in mcmc.options.model.parameters],
                      mcmc.cur_params(position=mcmc.initial_position))
    init_vals_str = 'Initial values:\n'
    init_vals_str += '\n'.join(['%s: %g' % (init_vals[i][0],
                                             init_vals[i][1])
                                 for i in range(0, len(init_vals))])
    print init_vals_str
    if do_report:
        rep.add_text(init_vals_str)

    # Plot "After" curves ------------
    # Set to last fit position
    if mixed_start is None:
        mixed_start = mcmc.options.nsteps / 2
    mixed_positions = mcmc.positions[mixed_start:,:]
    mixed_accepted_positions = mixed_positions[mcmc.accepts[mixed_start:]]
    last_position = mixed_accepted_positions[-1,:]

    x = nbd_timecourse(mcmc, last_position, nbd_observable)
    plt.figure()
    plot_data(nbd_site)
    plt.plot(mcmc.options.tspan, x, 'b', label=nbd_site)
    plt.legend(loc='lower right')
    plt.title('Final accepted position')
    plt.show()
    if do_report:
        rep.add_current_figure()

    # Print final parameter values
    last_fit_params = mcmc.cur_params(position=last_position)
    last_vals = zip([p.name for p in mcmc.options.model.parameters],
                       last_fit_params)
    last_vals_str = 'Final values:\n'
    last_vals_str += '\n'.join(['%s: %g' % (last_vals[i][0],
                                            last_vals[i][1])
                                for i in range(0, len(last_vals))])
    print last_vals_str
    if do_report:
        rep.add_text(last_vals_str)

    # Plot sampling of fits # FIXME should be real sampling------
    plt.figure()
    plot_data(nbd_site)
    max_position_index = len(mixed_accepted_positions) - 1
    num_to_plot = min(num_samples, max_position_index)
    for i in range(num_to_plot):
        cur_position = mixed_accepted_positions[max_position_index - i,:]
        x = nbd_timecourse(mcmc, cur_position, nbd_observable)
        plt.plot(mcmc.options.tspan, x, 'b', label=nbd_site, alpha=0.02)
    plt.title("Last %d accepted positions" % num_to_plot)
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
               prop={'size':7})
    plt.show()
    if do_report:
        rep.add_current_figure()

    # Add code to report
    if do_report:
        # Add the code for the fitting (this file)
        rep.add_python_code('nbd_mcmc_pysb.py')
       
        # Add the code for the model
        rep.add_python_code('models/core.py')
        rep.add_python_code('models/one_cpt.py')

        # Write the report 
        rep.write_report(basename)

def plot_data(nbd_site):
    alpha = 0.5
    if nbd_site == 'c3':
        plt.plot(nbd.time_other, nbd_avgs[0], 'r.', label='c3 data',
                 alpha=alpha)
    elif nbd_site == 'c62':
        plt.plot(nbd.time_other, nbd_avgs[1], 'g.', label='c62 data',
                 alpha=alpha)
    #plt.plot(nbd.time_other, nbd_avgs[2], 'b.', label='c120 data', alpha=alpha)
    #plt.plot(nbd.time_other, nbd_avgs[3], 'm.', label='c122 data', alpha=alpha)
    #plt.plot(nbd.time_other, nbd_avgs[4], 'k.', label='c126 data', alpha=alpha)

# A function to generate the likelihood function
def get_likelihood_function(nbd_site, nbd_observable):
    """Returns a likelihood function for the specified NBD site."""
    if nbd_site == 'c3':
        data_index = 0
    elif nbd_site == 'c62':
        data_index = 1
    else:
        raise Exception('Invalid value for nbd_site!')

    def likelihood(mcmc, position):
        yout = mcmc.simulate(position, observables=True)

        # TODO Need to be able to get the indices from the model so that
        # they're not hardcoded
        params = mcmc.cur_params(position)

        timecourse = ((yout[nbd_observable] /
                       mcmc.options.model.parameters['Bax_0'].value)
                      * params[3])

        return numpy.sum((nbd_avgs[data_index] - timecourse)**2 /
                         (2 * nbd_stds[data_index]**2))

    return likelihood

def step(mcmc):
    """The function to call at every iteration. Currently just prints
    out a few progress indicators.
    """
    window = mcmc.options.accept_window

    local_acc = numpy.sum(mcmc.accepts[(mcmc.iter - window):mcmc.iter]) / \
                          float(window)

    if mcmc.iter % 20 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  loc_acc=%-.3f  ' \
              'glob_acc=%-.3f  lkl=%g  prior=%g  post=%g' % \
              (mcmc.iter, mcmc.sig_value, mcmc.T,
               local_acc,
               mcmc.acceptance/(mcmc.iter+1.), mcmc.accept_likelihood,
               mcmc.accept_prior, mcmc.accept_posterior)

# Chain handling helper function
# ==============================
def import_mcmc_groups(filenames):
    """Loads the chains into groups representing multiple runs of the same model.

    Assumes that the pickle filenames are structured as

        basename = '%s_%s_%s_%d_s%d.pck' % (model, nbd_site, nbd_observable,
                                        nsteps, random_seed)

    With the suffix ``.xxx`` separated by a dot and the seed coming last in
    the underscore-separated arguments.

    Parameter
    ---------
    filenames : list of strings
        List of strings representing the chain filenames to be sorted into
        groups, e.g., of the type returned by ``glob.glob()``.

    Returns
    -------
    dict of lists of MCMC filenames. The keys in the dict are the filename
    prefixes that represent the arguments to the MCMC procedure (e.g.,
    ``tard_c3_iBax_4000``. Each dict entry contains a list of MCMC filenames
    associated with those run conditions.
    """

    mcmc_groups = {}

    for filename in filenames:
        # Split off the suffix from the filename
        (prefix, suffix) = filename.split('.')
        # Separate the filename into the final argument identifying the
        # random seed, and everything that comes before it:
        (mcmc_args, seed) = prefix.rsplit('_', 1)
        # Check to see if we've imported another chain of this type already.
        # If so, add the current chain to the list of chains for this group:
        if mcmc_args in mcmc_groups:
            mcmc_groups[mcmc_args].append(filename)
        # If not, create a new entry in the dict containing this MCMC
        else:
            mcmc_groups[mcmc_args] = [filename]

    return mcmc_groups

# Main function
# =============

if __name__ == '__main__':
    # Set the type of model to be built here
    from tbidbaxlipo.models.one_cpt import Builder

    # Keyword args are set at the command line as e.g., key=val
    # and subsequently split at the equals sign
    kwargs = dict([arg.split('=') for arg in sys.argv[1:]])

    # We set these all to None so later on we can make sure they were
    # properly initialized.
    nbd_site = None
    nbd_observable = None
    random_seed = None
    model = None

    print "Keyword arguments: "
    print kwargs

    # Before we begin, we make sure we have all the keyword arguments that
    # we are going to need.
    if 'nbd_site' not in kwargs or \
            'random_seed' not in kwargs or \
            'nsteps' not in kwargs or \
            'model' not in kwargs or \
            'nbd_observable' not in kwargs:
        raise Exception('One or more needed arguments was not specified! ' \
                'Arguments must include nbd_site, random_seed, model, ' \
                'nbd_observable and nsteps.')

    # Because the NBD site(s) to fit has to be specified when the model
    # builder object is created, we get this arg first.
    if kwargs['nbd_site'] not in nbd_site_names:
        raise Exception('%s is not an allowed NBD site!' %
                        kwargs['nbd_site'])
    else:
        nbd_site = kwargs['nbd_site']

    # Now that we have the NBD site we can instantiate the Builder:
    builder = Builder(nbd_sites=[nbd_site])

    # The observable associated with the NBD site signal also has to be
    # specified:
    observables = [o.name for o in builder.model.observables]

    if kwargs['nbd_observable'] not in observables:
        raise Exception('%s is not an allowed NBD observable!' %
                        kwargs['nbd_observable'])
    else:
        nbd_observable = kwargs['nbd_observable']

    # ...and then we get the model, which is specified as a string from the
    # set seen below.
    if kwargs['model'] not in model_names:
        raise Exception('%s is not an allowed model!' %
                        kwargs['model'])
    else:
        # Here we use a bit of Python trickery to avoid a long list of
        # if/elif statements: we append the model abbreviation to the
        # build_model function name and then eval it:
        build_model_cmd = 'builder.build_model_%s()' % kwargs['model']
        eval(build_model_cmd)
        model = builder.model

    # Set the random seed:
    random_seed = int(kwargs['random_seed'])

    # Set the number of steps:
    nsteps = int(kwargs['nsteps'])

    # A sanity check to make sure everything worked:
    if None in [model, nbd_site, random_seed, nsteps, nbd_observable]:
        raise Exception('Something went wrong! One of the arguments to ' \
                        'do_fit was not initialized properly.')

    # We programmatically build up the base filename from the provided args:
    basename = '%s_%s_%s_%d_s%d' % (kwargs['model'], nbd_site, nbd_observable,
                                    nsteps, random_seed)

    # Call do_fit to get things started, feeding in the prior,
    # estimate_params, and randomized initial values from the Builder
    # object:
    numpy.random.seed(random_seed)

    mcmc = do_fit(model,
                  get_likelihood_function(nbd_site, nbd_observable),
                  prior=builder.prior,
                  estimate_params=builder.estimate_params,
                  basename=basename,
                  random_seed=random_seed,
                  nsteps=nsteps,
                  initial_values=builder.random_initial_values())

    # When sampling is completed, make the figures:
    generate_figures(mcmc, nbd_site, nbd_observable, basename=basename)

    print "Done."
