import tbidbaxlipo.nbd_mcmc_pysb
from nbd_mcmc_pysb import model_names, nbd_site_names
import nbd_analysis as nbd
import sys
import numpy as np
import bayessb
import pickle

if __name__ == '__main__':
    # Set the type of model to be built here
    from tbidbaxlipo.models.one_cpt import Builder

    # Prepare the data
    # ================

    tspan = nbd.time_other
    nbd_avgs, nbd_stds = nbd.calc_norm_avg_std()

    # The c62 data is on a slightly different timescale and has two extra
    # datapoints. For simplicity, we simply remove them.
    nbd_avgs[1] = nbd_avgs[1][0:-2]
    nbd_stds[1] = nbd_stds[1][0:-2]

    # Parse the args
    # ==============
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

    # We set the random_seed here because it affects our choice of initial
    # values
    np.random.seed(random_seed)

    # Initialize the MCMC arguments
    opts = bayessb.MCMCOpts()
    opts.model = model
    opts.tspan = nbd.time_other
    opts.estimate_params = builder.estimate_params
    opts.initial_values = builder.random_initial_values()
    opts.nsteps = nsteps

    opts.sigma_step = 0.9
    opts.sigma_max = 50
    opts.sigma_min = 0.01
    opts.accept_rate_target = 0.23
    opts.accept_window = 200
    opts.sigma_adj_interval = 200
    opts.anneal_length = nsteps / 10
    opts.use_hessian = True
    opts.hessian_scale = 1
    opts.hessian_period = opts.nsteps / 10 #10
    opts.seed = random_seed
    mcmc = tbidbaxlipo.nbd_mcmc_pysb.NBD_MCMC(opts, nbd_avgs, nbd_stds,
                                nbd_site, nbd_observable, builder)

    mcmc.do_fit()

    # Pickle it
    output_file = open('%s.mcmc' % mcmc.get_basename(kwargs['model']), 'w')
    pickle.dump(mcmc, output_file)
    output_file.close()

    # When sampling is completed, make the figures:
    #generate_figures(mcmc, nbd_site, nbd_observable, basename=basename)

    print "Done."
