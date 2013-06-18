from tbidbaxlipo.mcmc.pore_mcmc import PoreMCMC
from tbidbaxlipo.mcmc.pore_mcmc import model_names
import sys
import numpy as np
import bayessb
import pickle

if __name__ == '__main__':
    # Set the type of model to be built here
    from tbidbaxlipo.models.lipo_sites import Builder

    # Allowable compartmentalization type names:
    cpt_types = ['one_cpt', 'lipo_sites']

    # Prepare the data
    # ================
    data = pickle.load(open('../data/130614_norm_timecourses_df.pck'))
    dataset_name = '130614'

    # Parse the args
    # ==============
    # Keyword args are set at the command line as e.g., key=val
    # and subsequently split at the equals sign
    kwargs = dict([arg.split('=') for arg in sys.argv[1:]])

    # Arguments that can take a list of values (e.g., the nbd sites or
    # nbd observables) are expected to be separated by this character:
    SEP_CHAR = '-'

    # We set these all to None so later on we can make sure they were
    # properly initialized.
    random_seed = None
    model = None

    print "Keyword arguments: "
    print kwargs

    # Before we begin, we make sure we have all the keyword arguments that
    # we are going to need.
    if 'random_seed' not in kwargs or \
       'cpt_type' not in kwargs or \
       'nsteps' not in kwargs or \
       'model' not in kwargs:
        raise Exception('One or more needed arguments was not specified! ' \
                'Arguments must include random_seed, model, cpt_type, ' \
                'and nsteps.')

    # First we get the compartmentalization type:
    if kwargs['cpt_type'] not in cpt_types:
        raise Exception('Allowable values for cpt_type are: one_cpt, '
                        'lipo_sites')

    # Get the appropriate builder
    if cpt_type == 'lipo_sites':
        builder = lipo_sites.Builder()
    elif cpt_type == 'one_cpt':
        builder = one_cpt.Builder()

    # Now we get the model, which is specified as a string from the
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

    # Set the initial temperature
    if 'T_init' in kwargs:
        T_init = float(kwargs['T_init'])
    else:
        T_init = 1

    # A sanity check to make sure everything worked:
    if None in [model, random_seed, nsteps]:
        raise Exception('Something went wrong! One of the arguments to ' \
                        'do_fit was not initialized properly.')

    # We set the random_seed here because it affects our choice of initial
    # values
    np.random.seed(random_seed)

    # Initialize the MCMC arguments
    opts = bayessb.MCMCOpts()
    opts.model = model
    opts.tspan = np.linspace(0, 6000, 196) # This should not be used
    opts.estimate_params = builder.estimate_params
    opts.initial_values = builder.random_initial_values()
    opts.nsteps = nsteps
    opts.norm_step_size = 0.01
    opts.sigma_step = 0
    #opts.sigma_max = 50
    #opts.sigma_min = 0.01
    #opts.accept_rate_target = 0.23
    #opts.accept_window = 100
    #opts.sigma_adj_interval = 200
    opts.anneal_length = 0
    opts.use_hessian = False #True # TRUE
    opts.hessian_scale = 1
    opts.hessian_period = opts.nsteps / 20 #10
    opts.seed = random_seed
    opts.T_init = 1 # T_init
    mcmc = PoreMCMC(opts, data, dataset_name, builder)

    mcmc.do_fit()

    # Pickle it
    output_file = open('%s.mcmc' % mcmc.get_basename(), 'w')
    pickle.dump(mcmc, output_file)
    output_file.close()

    # When sampling is completed, make the figures:
    #mcmc.generate_figures(report_name=mcmc.get_basename(),
    #                      mixed_start=0)

    print "Done."
