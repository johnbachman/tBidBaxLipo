import yaml
import pickle
import sys
from tbidbaxlipo.models import one_cpt
from tbidbaxlipo.models.nbd import multiconf
from tbidbaxlipo.util import emcee_fit
import numpy as np

if __name__ == '__main__':
    # Parameters for the job are specified in a YAML file, which should be
    # provided as an argument.
    # Check that a path to a YAML file was provided
    if len(sys.argv) < 3:
        print("Usage: %s yaml_file random_seed" % __file__)
        sys.exit()

    # Load the args from the YAML file
    with open(sys.argv[1]) as yaml_file:
        args = yaml.load(yaml_file)

    ### DATA
    # Import the module containing the data
    data_args = args['data']
    __import__(data_args['module'])
    data_module = sys.modules[data_args['module']]
    # Get the relevant variables from the module containing the data
    data_var = data_module.__dict__[data_args['data_var']]
    data_sigma_var = data_module.__dict__[data_args['data_sigma_var']]
    time_var = data_module.__dict__[data_args['time_var']]

    # Get the name of the variable containing the initial conditions vector,
    # which may not exist
    ic_var_name = data_args['initial_condition_var']
    if ic_var_name is None:
        ic_var = None
    else:
        ic_var = data_module.__dict__[ic_var_name]

    ### MODEL
    # Call the appropriate model-building macro
    # If there is a multiconf attribute, that trumps any other attribute
    # and determines that this is a multiconf model
    if 'multiconf' in args['model']:
        bd = multiconf.Builder()
        num_confs = args['model']['multiconf']
        norm_data = args['model']['normalized_nbd_data'] = True
        bd.build_model_multiconf(num_confs, 1, normalized_data=norm_data)
    else:
        bd = one_cpt.Builder()
        bd.build_model_from_dict(args['model'])

    # Set the initial conditions
    for ic_name, ic_value in args['global_initial_conditions'].iteritems():
        bd.model.parameters[ic_name].value = ic_value

    ### PARAMETERS TO FIT
    if args['global_params'] == 'all':
        bd.global_params = bd.estimate_params
        bd.local_params = []
    else:
        bd.global_params = [bd.model.parameters[p_name]
                            for p_name in args['global_params']]
        bd.local_params = [bd.model.parameters[p_name]
                           for p_name in args['local_params']]

    local_ic_name = args['local_initial_condition']
    if ic_var is None or local_ic_name is None:
        params = None
    else:
        params = {local_ic_name: ic_var}

    # Create the global fit instance
    gf = emcee_fit.GlobalFit(bd, time_var, data_var, data_sigma_var, params,
                             args['model_observable'])

    # Seed the random number generator
    random_seed = int(sys.argv[2])
    np.random.seed(random_seed)

    # Make the beta ladder
    ntemps = args['ntemps']
    betas = 10 ** np.linspace(0, args['highest_temp'], ntemps)

    ### RUN the sampler!
    sampler = emcee_fit.pt_mpi_sample(gf,
                                      ntemps,
                                      args['nwalkers'],
                                      args['nburnin'],
                                      args['nsample'],
                                      thin=args['thin'],
                                      betas=betas)

    # After sampling, get rid of the pool so we can pickle the sampler
    sampler.pool = None
    # The basename of the pickle file is based on the name of the .yaml file
    basename = sys.argv[1].split('.')[0]
    with open('%s_%d.mcmc' % (basename, random_seed), 'w') as f:
        pickle.dump((gf, sampler), f)

    # Done
    sys.exit()
