import yaml
import pickle
import sys
from tbidbaxlipo.models.one_cpt import Builder
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
    ic_var = data_module.__dict__[data_args['initial_condition_var']]
    time_var = data_module.__dict__[data_args['time_var']]

    ### MODEL
    # Call the appropriate model-building macro
    bd = Builder()
    model_macro = getattr(bd, args['model_macro'])
    model_macro()
    # Set the initial conditions
    for ic_name, ic_value in args['global_initial_conditions'].iteritems():
        bd.model.parameters[ic_name].value = ic_value

    ### PARAMETERS TO FIT
    bd.global_params = [bd.model.parameters[p_name]
                        for p_name in args['global_params']]
    bd.local_params = [bd.model.parameters[p_name]
                       for p_name in args['local_params']]
    params = {args['local_initial_condition']: ic_var}

    # Create the global fit instance
    gf = emcee_fit.GlobalFit(bd, time_var, data_var, params,
                             args['model_observable'])

    # Seed the random number generator
    np.random.seed(int(sys.argv[2]))

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

