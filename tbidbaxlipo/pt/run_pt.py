import yaml
import pickle
import sys
import os
from tbidbaxlipo.util import emcee_fit
import numpy as np

def run(args, mcmc_filename, random_seed, pos_filename, mpi=True):
    """Takes dict of args and initializes PT run.

    The args can come from a YAML file or be passed in by a script.
    """

    gf = emcee_fit.global_fit_from_args(args)

    # Seed the random number generator
    np.random.seed(random_seed)

    # Make the beta ladder
    ntemps = args['ntemps']
    betas = 10 ** np.linspace(0, args['highest_temp'], ntemps)

    ### RUN the sampler!
    # If the file does exist, load the position as the start position
    if os.path.isfile(pos_filename):
        # Open the position/random state file
        with open(pos_filename) as f:
            (pos, rs) = pickle.load(f)
        # Set the global random state from the saved value
        np.random.set_state(rs)
        print("Continuing run with position and random state found in "
              "position file %s" % pos_filename)
    # If the file does not exist
    else:
        print("Position file %s does not exist, it will be created." %
              pos_filename)
        pos = None

    # Call the appropriate sampling func depending on whether we're using
    # MPI
    if mpi:
        sample_func = emcee_fit.pt_mpi_sample
    else:
        sample_func = emcee_fit.pt_sample
    sampler = sample_func(gf, ntemps, args['nwalkers'], args['nburnin'],
                          args['nsample'], thin=args['thin'], betas=betas,
                          pos_filename=pos_filename, pos=pos)

    # After sampling, get rid of the pool so we can pickle the sampler
    sampler.pool = None

    with open(mcmc_filename, 'w') as f:
        pickle.dump((gf, sampler), f)

if __name__ == '__main__':
    # Parameters for the job are specified in a YAML file, which should be
    # provided as an argument.
    # Check that a path to a YAML file was provided
    if len(sys.argv) < 4:
        print("Usage: %s yaml_file random_seed pos_filename [opt: nompi]" % \
               __file__)
        sys.exit()

    # Load the args from the YAML file
    with open(sys.argv[1]) as yaml_file:
        args = yaml.load(yaml_file)

    # The basename of the pickle file is based on the name of the .yaml file
    basedir = os.path.dirname(sys.argv[1])
    basename = os.path.basename(sys.argv[1])
    basename = basename.split('.')[0]
    mcmc_filename = os.path.join(basedir, '%s.mcmc' % basename)
    # Get the random seed and the filename to pickle the current position to
    random_seed = int(sys.argv[2])
    pos_filename = sys.argv[3]
    use_mpi = False if 'nompi' in sys.argv else True
    run(args, mcmc_filename, random_seed, pos_filename, mpi=use_mpi)
