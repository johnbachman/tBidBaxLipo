import sys
import pickle
import numpy as np
from matplotlib import pyplot as plt

from preprocess_data import data_to_fit, bg_time, lipo_concs_to_fit
from tbidbaxlipo.util import set_fig_params_for_publication, \
                             emcee_fit, format_axis
from pysb import *
from tbidbaxlipo.models.nbd.multiconf import Builder

bd = Builder()
bd.build_model_multiconf(num_confs=2, c0_scaling=1, normalized_data=True,
                        reversible=False)

def run_pt_mcmc(time, data, lipo_concs, filename=None):
    y = [data[-2]]
    bd.global_params = [bd['c1_scaling'], bd['c0_to_c1_k']]
    bd.local_params = []
    params = {}
    gf = emcee_fit.GlobalFit(bd, time, y, params, 'NBD')
    num_temps = 10
    num_walkers = 100
    burn_steps = 20
    sample_steps = 50
    sampler = emcee_fit.pt_mpi_sample(gf, num_temps, num_walkers, burn_steps,
                                      sample_steps)
    # Get rid of the pool so that the sampler can be pickled
    sampler.pool = None
    if filename is not None:
        with open(filename, 'w') as f:
            pickle.dump((gf, sampler), f)

    return (gf, sampler)

def plot_chain(gf, sampler, num_samples=500):
    # Check convergence
    plt.figure('Chains')
    plt.plot(sampler.lnprobability[0].T)
    plt.xlabel('MCMC Step')
    plt.ylabel('log(Posterior)')
    plt.title('Chain convergence')
    # Plot the data
    plt.figure()
    plt.plot(gf.time, gf.data[0])
    # Plot sample trajectories
    for i in range(num_samples):
        # Index 1 gives is the total number of steps (walkers * samples)
        s_ix = np.random.randint(sampler.flatchain.shape[1])
        # We're only interested in the lowest temp chain, chain 0;
        # Get the set of parameters at the randomly chosen index, s_ix
        gf.plot_func(sampler.flatchain[0, s_ix], alpha=0.1)
    # Triangle plot
    triangle.corner(sampler.flatchain[0])

if __name__ == '__main__':
    import triangle
    plt.ion()

    usage_msg =  "\nUsage:\n"
    usage_msg += "To run the fits and save chain to a pickle file:\n"
    usage_msg += "     python %s sample output_filename\n" % sys.argv[0]
    usage_msg += "To plot results from a pickled chain:\n"
    usage_msg += "     python %s plot input_filename\n" % sys.argv[0]

    # Check command-line args
    if len(sys.argv) < 3:
        print usage_msg
        sys.exit()

    # Sample
    pck_filename = sys.argv[2]
    if sys.argv[1] == 'sample':
        (gf, sampler) = run_pt_mcmc(bg_time, data_to_fit, lipo_concs_to_fit,
                                    filename=pck_filename)
    # Plot
    elif sys.argv[1] == 'plot':
        with open(pck_filename) as f:
            (gf, sampler) = pickle.load(f)
        plot_chain(gf, sampler)
    else:
        print usage_msg
        sys.exit()
