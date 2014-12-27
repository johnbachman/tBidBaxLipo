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

def fit_with_2conf(time, data, lipo_concs):
    y = [data[-2]]
    bd.global_params = [bd['c1_scaling'], bd['c0_to_c1_k']]
    bd.local_params = []
    params = {}
    gf = emcee_fit.GlobalFit(bd, time, y, params, 'NBD')
    num_temps = 10
    num_walkers = 100
    burn_steps = 10
    sample_steps = 50
    sampler = emcee_fit.pt_mpi_sample(gf, num_temps, num_walkers, burn_steps,
                                      sample_steps)
    return (gf, sampler)

def plot_sample_fits(gf, sampler, num_samples=500):
    # Plot the data
    plt.figure()
    plt.plot(gf.time, gf.data[0])
    # Plot sample trajectories
    for i in range(num_samples):
        s_ix = np.random.randint(sampler.flatchain.shape[0])
        gf.plot_func(sampler.flatchain[s_ix]

def plot_posterior(sampler):
    plt.figure()
    plt.plot(sampler.lnprobability.T)

if __name__ == '__main__':
    import triangle
    plt.ion()

    (gf, sampler) = fit_with_2conf(bg_time, data_to_fit, lipo_concs_to_fit)

    with open('exp_fit_mcmc.mcmc', 'w') as f:
        pickle.dump((gf, sampler), f)


