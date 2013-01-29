import bayessb.convergence
import pickle
import numpy as np
import matplotlib.pyplot as plt
import glob

#chain_file_list = [
#    'nbd_mcmc_linear_2000_steps_seed_0.pck',
#    'nbd_mcmc_linear_2000_steps_seed_1.pck',
#    'nbd_mcmc_linear_2000_steps_seed_2.pck',
#    'nbd_mcmc_linear_2000_steps_seed_4.pck',
#    'nbd_mcmc_linear_2000_steps_seed_5.pck',
#    'nbd_mcmc_linear_2000_steps_seed_7.pck',
#    'nbd_mcmc_linear_2000_steps_seed_8.pck',
#    'nbd_mcmc_linear_2000_steps_seed_9.pck',
#    'nbd_mcmc_linear_2000_steps_seed_10.pck',
#    'nbd_mcmc_linear_2000_steps_seed_12.pck',
#    'nbd_mcmc_linear_2000_steps_seed_13.pck',
#    'nbd_mcmc_linear_2000_steps_seed_15.pck',
#    'nbd_mcmc_linear_2000_steps_seed_16.pck',
#    'nbd_mcmc_linear_2000_steps_seed_18.pck']

chain_file_list = glob.glob('nbd_mcmc_parallel_*steps_seed*.pck')

chain_set = []
total_steps = 0

for file in chain_file_list:
    f = open(file, 'r')
    mcmc = pickle.load(f)
    mixed_start = mcmc.options.nsteps / 2
    total_steps += len(mcmc.positions[mixed_start:, :])
    chain_set.append(mcmc.positions[mixed_start:, :])

cc = bayessb.convergence.convergence_criterion(chain_set)
print cc

#print chain_set[0].shape[1]
pooled_chain = np.empty((total_steps, chain_set[0].shape[1]))

# Combine the steps together
# TODO: Implement this as a function for MCMC
start = 0
for i in range(0, len(chain_set)):
    num_steps = chain_set[i].shape[0]
    pooled_chain[start:(start + num_steps), :] = chain_set[i]
    start += num_steps

plt.ion()

num_params = len(mcmc.options.estimate_params)
for i in range(num_params):
    plt.figure()
    for chain in chain_set:
        plt.plot(chain[:,i])
    plt.show()

