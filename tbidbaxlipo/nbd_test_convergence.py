import bayessb.convergence
import pickle
import numpy as np
import matplotlib.pyplot as plt
import glob

chain_file_list = glob.glob('nbd_mcmc_parallel_*steps_seed*.pck')

mcmc_set = []
#total_steps = 0

#mixed_start = 0

for file in chain_file_list:
    f = open(file, 'r')
    mcmc = pickle.load(f)
    mcmc_set.append(mcmc)

    #total_steps += len(mcmc.positions[mixed_start:, :])
    #chain_set.append(mcmc.positions[mixed_start:, :])

cc = bayessb.convergence.convergence_criterion(mcmc_set, mask=2000, thin=5)
print cc

#print chain_set[0].shape[1]
#pooled_chain = np.empty((total_steps, chain_set[0].shape[1]))

# Combine the steps together
# TODO: Implement this as a function for MCMC
#start = 0
#for i in range(0, len(chain_set)):
#    num_steps = chain_set[i].shape[0]
#    pooled_chain[start:(start + num_steps), :] = chain_set[i]
#    start += num_steps

plt.ion()

num_params = len(mcmc.options.estimate_params)
for i in range(num_params):
    plt.figure()
    #for mcmc in mcmc_set:
    #    plt.plot(mcmc.positions[:,i])
    for j, mcmc in enumerate(mcmc_set):
        mixed_positions = mcmc.positions[mixed_start:,:]
        mixed_accepts = mixed_positions[mcmc.accepts[mixed_start:]]
        thinned_accepts = mixed_accepts[::thin]
        accept_steps = steps[mcmc.accepts[mixed_start:]]
        thinned_accept_steps = accept_steps[::thin]
        plt.plot(thinned_accept_steps, thinned_accepts[:,i], label=chain_file_l
    plt.legend(loc='lower left', prop={'size':7})
    plt.show()

