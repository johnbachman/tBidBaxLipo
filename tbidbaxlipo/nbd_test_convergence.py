import bayessb.convergence
import pickle
import numpy as np
import matplotlib.pyplot as plt

chain_file_list = [
    'nbd_mcmc_c3_1000_steps_seed_1.pck',
    'nbd_mcmc_c3_1000_steps_seed_2.pck',
    'nbd_mcmc_c3_1000_steps_seed_3.pck']

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

plt.figure()
plt.plot(chain_set[0][:,0])
plt.plot(chain_set[1][:,0])
plt.plot(chain_set[2][:,0])
plt.show()

plt.figure()
plt.plot(chain_set[0][:,1])
plt.plot(chain_set[1][:,1])
plt.plot(chain_set[2][:,1])
plt.show()

