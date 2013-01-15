import bayessb.convergence
import pickle

chain_file_list = [
    'nbd_mcmc_100_steps_seed_0.pck',
    'nbd_mcmc_100_steps_seed_1.pck',
    'nbd_mcmc_100_steps_seed_2.pck',
    'nbd_mcmc_100_steps_seed_3.pck',
    'nbd_mcmc_100_steps_seed_4.pck']

chain_set = []
for file in chain_file_list:
    f = open(file, 'r')
    mcmc = pickle.load(f)

    mixed_start = mcmc.options.nsteps / 2
    chain_set.append(mcmc.positions[mixed_start:, :])

cc = bayessb.convergence.convergence_criterion(chain_set)
print cc
