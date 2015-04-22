import numpy as np
from tbidbaxlipo.models.nbd.multiconf import Builder
from tbidbaxlipo.util import emcee_fit
from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
import pickle
import sys

time_var = df[('Bid', 'NBD', '126', 1, 'TIME')].values
data_var = np.zeros((1, 1, len(time_var)))
data_var[0, 0, :] = df[('Bid', 'NBD', '126', 1, 'VALUE')].values

model_observable = ['NBD']

# data_sigma_var = np.zeros((1, 1))
# data_sigma_var[0, 0] = calc_err_var_cubic(data_var[0, 0, :], last_n_pts=80)
#54C: data_sigma_var = np.array([[0.014037]])
# 126C:
data_sigma_var = np.array([[0.021085]])

highest_temp = -2
nburnin = 10
nsample = 20
ntemps = 3
nwalkers = 200
thin = 1
ic_params = None

bd = Builder()
bd.build_model_multiconf(2, 1, normalized_data=True)

# Set the variables in the bd
bd.global_params = bd.estimate_params
bd.local_params = []

# Create the globalfit instance
gf = emcee_fit.GlobalFit(bd, time_var, data_var, data_sigma_var, ic_params,
                         model_observable)

random_seed = 1
np.random.seed(random_seed)

betas = 10 ** np.linspace(0, highest_temp, ntemps)
sampler = emcee_fit.pt_mpi_sample(gf, ntemps, nwalkers, nburnin, nsample,
                                  thin=thin, betas=betas)
sampler.pool = None

with open('test.mcmc', 'w') as f:
    pickle.dump((gf, sampler), f)

