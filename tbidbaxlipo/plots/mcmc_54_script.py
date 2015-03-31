import numpy as np
from tbidbaxlipo.models.nbd.multiconf import Builder
from tbidbaxlipo.util import emcee_fit
from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
import pickle
import sys

time_var = df[('Bid', 'NBD', '54', 1, 'TIME')].values
data_54 = np.zeros((1, 1, len(time_var)))
data_54[0, 0, :] = df[('Bid', 'NBD', '54', 1, 'VALUE')].values

model_observable = ['NBD']

# data_54_sigma = np.zeros((1, 1))
# data_54_sigma[0, 0] = calc_err_var_cubic(data_54[0, 0, :], last_n_pts=80)
data_54_sigma = np.array([[0.014037]])

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
gf = emcee_fit.GlobalFit(bd, time_var, data_54, data_54_sigma, ic_params,
                         model_observable)

random_seed = 1
np.random.seed(random_seed)

betas = 10 ** np.linspace(0, highest_temp, ntemps)
sampler = emcee_fit.pt_mpi_sample(gf, ntemps, nwalkers, nburnin, nsample,
                                  thin=thin, betas=betas)
sampler.pool = None

with open('test.mcmc', 'w') as f:
    pickle.dump((gf, sampler), f)

