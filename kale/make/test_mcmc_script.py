import numpy as np
from tbidbaxlipo.models.nbd.multiconf import Builder
from tbidbaxlipo.util import emcee_fit
from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
import tbidbaxlipo.plots.bid_bim_nbd_release.preprocess_data as ppd
import pickle
import sys

activator = 'Bid'
nbd_site = '126'
observable = 'NBD'
rep_num = 1

# Get the preprocessed data, along with the experimental error
name_pattern = '%s_%s_%s' % (activator, observable, nbd_site)
time_var = ppd.__dict__['time_%s_r%d' % (name_pattern, rep_num)]
data_var = ppd.__dict__['data_%s_r%d' % (name_pattern, rep_num)]
data_sigma_var = ppd.__dict__['data_sigma_%s' % name_pattern]

model_observable = ['NBD']

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

