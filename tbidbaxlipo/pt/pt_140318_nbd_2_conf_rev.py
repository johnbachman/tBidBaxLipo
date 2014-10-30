import sys
import pickle
from tbidbaxlipo.models import one_cpt
from tbidbaxlipo.util import emcee_fit
from tbidbaxlipo.plots.layout_140318 import data_to_fit, lipo_concs_to_fit, \
                                            bg_time
import numpy as np

if len(sys.argv) < 2:
    print "Please include the random seed as an argument."
    sys.exit()

random_seed = int(sys.argv[1])
np.random.seed(random_seed)

params_dict = {'Bax_0': 185.,
               'c1_scaling': 4.85,
               'basal_Bax_kf': 1.73e-3,
               'Bax_transloc_kr': 1.4835e-1,
               'Bax_transloc_kf': 1e-2,
               'iBax_reverse_k': 1e-3}
bd = one_cpt.Builder(params_dict=params_dict)
bd.build_model_nbd_2_conf_rev()
bd.global_params = (bd['c1_scaling'],
                    bd['Bax_transloc_kf'],
                    bd['Bax_transloc_kr'],
                    bd['basal_Bax_kf'],
                    bd['iBax_reverse_k'])
bd.local_params = []
params = {'Vesicles_0': lipo_concs_to_fit}
gf = emcee_fit.GlobalFit(bd, bg_time, data_to_fit, params, 'NBD')
ntemps = 25
betas = 10 ** np.linspace(0, -6, ntemps)
sampler = emcee_fit.pt_mpi_sample(gf, ntemps, 300, 10, 250)

# Get rid of the pool so we can pickle the sampler
sampler.pool = None
with open('pt_140318_nbd_2_conf_rev_%d.pck' % random_seed, 'w') as f:
    pickle.dump((gf, sampler), f)

sys.exit()
