import os
from tbidbaxlipo.pt import run_pt
from tbidbaxlipo.util.emcee_fit import global_fit_from_args
from nose.tools import ok_
import numpy as np

pt_args = {
    'data': {
        'module': 'tbidbaxlipo.plots.x140320_NBD_Bax_BimBH3_unlab_Bax_titration.' \
                  'preprocess_data',
        'data_var': 'data_to_fit',
        'data_sigma_var': 'data_sigma',
        'initial_condition_var': 'bax_concs_to_fit',
        'time_var': 'time',
    },
    'model': {
        'baxtranslocation': 1,
        'activation': 1,
        'nbd': 1,
        'bleach': 1,
    },
    'model_observable': ['NBD_bleach'],
    'global_initial_conditions': {
        'tBid_0': 0.,
        'c0_scaling': 1.,
        'Bax_NBD_0': 96.,
        'Vesicles_0': 1.9,
    },
    'local_initial_condition': 'Bax_0',
    'global_params': 'all',
    'local_params': [],
    'ntemps': 2,
    'highest_temp': -2,
    'nwalkers': 10,
    'nburnin': 1,
    'nsample': 1,
    'thin': 1,
}

def test_run_pt():
    sampler = run_pt.run(pt_args, 'test.mcmc', 1, 'test.pos', mpi=False)
    # Clean up
    if os.path.isfile('test.mcmc'):
        os.unlink('test.mcmc')
    if os.path.isfile('test.pos'):
        os.unlink('test.pos')

if __name__ == '__main__':
    test_run_pt()

