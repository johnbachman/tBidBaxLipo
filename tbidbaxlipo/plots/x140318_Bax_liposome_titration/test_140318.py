import os
from tbidbaxlipo.pt import run_pt
from tbidbaxlipo.util.emcee_fit import global_fit_from_args
from nose.tools import ok_
import numpy as np

pt_args = {
    'data': {
        'module': 'tbidbaxlipo.plots.x140318_Bax_liposome_titration.' \
                  'preprocess_data',
        'data_var': 'data_to_fit',
        'data_sigma_var': 'data_sigma',
        'initial_condition_var': 'lipo_concs_to_fit',
        'time_var': 'bg_time',
    },
    'model': {
        'baxtranslocation': 1,
        'activation': 1,
        'nbd': 1,
    },
    'model_observable': ['NBD'],
    'global_initial_conditions': {
        'tBid_0': 0,
        'c0_scaling': 1.0,
        'Bax_NBD_0': 185.0,
    },
    'local_initial_condition': 'Vesicles_0',
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
    globals().update(locals())

def test_ysim_nans():
    """This test is here because I accidentally set Bax_0, rather than
    Bax_NBD_0, as my initial condition. This led to division by 0, and caused
    the solver results to be filled with NaNs. This checks to make sure that
    the model can be run with these arguments without yielding NaNs."""
    gf = global_fit_from_args(pt_args)
    gf.solver.run()
    ok_(not np.any(np.isnan(gf.solver.yexpr['NBD'])))

if __name__ == '__main__':
    test_ysim_nans()
    test_run_pt()

