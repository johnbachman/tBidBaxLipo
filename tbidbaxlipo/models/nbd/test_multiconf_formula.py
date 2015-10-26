from tbidbaxlipo.models.nbd.multiconf import Builder
from pysb.integrate import Solver
from matplotlib import pyplot as plt
import numpy as np
from nose.tools import ok_

def test_2conf_irrev_formula():
    params_dict = {'c1_scaling': 5}
    bd = Builder(params_dict=params_dict)
    bd.build_model_multiconf(2, 1., reversible=False, normalized_data=True)
    t = np.linspace(0, 4000, 100)
    sol = Solver(bd.model, t)
    sol.run()
    nbd_sol = sol.yexpr['NBD']
    nbd_func = bd.obs_func(t)
    ok_(np.allclose(nbd_sol, nbd_func),
        'Integrated NBD does not match closed form NBD for 2conf model')

def test_2conf_rev_formula():
    params_dict = {'c1_scaling': 5}
    bd = Builder(params_dict=params_dict)
    bd.build_model_multiconf(2, 1., reversible=True, normalized_data=True)
    t = np.linspace(0, 4000, 100)
    sol = Solver(bd.model, t)
    sol.run()
    nbd_sol = sol.yexpr['NBD']
    nbd_func = bd.obs_func(t)
    ok_(np.allclose(nbd_sol, nbd_func),
        'Integrated NBD does not match closed form for 2conf reversible model')

def test_3conf_irrev_formula_k1_k2_different():
    params_dict = {'c1_scaling': 8,
                   'c2_scaling': 2,
                   'c0_to_c1_k': 0.05,
                   'c1_to_c2_k': 0.005
            }
    bd = Builder(params_dict=params_dict)
    bd.build_model_multiconf(3, 1., reversible=False, normalized_data=True)
    t = np.linspace(0, 4000, 100)
    sol = Solver(bd.model, t)
    sol.run()
    nbd_sol = sol.yexpr['NBD']
    nbd_func = bd.obs_func(t)
    ok_(np.allclose(nbd_sol, nbd_func),
        'Integrated NBD does not match closed form NBD for 3conf model')

def test_3conf_irrev_formula_k1_k2_same():
    params_dict = {'c1_scaling': 8,
                   'c2_scaling': 2,
                   'c0_to_c1_k': 0.05,
                   'c1_to_c2_k': 0.05
            }
    bd = Builder(params_dict=params_dict)
    bd.build_model_multiconf(3, 1., reversible=False, normalized_data=True)
    t = np.linspace(0, 4000, 100)
    sol = Solver(bd.model, t)
    sol.run()
    nbd_sol = sol.yexpr['NBD']
    nbd_func = bd.obs_func(t)
    ok_(np.allclose(nbd_sol, nbd_func),
        'Integrated NBD does not match closed form NBD for 3conf model')

if __name__ == '__main__':
    pass
