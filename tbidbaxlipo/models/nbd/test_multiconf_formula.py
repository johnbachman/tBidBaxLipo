from tbidbaxlipo.models.nbd.multiconf import Builder
from pysb.integrate import Solver
from matplotlib import pyplot as plt
import numpy as np
from nose.tools import ok_, raises

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

def test_4conf_irrev_formula_k1_k2_k3_different():
    params_dict = {'c1_scaling': 8,
                   'c2_scaling': 2,
                   'c3_scaling': 4,
                   'c0_to_c1_k': 0.05,
                   'c1_to_c2_k': 0.005,
                   'c2_to_c3_k': 0.001
            }
    bd = Builder(params_dict=params_dict)
    bd.build_model_multiconf(4, 1., reversible=False, normalized_data=True)
    t = np.linspace(0, 4000, 100)
    sol = Solver(bd.model, t)
    sol.run()
    nbd_sol = sol.yexpr['NBD']
    nbd_func = bd.obs_func(t)
    ok_(np.allclose(nbd_sol, nbd_func),
        'Integrated NBD does not match closed form NBD for 4conf model')

@raises(ValueError)
def test_4conf_irrev_formula_k1_k2_same():
    params_dict = {'c1_scaling': 8,
                   'c2_scaling': 2,
                   'c3_scaling': 4,
                   'c0_to_c1_k': 0.05,
                   'c1_to_c2_k': 0.05,
                   'c2_to_c3_k': 0.001
            }
    bd = Builder(params_dict=params_dict)
    bd.build_model_multiconf(4, 1., reversible=False, normalized_data=True)
    t = np.linspace(0, 4000, 100)
    bd.obs_func(t)

def test_5conf_irrev_formula_k1_k2_k3_k4_different():
    params_dict = {'c1_scaling': 8,
                   'c2_scaling': 2,
                   'c3_scaling': 6,
                   'c4_scaling': 1.5,
                   'c0_to_c1_k': 0.05,
                   'c1_to_c2_k': 0.005,
                   'c2_to_c3_k': 0.001,
                   'c3_to_c4_k': 0.002,
            }
    bd = Builder(params_dict=params_dict)
    bd.build_model_multiconf(5, 1., reversible=False, normalized_data=True)
    t = np.linspace(0, 4000, 100)
    sol = Solver(bd.model, t)
    sol.run()
    nbd_sol = sol.yexpr['NBD']
    nbd_func = bd.obs_func(t)
    plt.ion()
    plt.figure()
    plt.plot(t, nbd_sol)
    plt.plot(t, nbd_func)
    ok_(np.allclose(nbd_sol, nbd_func),
        'Integrated NBD does not match closed form NBD for 5conf model')

if __name__ == '__main__':
    test_5conf_irrev_formula_k1_k2_k3_k4_different()
