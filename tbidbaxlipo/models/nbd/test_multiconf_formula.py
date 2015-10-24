from tbidbaxlipo.models.nbd.multiconf import Builder
from pysb.integrate import Solver
from matplotlib import pyplot as plt
import numpy as np
from nose.tools import ok_

def test_2conf_formula():
    params_dict = {'c1_scaling': 5}
    bd = Builder(params_dict=params_dict)
    bd.build_model_multiconf(2, 1., reversible=False, normalized_data=True)
    t = np.linspace(0, 4000, 100)
    sol = Solver(bd.model, t)
    sol.run()
    nbd_sol = sol.yexpr['NBD']
    nbd_func = bd.nbd_func(t)
    ok_(np.allclose(nbd_sol, nbd_func),
        'Integrated NBD does not match closed form NBD for 2conf model')

if __name__ == '__main__':
    test_2conf_formula()

