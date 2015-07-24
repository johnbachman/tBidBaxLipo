from nose.tools import ok_
from tbidbaxlipo.models import one_cpt, lipo_sites
from pysb.integrate import Solver
import numpy as np
from matplotlib import pyplot as plt

def test_lipo_sites_builder():
    """lipo_sites_builder builds and runs without error."""
    lsb = lipo_sites.Builder()
    lsb.translocate_Bax()
    t = np.linspace(0, 100)
    sol = Solver(lsb.model, t)
    sol.run()

def test_one_cpt_lipo_sites_comparison():
    """Simulation of one_cpt and the lipo_sites model with Bax_0 numbers
    of sites per liposome should give the same result."""
    # A simulation of
    params_dict = {'Bax_0': 10., 'sites_per_liposome': 1000.}
    ocb = one_cpt.Builder(params_dict=params_dict)
    lsb = lipo_sites.Builder(params_dict=params_dict)
    t = np.linspace(0, 100, 1000)
    plt.ion()
    y_list = []
    for builder in (ocb, lsb):
        builder.translocate_Bax()
        sol = Solver(builder.model, t)
        sol.run()
        y_list.append(sol.yobs['mBax'])
    # The two arrays have a maximal absolute difference of 1.481e-3
    ok_(np.allclose(y_list[0], y_list[1], atol=1.5e-3),
        'one_cpt and lipo_sites simulations should be approximately equal '
        'under these conditions.')

if __name__ == '__main__':
    #test_lipo_sites_builder()
    test_one_cpt_lipo_sites_comparison()
