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
    ocb = one_cpt.Builder()
    lsb = lipo_sites.Builder()

if __name__ == '__main__':
    test_lipo_sites_builder()
    #test_one_cpt_lipo_sites_comparison()
