from pysb import *
import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.util.fitting import fit, mse
from tbidbaxlipo.util import color_iter
from tbidbaxlipo.models import one_cpt
from matplotlib.font_manager import FontProperties
from pysb.integrate import odesolve, Solver

Solver._use_inline = True

class Builder(one_cpt.Builder):

    def translocate_Bax(self):
        print("lipo_sites: translocate_Bax()")

        Bax_transloc_kf = self.parameter('Bax_transloc_kf', 1e-2)
        Bax_transloc_kr = self.parameter('Bax_transloc_kr', 1e-1)

        Bax = self['Bax']
        Vesicles = self['Vesicles']
        solution = self['solution']
        ves = self['ves']

        self.rule('Bax_translocates_sol_to_ves',
             Bax(loc='c', lipo=None) ** solution +
             Vesicles(bax=None) ** solution <>
             Bax(loc='m', lipo=1) ** ves % Vesicles(bax=1) ** solution,
             Bax_transloc_kf, Bax_transloc_kr)



