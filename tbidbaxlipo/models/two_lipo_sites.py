from pysb import *
import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.util.fitting import fit, mse
from tbidbaxlipo.util import color_iter
from tbidbaxlipo.models import one_cpt, core
from matplotlib.font_manager import FontProperties
from pysb.integrate import odesolve, Solver
from bayessb.priors import Normal

Solver._use_inline = True

class Builder(one_cpt.Builder):

    def __init__(self, params_dict=None, nbd_sites=None):
        """Differs from one_cpt.__init__ only in that it allows the
        concentration of vesicles, Vesicles_0, to be estimated."""
        # Sets self.model = Model(), and self.param_dict
        core.Builder.__init__(self, params_dict=params_dict)

        self.declare_monomers()
        Lipo_sat = self.monomer('Lipo_sat', ['bax'])

        # COMPARTMENTS
        solution = self.compartment('solution', dimension=3, parent=None)
        self.compartment('ves', dimension=2, parent=solution)

        # INITIAL CONDITIONS
        Vesicles_0 = self.parameter('Vesicles_0', 5, estimate=False)
        Lipo_sat_factor = self.parameter('Lipo_sat_factor', 10, estimate=True)
        #Lipo_nonsat_factor = self.parameter('Lipo_nonsat_factor', 0.5,
        #                                  estimate=True)
        Lipo_sat_0 = self.expression('Lipo_sat_0', Vesicles_0 * Lipo_sat_factor)
        #Lipo_nonsat_0 = self.expression('Lipo_nonsat_0',
        #                               Vesicles_0 * Lipo_nonsat_factor)
        self.parameter('tBid_0', 20, estimate=False)
        self.parameter('Bax_0', 100, estimate=False)

        tBid = self['tBid']
        Bax = self['Bax']
        Vesicles = self['Vesicles']

        self.initial(tBid(loc='c', bh3=None) ** solution, self['tBid_0'])
        self.initial(Bax(loc='c', bh3=None, a6=None, lipo=None,
                     c3='s', c62='s', c120='s', c122='s', c126='s', c184='s')
                     ** solution, self['Bax_0'])
        self.initial(Vesicles(bax=None) ** solution, self['Vesicles_0'])
        self.initial(Lipo_sat(bax=None) ** solution, Lipo_sat_0)

        # OBSERVABLES
        self.observable('ctBid', tBid(loc='c'))
        self.observable('mtBid', tBid(loc='m'))
        self.observable('cBax', Bax(loc='c'))
        self.observable('mBax', Bax(loc='m'))
        self.observable('iBax', Bax(loc='i'))
        self.observable('aBax', Bax(loc='a'))
        self.observable('tBidBax', tBid(bh3=1) % Bax(bh3=1))
        self.observable('Bax2', Bax(bh3=1) % Bax(bh3=1))
        self.observable('Baxbh3', Bax(bh3=1))
        self.observable('Bax4',
             MatchOnce(Bax(loc='i', bh3=1, a6=3) % Bax(loc='i', bh3=1, a6=4) % 
                       Bax(loc='i', bh3=2, a6=3) % Bax(loc='i', bh3=2, a6=4)))

        # SCALING PARAMETERS
        if nbd_sites is not None:
            self.declare_nbd_scaling_parameters(nbd_sites)

    def translocate_Bax(self):
        print("two_lipo_sites: translocate_Bax()")

        Bax_transloc_sat_kf = self.parameter('Bax_transloc_sat_kf', 1e-2)
        Bax_transloc_sat_kr = self.parameter('Bax_transloc_sat_kr', 1e-1)
        Bax_transloc_nonsat_kf = self.parameter('Bax_transloc_nonsat_kf', 1e-2)
        Bax_transloc_nonsat_kr = self.parameter('Bax_transloc_nonsat_kr', 1e-1)

        Bax = self['Bax']
        Lipo_sat = self['Lipo_sat']
        solution = self['solution']
        Vesicles = self['Vesicles']
        ves = self['ves']

        self.rule('Bax_translocates_to_sat_sites',
             Bax(loc='c', lipo=None) ** solution +
             Lipo_sat(bax=None) ** solution <>
             Bax(loc='m', lipo=1) ** ves % Lipo_sat(bax=1) ** solution,
             Bax_transloc_sat_kf, Bax_transloc_sat_kr)

        self.rule('Bax_translocates_to_nonsat_sites',
             Bax(loc='c', lipo=None) ** solution +
             Vesicles() ** solution <>
             Bax(loc='m', lipo=None) ** ves + Vesicles() ** solution,
             Bax_transloc_nonsat_kf, Bax_transloc_nonsat_kr)

    def basal_Bax_activation(self, reversible=False):
        print("two_lipo_sites: basal_Bax_activation, reversible=%s" % reversible)
        Bax = self['Bax']
        basal_Bax_sat_kf = self.parameter('basal_Bax_sat_kf', 2e-3)
        basal_Bax_nonsat_kf = self.parameter('basal_Bax_nonsat_kf', 2e-3)

        self.rule('basal_Bax_activation_sat',
                Bax(bh3=None, loc='m', lipo=1) >>
                Bax(bh3=None, loc='i', lipo=1),
                basal_Bax_sat_kf)
        self.rule('basal_Bax_activation_nonsat',
                Bax(bh3=None, loc='m', lipo=None) >>
                Bax(bh3=None, loc='i', lipo=None),
                basal_Bax_nonsat_kf)

        if reversible:
            basal_Bax_sat_kr = self.parameter('basal_Bax_sat_kr', 1e-2)
            basal_Bax_nonsat_kr = self.parameter('basal_Bax_nonsat_kr', 1e-2)

            self.rule('basal_Bax_activation_sat',
                    Bax(bh3=None, loc='i', lipo=1) >>
                    Bax(bh3=None, loc='m', lipo=1),
                    basal_Bax_sat_kr)
            self.rule('basal_Bax_activation_nonsat',
                    Bax(bh3=None, loc='i', lipo=None) >>
                    Bax(bh3=None, loc='m', lipo=None),
                    basal_Bax_nonsat_kr)

    """
    def Bax_reverses(self):
        print('two_lipo_sites: Bax_reverses')

        Bax = self['Bax']
        sol = self['solution']
        ves = self['ves']
        krev = self.parameter('iBax_reverse_k', 1e-3)
        self.rule('iBax_nonsat_reverses',
                Bax(loc='i'
    """


