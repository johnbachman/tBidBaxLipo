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

        # COMPARTMENTS
        solution = self.compartment('solution', dimension=3, parent=None)
        self.compartment('ves', dimension=2, parent=solution)

        # INITIAL CONDITIONS
        self.parameter('Vesicles_0', 5, prior=Normal(1., 2.))
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
        print("peri_lipo_sites: translocate_Bax()")

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

    def translocate_convert_Bax(self):
        print("peri_lipo_sites: translocate_convert_Bax()")

        Bax_transloc_kf = self.parameter('Bax_transloc_kf', 1e-2)
        Bax_transloc_kr = self.parameter('Bax_transloc_kr', 1e-1)
        Bax_transloc_kcat = self.parameter('Bax_transloc_kcat', 1e-1)
        Bax_dissoc_kr = self.parameter('Bax_dissoc_kr', 1e-1)

        Bax = self['Bax']
        Vesicles = self['Vesicles']
        solution = self['solution']
        ves = self['ves']

        self.rule('Bax_translocates_sol_to_ves',
             Bax(loc='c', lipo=None) ** solution +
             Vesicles(bax=None) ** solution <>
             Bax(loc='c', lipo=1) ** ves % Vesicles(bax=1) ** solution,
             Bax_transloc_kf, Bax_transloc_kr)
        self.rule('Bax_associates_with_mito',
             Bax(loc='c', lipo=1) ** ves % Vesicles(bax=1) ** solution >>
             Bax(loc='m', lipo=None) ** ves + Vesicles(bax=None) ** solution,
             Bax_transloc_kcat)
        self.rule('Bax_dissociates_ves_to_sol',
             Bax(loc='m', lipo=None) ** ves >>
             Bax(loc='c', lipo=None) ** solution,
             Bax_dissoc_kr)

    def basal_Bax_activation(self, reversible=False):
        print "peri_lipo_sites: basal_Bax_activation, reversible=%s" % reversible

        Bax = self['Bax']
        Vesicles = self['Vesicles']
        solution = self['solution']
        ves = self['ves']

        basal_Bax_kf = self.parameter('basal_Bax_kf', 2e-3)
        self.rule('basal_Bax_activation',
                  Bax(bh3=None, loc='m', lipo=1) ** ves %
                  Vesicles(bax=1) ** solution >>
                  Bax(bh3=None, loc='i', lipo=None) ** ves +
                  Vesicles(bax=None) ** solution,
                  basal_Bax_kf)

    def secondary_translocation(self):
        print("peri_lipo_sites: secondary_translocation()")

        Bax_sec_transloc_kf = self.parameter('Bax_sec_transloc_kf', 1e-2)
        Bax_sec_transloc_kr = self.parameter('Bax_sec_transloc_kr', 1e-1)

        Bax = self['Bax']
        Vesicles = self['Vesicles']
        solution = self['solution']
        ves = self['ves']

        self.rule('Bax_sec_translocates_sol_to_ves',
             Bax(loc='m', lipo=1) ** ves +
             Bax(loc='c', lipo=None) ** solution >>
             Bax(loc='m', lipo=1) ** ves +
             Bax(loc='m', lipo=None) ** ves,
             Bax_sec_transloc_kf,
             Bax_sec_transloc_kr)

