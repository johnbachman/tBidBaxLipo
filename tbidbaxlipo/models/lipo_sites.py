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
        self.parameter('Vesicles_0', 20, prior=Normal(1., 2.))
        self.parameter('tBid_0', 20, prior=None)
        self.parameter('Bax_0', 100, prior=None)

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
        print("lipo_sites: translocate_Bax()")

        Bax_transloc_kf = self.parameter('Bax_transloc_kf', 1e-2,
                            prior=Normal(-2, 2))
        Bax_transloc_kr = self.parameter('Bax_transloc_kr', 1e-1,
                            prior=Normal(-1, 2))

        Bax = self['Bax']
        Vesicles = self['Vesicles']
        solution = self['solution']
        ves = self['ves']

        self.rule('Bax_translocates_sol_to_ves',
             Bax(loc='c', lipo=None) ** solution +
             Vesicles(bax=None) ** solution <>
             Bax(loc='m', lipo=1) ** ves % Vesicles(bax=1) ** solution,
             Bax_transloc_kf, Bax_transloc_kr)

    def basal_Bax_activation_nonsat(self, reversible=False):
        print "lipo_sites: basal_Bax_activation=%s" % reversible
        Bax = self['Bax']
        Vesicles = self['Vesicles']
        basal_Bax_kf = self.parameter('basal_Bax_kf', 2e-3,
                            prior=Normal(-3, 2))
        self.rule('basal_Bax_activation',
                  Bax(bh3=None, loc='m', lipo=1) % Vesicles(bax=1) >>
                  Bax(bh3=None, loc='i', lipo=None) + Vesicles(bax=None),
                  basal_Bax_kf)
