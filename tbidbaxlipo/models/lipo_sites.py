from pysb import *
import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.util.fitting import fit, mse
from tbidbaxlipo.util import color_iter
from tbidbaxlipo.models import one_cpt, core
from matplotlib.font_manager import FontProperties
from pysb.integrate import odesolve, Solver
from tbidbaxlipo.priors import Normal, UniformLinear

Solver._use_inline = True

class Builder(core.Builder):

    def __init__(self, params_dict=None, nbd_sites=None):
        """Differs from one_cpt.__init__ only in that it allows the
        concentration of vesicles, Vesicles_0, to be estimated."""
        # Sets self.model = Model(), and self.param_dict
        super(Builder, self).__init__(params_dict=params_dict)

        # COMPARTMENTS
        self.cpt_list = ['ves']

        # By default, these parameters are not estimated
        self.parameter('Vesicles_0', 5)
        self.parameter('tBid_0', 20, prior=None)
        self.parameter('Bax_0', 100, prior=None)
        self.parameter('sites_per_liposome', 10., prior=UniformLinear(0, 3))
        self.expression('Vesicles_norm_0',
                        self['Vesicles_0'] * self['sites_per_liposome'])

        self.declare_components()

    def temp(self):
        # Call the lipo-sites-specific implementation of declare_components
        self.declare_components()

        self.initial(tBid(loc='c', bh3=None) ** solution, self['tBid_0'])
        self.initial(Bax(loc='c', bh3=None, a6=None, lipo=None,
                     c3='s', c62='s', c120='s', c122='s', c126='s', c184='s')
                     ** solution, self['Bax_0'])
        self.initial(Vesicles(bax=None) ** solution, self['Vesicles_norm_0'])

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
        if core.LOGGING:
            print("lipo_sites: translocate_Bax()")

        Bax_transloc_kf = self.parameter('Bax_transloc_kf', 1e-2,
                            factor=(1/float(self['sites_per_liposome'].value)),
                            prior=Normal(-3, 1))
        Bax_transloc_kr = self.parameter('Bax_transloc_kr', 1e-1,
                            prior=Normal(-3, 1))

        Bax_mono = self['Bax'](bh3=None, a6=None)
        Vesicles = self['Vesicles']

        cpt_name = self.cpt_list[0]
        self.rule('Bax_mono_translocates_sol_to_%s' % cpt_name,
             Bax_mono(cpt='sol', conf='aq', lipo=None) + Vesicles(bax=None) >>
             Bax_mono(cpt=cpt_name, conf='mem', lipo=1) % Vesicles(bax=1),
             Bax_transloc_kf)
        self.rule('Bax_mono_translocates_%s_to_sol' % cpt_name,
             Bax_mono(cpt=cpt_name, conf='mem', pore='n', lipo=1) %
             Vesicles(bax=1) >>
             Bax_mono(cpt='sol', conf='aq', pore='n', lipo=None) +
             Vesicles(bax=None),
             Bax_transloc_kr)

    def translocate_tBid(self):
        if core.LOGGING:
            print("lipo_sites: translocate_tBid()")

        tBid_transloc_kf = self.parameter('tBid_transloc_kf', 1e-2,
                            factor=(1/float(self['sites_per_liposome'].value)),
                            prior=Normal(-3, 1))
        tBid_transloc_kr = self.parameter('tBid_transloc_kr', 1e-1,
                            prior=Normal(-1, 2))

        tBid = self['tBid']
        Vesicles = self['Vesicles']

        cpt_name = self.cpt_list[0]
        self.rule('tBid_translocates_sol_to_%s' % cpt_name,
             tBid(cpt='sol', conf='aq', lipo=None) + Vesicles(tbid=None) >>
             tBid(cpt=cpt_name, conf='mem', lipo=1) % Vesicles(tbid=1),
             tBid_transloc_kf)
        self.rule('tBid_translocates_%s_to_sol' % cpt_name,
             tBid(cpt=cpt_name, conf='mem', lipo=1) % Vesicles(tbid=1) >>
             tBid(cpt='sol', conf='aq', lipo=None) + Vesicles(tbid=None),
             tBid_transloc_kr)

    def basal_Bax_activation_nonsat(self, reversible=False):
        if core.LOGGING:
            print "lipo_sites: basal_Bax_activation=%s" % reversible
        Bax = self['Bax']
        Vesicles = self['Vesicles']
        basal_Bax_kf = self.parameter('basal_Bax_kf', 2e-3,
                            prior=Normal(-3, 2))
        self.rule('basal_Bax_activation',
                  Bax(bh3=None, loc='m', lipo=1) % Vesicles(bax=1) >>
                  Bax(bh3=None, loc='i', lipo=None) + Vesicles(bax=None),
                  basal_Bax_kf)
