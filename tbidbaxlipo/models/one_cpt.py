"""
Implementation of most basic membrane compartmentalization model, with a
solution compartment and a single homogeneous membrane compartment.

Requirements of the model (these should be incorporated into tests):

- partitioning of tBid into membranes should occur rapidly and lead to nearly
  all (> 90-95%) of tBid in the membrane at equilibrium

Things that are unknown:

- Does the partitioning of tBid and/or Bax to the membrane saturate, or is the
  leveling off due to non-ideal partitioning or aggregation?

The way this is written,

- Bax equilibrates equally between the cytosolic and peripherally bound state
  to both full and empty liposomes.
- But, once in the "inserted" state, it can't dissociate from the vesicle!!
  (i.e., there's no reversal either to the cytosolic or peripherally bound
  states).
- Further, there's no reversal from the porated back to the inserted
  state--instead Bax goes from porated to the peripherally bound state on
  an empty vesicle.

There is a tradeoff between k_pore_rev and k_eflx. If k_pore_rev is slow, this
means that pores, once formed, are very stable, and so the trend is for almost
all Bax to end up as oligomerized. Once pore-state Bax reaches steady state,
dye release occurs basically linearly afterwards.

Also note that in this model, once Bax is in the inserted state, it cannot
revert back to the peripherally bound state without first passing through the
pore state. Also note that it cannot dissociate from the membrane unless it is
in the 'm' (peripherally bound) state.

There is a cycle in that Bax can go from m to m via::

    mtBid + mBax <-> mtBid:mBax (<)-> mtBid + iBax (<)-> pBax
                  K1               K2                 K3

Therefore we need to have
tBid_Bax_kf * tBid_Bax_kc * k_pore = tBid_Bax_kr * (tBid_iBax_kf) * k_pore_rev

The product of the rate constants in both directions should be the same--even
though dye is released in this process, one would imagine that the process
would proceed exactly the same even without dye.

If the whole chain is at equilibrium, then each individual link must be at
equilibrium.

Behavior of the model:---------

- Why does mftBid not come down as the number of empty vesicles goes up?
  Theoretically, as tBid cycles back and forth to the cytoplasm, it should
  unbind, return to dye='none', and then be dye='e' after rebinding an
  empty vesicle.

In the Almeida model, the pore is presumed to have a particular average
lifetime, and the transition from P_lf to P* is determined by some first order
rate, and then The rate of the efflux reaction is proportional to the amount of
pBax(full), the lifetime of the transition between pBax(full) -> pBax(empty),
along with the efflux rate constant.

So imagine that there are precisely 100 (V) vesicles, and 100 Bax molecules.
Now suppose that at a given time there is precisely 1 Bax molecule in the pore
state on a full vesicle, and this is sufficient to trigger dye release. Further
suppose that pore formation is irreversible.  Within that timestep (say, 1
second), you would expect that the pBax(full) concentration would go to 0, and
the CF (efflux) function would change by precisely 1%.

"""

from pysb import *
import numpy as np
from matplotlib import pyplot as plt
from tbidbaxlipo.util.fitting import fit, mse
from tbidbaxlipo.util import color_iter
from tbidbaxlipo.models import core
from matplotlib.font_manager import FontProperties
from pysb.integrate import odesolve, Solver
from bayessb.priors import Normal

Solver._use_inline = True

class Builder(core.Builder):

    def within_compartment_rsf(self):
        return 1.0 / (self['Vesicles_0'].value)

    def __init__(self, params_dict=None, nbd_sites=None, scaling_factor=None):
        # Sets self.model = Model(), and self.param_dict
        core.Builder.__init__(self, params_dict=params_dict)

        self.declare_monomers()

        # INITIAL CONDITIONS
        self.parameter('Vesicles_0', 5, prior=None)
        self.parameter('tBid_0', 20, prior=None)
        self.parameter('Bax_0', 100, prior=None)

        tBid = self['tBid']
        Bax = self['Bax']
        Vesicles = self['Vesicles']
        Pores = self['Pores']

        self.initial(tBid(cpt='sol', conf='aq', bh3=None), self['tBid_0'])
        self.initial(Bax(cpt='sol', conf='aq', bh3=None, a6=None,
                         lipo=None, pore='n'),
                     self['Bax_0'])
        self.initial(Vesicles(bax=None), self['Vesicles_0'])

        # OBSERVABLES
        self.observable('ctBid', tBid(cpt='sol'))
        self.observable('mtBid', tBid(cpt='ves'))
        self.observable('cBax', Bax(cpt='sol'))
        self.observable('mBax', Bax(cpt='ves'))
        self.observable('iBax', Bax(conf='ins'))
        self.observable('tBidBax', tBid(bh3=1) % Bax(bh3=1))
        self.observable('Bax2', Bax(bh3=1) % Bax(bh3=1))
        self.observable('Baxbh3', Bax(bh3=1))
        self.observable('Bax4',
             MatchOnce(Bax(conf='ins', bh3=1, a6=3) %
                       Bax(conf='ins', bh3=1, a6=4) %
                       Bax(conf='ins', bh3=2, a6=3) %
                       Bax(conf='ins', bh3=2, a6=4)))

    # MODEL MACROS
    def dye_release_one_step_catalytic(self):
        """Activated Bax acts as an enzyme that converts liposomes F->E.

        There is no binding assumed--it's a catalytic mechanism.
        """
        pass

    def translocate_tBid_Bax(self):
        print("one_cpt: translocate_tBid_Bax()")
        self.translocate_tBid()
        self.translocate_Bax()

    def translocate_tBid(self):
        print("one_cpt: translocate_tBid()")

        tBid_transloc_kf = self.parameter('tBid_transloc_kf', 1e-1,
                                          prior=None)
        tBid_transloc_kr = self.parameter('tBid_transloc_kr', 1e-1,
                                          prior=Normal(-1, 2))

        Vesicles = self['Vesicles']
        tBid = self['tBid']

        self.rule(
             'tBid_translocates_sol_to_ves',
             tBid(cpt='sol') + Vesicles() >> tBid(cpt='ves') + Vesicles(),
             tBid_transloc_kf)
        self.rule(
             'tBid_translocates_ves_to_sol',
             tBid(cpt='ves', bh3=None) >> tBid(cpt='sol', bh3=None),
             tBid_transloc_kr)

    def translocate_Bax(self):
        print("one_cpt: translocate_Bax()")

        Bax_transloc_kf = self.parameter('Bax_transloc_kf', 1e-2,
                            prior=Normal(-3, 1))
        Bax_transloc_kr = self.parameter('Bax_transloc_kr', 1e-1,
                            prior=Normal(-1, 1))

        Bax_mono = self['Bax'](bh3=None, a6=None)
        Vesicles = self['Vesicles']

        self.rule('Bax_mono_translocates_sol_to_ves',
             Bax_mono(cpt='sol') + Vesicles() >>
             Bax_mono(cpt='ves') + Vesicles(),
             Bax_transloc_kf)
        self.rule('Bax_mono_translocates_ves_to_sol',
             Bax_mono(cpt='ves') >> Bax_mono(cpt='sol'),
             Bax_transloc_kr)

    def translocate_Bax_dimers(self):
        print("one_cpt: translocate_Bax_dimers()")

        param_names = [p.name for p in self.model.parameters]
        if 'Bax_transloc_kf' not in param_names or \
           'Bax_transloc_kr' not in param_names:
            raise Exception("translocate_Bax motif must be called first "
                            "to declare translocation rate parameters.")
        Bax_transloc_kf = self['Bax_transloc_kf']
        Bax_transloc_kr = self['Bax_transloc_kr']
        ves = self['ves']
        solution = self['solution']
        Vesicles = self['Vesicles']
        Bax = self['Bax']

        self.rule('Bax_dimer_translocates_sol_to_%s' % ves.name,
                Bax(loc='c', bh3=1, a6=None) ** solution %
                Bax(loc='c', bh3=1, a6=None) ** solution +
                Vesicles() ** solution >>
                Bax(loc='m', bh3=1, a6=None) ** ves %
                Bax(loc='m', bh3=1, a6=None) ** ves +
                Vesicles() ** solution,
                Bax_transloc_kf)
        self.rule('Bax_dimer_translocates_%s_to_sol' % ves.name,
                Bax(loc='m', bh3=1, a6=None) ** ves %
                Bax(loc='m', bh3=1, a6=None) ** ves >>
                Bax(loc='c', bh3=1, a6=None) ** solution %
                Bax(loc='c', bh3=1, a6=None) ** solution,
                Bax_transloc_kf)

    def pores_from_Bax_monomers(self, bax_loc_state='i', reversible=False):
        """Basically a way of counting the time-integrated amount of
           forward pore formation.
        """
        print("one_cpt: pores_from_Bax_monomers()")

        pore_formation_rate_k = self.parameter('pore_formation_rate_k', 1e-4,
                prior=Normal(-4, 1))

        Bax = self['Bax']
        Pores = self['Pores']
        ves = self['ves']
        solution = self['solution']

        self.rule('Pores_From_Bax_Monomers', 
             Bax(loc=bax_loc_state) >> Bax(loc='p') + Pores(),
             pore_formation_rate_k)

        if reversible:
            pore_reverse_rate_k = self.parameter('pore_reverse_rate_k', 1e-5,
                    prior=Normal(-5,2))

            self.rule('Pores_reverse',
                 Bax(loc='p') ** ves >> Bax(loc='c') ** solution,
                 pore_reverse_rate_k)

        # Pore formation
        self.observable('pBax', Bax(loc='p'))
        self.observable('pores', Pores())

    def pores_from_Bax_dimers(self, bax_loc_state='i', reversible=False):
        """Basically a way of counting the time-integrated amount of
           forward pore formation."""
        print("one_cpt: pores_from_Bax_dimers()")

        Bax = self['Bax']
        Pores = self['Pores']
        ves = self['ves']
        solution = self['solution']

        pore_formation_rate_k = self.parameter('pore_formation_rate_k', 1e-3,
                                               prior=Normal(-3, 2))

        self.rule('Pores_From_Bax_Dimers',
             Bax(loc=bax_loc_state, bh3=1, a6=None) %
             Bax(loc=bax_loc_state, bh3=1, a6=None) >>
             Bax(loc='p', bh3=1, a6=None) % Bax(loc='p', bh3=1, a6=None) +
             Pores(),
             pore_formation_rate_k)

        if reversible:
            pore_reverse_rate_k = self.parameter('pore_reverse_rate_k', 1e-2,
                                                 prior=Normal(-2, 2))
            self.rule('Pores_reverse',
                 Bax(loc='p', bh3=1, a6=None) ** ves %
                 Bax(loc='p', bh3=1, a6=None) ** ves >>
                 Bax(loc='c', bh3=None, a6=None) ** solution +
                 Bax(loc='c', bh3=None, a6=None) ** solution,
                 pore_reverse_rate_k)

        # Pore formation
        self.observable('pBax', Bax(loc='p'))
        self.observable('pores', Pores())

    def pores_from_Bax_tetramers(self, bax_loc_state='m', reversible=False):
        """Basically a way of counting the time-integrated amount of
           forward pore formation."""
        print("one_cpt: pores_from_Bax_tetramers()")

        Bax = self['Bax']
        Pores = self['Pores']
        ves = self['ves']
        solution = self['solution']

        pore_formation_rate_k = self.parameter('pore_formation_rate_k', 1e-3,
                                               prior=Normal(-3, 2))

        self.rule('Pores_From_Bax_Tetramers', 
             MatchOnce(Bax(loc=bax_loc_state, bh3=1, a6=3) %
             Bax(loc=bax_loc_state, bh3=1, a6=4) %
             Bax(loc=bax_loc_state, bh3=2, a6=3) %
             Bax(loc=bax_loc_state, bh3=2, a6=4)) >>
             MatchOnce(Bax(loc='p', bh3=1, a6=3) %
             Bax(loc='p', bh3=1, a6=4) %
             Bax(loc='p', bh3=2, a6=3) %
             Bax(loc='p', bh3=2, a6=4)) +
             Pores(),
             pore_formation_rate_k)

        if reversible:
            pore_reverse_rate_k = self.parameter('pore_reverse_rate_k', 1e-2,
                                                 prior=Normal(-2, 2))

            self.rule('Pores_reverse',
                MatchOnce(Bax(loc='p', bh3=1, a6=3) ** ves %
                Bax(loc='p', bh3=1, a6=4) ** ves %
                Bax(loc='p', bh3=2, a6=3) ** ves %
                Bax(loc='p', bh3=2, a6=4) ** ves) >>
                MatchOnce(Bax(loc='c', bh3=1, a6=3) ** solution %
                Bax(loc='c', bh3=1, a6=4) ** solution %
                Bax(loc='c', bh3=2, a6=3) ** solution %
                Bax(loc='c', bh3=2, a6=4) ** solution),
                pore_reverse_rate_k)

        # Pore formation
        self.observable('pBax', Bax(loc='p'))
        self.observable('pores', Pores())

    # RUNNING THE MODEL
    def run_model(self, tmax=12000, figure_ids=[0,1]):

        t = np.linspace(0, tmax, 1000)
        x = odesolve(self.model, t)

        ci = color_iter()

        tBid_0 = self['tBid_0']
        Bax_0 = self['Bax_0']
        Vesicles_0 = self['Vesicles_0']

        plt.ion()

        # Translocation
        plt.figure(figure_ids[0])
        plt.plot(t, (x['ctBid'])/tBid_0.value, label='ctBid', color=ci.next())
        plt.plot(t, (x['mtBid'])/tBid_0.value, label='mtBid', color=ci.next())
        plt.plot(t, (x['cBax'])/Bax_0.value, label='cBax', color=ci.next())
        plt.plot(t, (x['mBax'])/Bax_0.value, label='mBax', color=ci.next())

        # Activation
        #plt.plot(t, ((x['Baxc62']/Bax_0.value) *
        #             self.model.parameters['c62_scaling'].value),
        #         label='Baxc62', color=ci.next())
        plt.plot(t, x['iBax']/Bax_0.value, label='iBax', color=ci.next())
        plt.plot(t, x['tBidBax']/tBid_0.value, label='tBidBax', color=ci.next())

        # Pore formation
        plt.plot(t, (x['pBax'])/Bax_0.value, label='pBax', color=ci.next())

        ci = color_iter()

        # Dye release expected via poisson assumption ------
        pores_per_ves = x['pores']/Vesicles_0.value
        dye_release = (1 - np.exp(-pores_per_ves)) # Per Schwarz's analysis
        plt.plot(t, dye_release, color=ci.next(), label='dye release')

        fontP = FontProperties()
        fontP.set_size = ('small')
        plt.legend(loc='upper center', prop=fontP, ncol=5,
                   bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=True)

        plt.show()

        #legend()
        #xlabel("Time (seconds)")
        #ylabel("Normalized Concentration")

        # Plot pores/vesicle in a new figure ------
        ci = color_iter()
        plt.figure(figure_ids[1])
        plt.plot(t, (x['pores']/Vesicles_0.value), label='pores/ves',
                 color=ci.next())
        plt.legend(loc='upper center', prop=fontP, ncol=1,
                   bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=True)
        plt.show()
        #legend()
        #xlabel("Time (seconds)")
        #ylabel("Pores/vesicle")

        #xlabel("Time (seconds)")
        #ylabel("Expected dye release")
        #title("Dye release expected based on avg pores/ves")

        #plot(t, x['metBid']/tBid_0.value, label='metBid')
        #plot(t, x['mtBid']/tBid_0.value, label='mtBid')
        #plot(t, x['mftBid']/tBid_0.value, label='mftBid')
        #plot(t, x['mfBax']/Bax_0.value, label='mfBax')
        #plot(t, x['meBax']/Bax_0.value, label='meBax')
        #plot(t, x['tBidBax']/Bax_0.value, label='tBidBax')
        #plot(t, x['tBidiBax']/Bax_0.value, label='tBidiBax')
        #plot(t, (x['iBax']+x['pBax'])/Bax_0.value, label='ipBax')
        #plot(t, (x['eiBax'])/Bax_0.value, label='eiBax')
        #plot(t, (2*x['ePore'])/Bax_0.value, label='ePore')
        #plot(t, (2*x['fPore'])/Bax_0.value, label='fPore')
        #plot(t, (x['pBax'])/Bax_0.value, label='pBax')
        #plot(t, (2*(x['Bax2']+(2*x['Bax4'])))/Bax_0.value, label='Bax2')
        #plot(t, (2*x['Bax2'])/Bax_0.value, label='Bax2')
        #plot(t, (4*x['Bax4'])/Bax_0.value, label='Bax4')
        #plot(t, x['pBax']/Bax_0.value, label='pBax')
        #plot(t, x['eVes']/Vesicles_0.value, label='eVes')

        #figure()
        #plot(t, (x['pores']/model.parameters['NUM_VESICLES'].value),
        #label='pores')
        #xlabel("Time (seconds)")
        #ylabel("Pores per vesicle")

        return [t, x]

    #  COMMENTS
