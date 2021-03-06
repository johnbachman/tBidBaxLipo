"""
Basic implementation of interactions and mechanisms for tBid/Bax interactions
with each other and with membranes. Contains macros to build up a variety
of alternative models of tBid/Bax interactions.

Code organization
=================

Models of tBid/Bax interactions draw on the same pool of mechanistic facts, but
may differ in their implementation of how interactions with membranes are
modeled. As a result, some rules (e.g., those that specify how proteins
translocate to membranes) have to be written in a compartment-specific way,
whereas others can be written in a generic way. To manage this degree of
complexity, macros implementing tBid/Bax interactions are implemented here in a
generic way (e.g.,
:py:meth:`tbidbaxlipo.models.core.Builder.tBid_activates_Bax`); if
compartment-specific modifications are necessary, these functions are overriden
in child classes.

The model builder classes contained in this file and also in
:py:mod:`tbidbaxlipo.models.one_cpt` and :py:mod:`tbidbaxlipo.models.n_cpt` are
used to manage the process of building up alternative models. Each of these
modules contains a class `Builder` that contains an instance of a PySB model.
The classes also contain a series of methods implementing small sub-pieces of
mechanism that can be termed "motifs". These motifs can be recombined in
different ways to create different models.

Individual mechanisms implemented
=================================

tBid Binding to Bax
-------------------

- :py:meth:`tbidbaxlipo.models.core.Builder.tBid_activates_Bax`

- tBid binds to Bax at the a6 (rear) site to trigger conformational change.

- tBid binds to Bax at the bh3 (front) site to trigger conformational change
  (perhaps more appropriate for Bak than Bax).

- tBid binds first at the rear site, then once conformational change is
  triggered, to the BH3 site (groove).

tBid Dissociation/Bax conformational Change
-------------------------------------------

- tBid binding at the rear triggers a conformational change in Bax that makes
  the BH3 exposed but inhibits further binding by tBid at that site.

- tBid binding at the rear triggers a conformational change that exposes the
  BH3 but allows tBid to stay bound

- Perhaps Bax BH3:groove dimerization triggers the conformational change in 
  the rear site that prevents further tBid binding, but tBid can stay bound
  until then.

It will be important to see how these alternatives affect predictions about
tBid/Bax interactions as they related to Andrews FRET experiments.

Bax Insertion
-------------

Many alternatives are possible in terms of insertion of Bax helices, including
N-terminus first or C-terminus first, or a5-a6 insertion preceding activation
or dimerization. Perhaps the way to model this is to represent these as
distinct insertion events that can be inserted into the pathway in any order--
however this will require some trickery to get the preconditions of the rules
to work correctly.

Possible conditioning could include the bound states, either to a tBid, to
another Bax in a dimer, or to Baxes at the tetramerization site; it could
also obviously include which other insertion events have already taken place.

Bax Oligomerization
-------------------

- Apparent consensus model: Bax (and Bak) dimerize first by the BH3:groove
  interface, which makes the a6:a6 interface accessible, by conformational
  change. They then bind by the a6:a6 interface to form higher-order oligomers.
  The oligomers are not bounded in size.

- The "daisy-chain" model, perhaps one in which the Bax/Bak BH3 engages the
  rear pocket of an adjacent Bax/Bak in the chain.

- The "daisy-chain" model with sequential BH3:groove interfaces that are not
  reciprocal 

Bcl-XL and Retrotranslocation
-----------------------------

- Basic model: BclXL/Bcl2 binds the 

Other Notes
===========

Note on vesicle concentration SATSOURA: It was found that the measured r was
often higher than the estimated r, with large variations observed from
preparation-to-preparation.  This discrepancy was mainly the result of a lipid
concentration that was lower than expected, showing that during the
resuspension and extrusion steps, not unexpectedly a significant amount of
lipid was not incorporated into liposomes. This is one reason why different
repeats of the same experiment might give different values for bound Bax, since
this value depends on r, and the actual r can vary from
experiment-to-experiment.

SATSOURA: Immunoblotting suggests that in the absence of tBid, less than ~10%
of Bax or EGFP-Bax is stably bound to liposomes (data not shown and [50]).
However, this does not rule out the possibility that some Bax may transiently
bind to membranes in the absence of tBid, an interac- tion that could be
disrupted during gel filtration. Indeed, the presence of lipid membranes alone
can induce a transient Bax conformational change [17], and in healthy cells a
small fraction of Bax is found loosely associated with mitochondrial membranes
(as shown by the fact that this fraction is carbonate extractable) [10,32]

SATSOURA: However, when observing the diffusion of the EGFP-Bax in the presence
of lipo- somes, but in the absence of tBid (Fig. 2F), it appeared that the
total amount of EGFP-Bax associated with the liposomes is always less than 5%
(for r>1), or that the interaction is too transient (lasting less than a few
ms) to be detected by FCS.

SATSOURA: In contrast, im- munoblotting against tBid showed that for all
proteins-to-liposome ra- tios explored, about 90% of tBid was associated with
the membrane (data not shown).

Simple model of tBid-Bax interaction, consisting only of tBid's activation of
Bax. The basic enzyme-substrate model.

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
  state--instead Bax goes from porated to the peripherally bound state on an
  empty vesicle.

There is a tradeoff between k_pore_rev and k_eflx. If k_pore_rev is slow, this
means that pores, once formed, are very stable, and so the trend is for almost
all Bax to end up as oligomerized. Once pore-state Bax reaches steady state,
dye release occurs basically linearly afterwards.

Also note that in this model, once Bax is in the inserted state, it cannot
revert back to the peripherally bound state without first passing through the
pore state. Also note that it cannot dissociate from the membrane unless it is
in the 'm' (peripherally bound) state.

There is a cycle in that Bax can go from m to m via

mtBid + mBax <-> mtBid:mBax (<)-> mtBid + iBax (<)-> pBax

Therefore we need to have
tBid_Bax_kf * tBid_Bax_kc * k_pore = tBid_Bax_kr * (tBid_iBax_kf) * k_pore_rev

The product of the rate constants in both directions should be the same--even
though dye is released in this process, one would imagine that the process
would proceed exactly the same even without dye.

If the whole chain is at equilibrium, then each individual link must be at
equilibrium.

Behavior of the model:

- Why does mftBid not come down as the number of empty vesicles goes up?
  Theoretically, as tBid cycles back and forth to the cytoplasm, it should
  unbind, return to dye='none', and then be dye='e' after rebinding an empty
  vesicle.

In the Almeida model, the pore is presumed to have a particular average
lifetime, and the transition from P_lf to P* is determined by some first order
rate, and then The rate of the efflux reaction is proportional to the amount of
pBax(full), the lifetime of the transition between pBax(full) -> pBax(empty),
along with the efflux rate constant.

So imagine that there are precisely 100 (V) vesicles, and 100 Bax molecules.
Now suppose that at a given time there is precisely 1 Bax molecule in the pore
state on a full vesicle, and this is sufficient to trigger dye release. Further
suppose that pore formation is irreversible. Within that timestep (say, 1
second), you would expect that the pBax(full) concentration would go to 0, and
the CF (efflux) function would change by precisely 1%.

Notes on Dye Release
====================

Should also set tBid state?

- How do I figure out how much/how fast to set tBid from f to e? 

- Using same rate as for lipos, the amount of tBid that transitions from f to e
  is proportional to the amount of lipo that transitions.

- So the idea is that tBid ** full_ves is spread equally across these liposomes
  and the decay rate of tBid f->e is the same as the decay rate of vesicles
  from f->e::

    Rule('tBid_Full_to_Empty',
         tBid() ** full_ves  + pore_forming_species ** full_ves >>
         tBid() ** empty_ves + pore_forming_species ** full_ves,
         lipo_eflx, move_connected=True)

Interestingly, when pore formation is irreversible, this model of dye release
(with dimers) produces the same steady state fraction of dye release,
regardless of the amount of vesicles (given the same rates for all of the Bax
and dye release steps). This must be because the decay of the vesicles from
full to empty essentially tracks the Bax dimerization rate???

Part of the problem here though is that the rates of tBid/Bax interaction are
not controlled for the amount of lipid there is. With fewer liposomes, there is
a greater tBid/Bax to liposome ratio, and hence encounter frequencies will be
higher.  So forward rates should be scaled by this concentration; they will
then be in terms of molar ratios of tBid/Bax to liposome/lipid.

In addition, lipid amounts will affect the translocation rates of tBid and Bax
to membranes. This is currently taken into account though in the translocation
steps.

To-do list
==========

.. todo:: List outline of the mechanistic alternatives before implementing them.

.. todo:: Think about how models could be built from rules or macros using tables.

.. todo:: Make it so macros can be invoked from within the motifs in this

paradigm.
"""

__author__ = 'johnbachman'

from pysb import *
import numpy as np
import pysb.builder
from tbidbaxlipo.priors import Normal, Uniform, UniformLinear
import sympy
import math
import itertools

LOGGING = False

class Builder(pysb.builder.Builder):
    """
    Class/constructor documentation here.
    """
    def __init__(self, params_dict=None, nbd_sites=None):
        # Sets self.model = Model(), and self.param_dict
        super(Builder, self).__init__(params_dict=params_dict)

    def declare_components(self):
        """Declares signatures for tBid, Bax, and Vesicles.

        Must be called after the child class constructor has initialized the
        list of compartments.
        """
        tBid = self.monomer('tBid', ['bh3', 'conf', 'cpt', 'lipo'],
                     {'conf': ['aq', 'mem'],
                      'cpt':  ['sol'] + self.cpt_list})
        Bax = self.monomer('Bax',
                     ['bh3', 'a6', 'lipo', 'conf', 'cpt', 'pore', 'dye'],
                     {'conf':  ['aq', 'mem', 'ins'],
                      'cpt': ['sol'] + self.cpt_list,
                      'pore': ['y', 'n'],
                      'dye': ['none', 'nbd', 'dac']})
        Vesicles = self.monomer('Vesicles', ['bax', 'tbid'])
        Pores = self.monomer('Pores', ['cpt'], {'cpt':self.cpt_list})

        self.initial(tBid(cpt='sol', conf='aq', bh3=None, lipo=None),
                     self['tBid_0'])
        self.initial(Vesicles(bax=None, tbid=None), self['Vesicles_norm_0'])

        # Initial condition for WT-Bax
        self.initial(Bax(cpt='sol', conf='aq', bh3=None, a6=None,
                         lipo=None, pore='n', dye='none'),
                     self['Bax_0'])
        # Initial condition for DAC-Bax
        self.parameter('Bax_DAC_0', 0., prior=None)
        self.initial(Bax(cpt='sol', conf='aq', bh3=None, a6=None,
                         lipo=None, pore='n', dye='dac'),
                     self['Bax_DAC_0'])
        # Initial condition for NBD-Bax
        self.parameter('Bax_NBD_0', 0., prior=None)
        self.initial(Bax(cpt='sol', conf='aq', bh3=None, a6=None,
                         lipo=None, pore='n', dye='nbd'),
                     self['Bax_NBD_0'])

        # OBSERVABLES
        self.observable('ctBid', tBid(conf='aq'))
        self.observable('mtBid', tBid(conf='mem'))

        self.observable('cBax', Bax(cpt='sol', conf='aq'))
        self.observable('cBax_NBD', Bax(cpt='sol', conf='aq', dye='nbd'))
        self.observable('mBax', Bax(conf='mem'))
        self.observable('mBax_NBD', Bax(conf='mem', dye='nbd'))
        self.observable('mBax_NBD_mono', Bax(conf='mem', dye='nbd', bh3=None))
        self.observable('mBax_mono', Bax(conf='mem', bh3=None))
        self.observable('iBax', Bax(conf='ins'))
        self.observable('iBax_NBD', Bax(conf='ins', dye='nbd'))
        self.observable('iBax_nopore', Bax(conf='ins', pore='n'))
        self.observable('iBax_nopore_NBD', Bax(conf='ins', pore='n', dye='nbd'))
        self.observable('iBax_mono', Bax(conf='ins', bh3=None))
        self.observable('iBax_mono_NBD', Bax(conf='ins', bh3=None, dye='nbd'))
        self.observable('tBidBax', tBid(bh3=1) % Bax(bh3=1))
        self.observable('tBidBax_NBD', tBid(bh3=1) % Bax(bh3=1, dye='nbd'))
        self.observable('tBidmBax_NBD', tBid(bh3=1) %
                                        Bax(bh3=1, conf='mem', dye='nbd'))
        self.observable('tBidiBax_NBD', tBid(bh3=1) %
                                        Bax(bh3=1, conf='ins', dye='nbd'))
        self.observable('Bax2', Bax(bh3=1) % Bax(bh3=1))
        self.observable('Bax2_NBD', Bax(bh3=1, dye='nbd') % Bax(bh3=1))
        self.observable('Bax4',
                        Bax(conf='ins', bh3=1, a6=3) %
                        Bax(conf='ins', bh3=1, a6=4) %
                        Bax(conf='ins', bh3=2, a6=3) %
                        Bax(conf='ins', bh3=2, a6=4))
        self.observable('Bax2FRET',
                        Bax(bh3=1, dye='dac') % Bax(bh3=1, dye='nbd'))
        # Pore formation
        self.observable('pBax', Bax(pore='y'))
        self.observable('pores', Pores())

        for cpt_name in self.cpt_list:
            self.observable('mtBid_%s' % cpt_name,
                            tBid(conf='mem', cpt=cpt_name))
            self.observable('mBax_%s' % cpt_name,
                            Bax(conf='mem', cpt=cpt_name))
            self.observable('iBax_%s' % cpt_name,
                            Bax(conf='ins', cpt=cpt_name))
            self.observable('tBidBax_%s' % cpt_name,
                            tBid(bh3=1, cpt=cpt_name) %
                            Bax(bh3=1, cpt=cpt_name))
            self.observable('Bax2_%s' % cpt_name,
                            Bax(bh3=1, cpt=cpt_name) %
                            Bax(bh3=1, cpt=cpt_name))
            self.observable('Bax4_%s' % cpt_name,
                            Bax(conf='ins', bh3=1, a6=3, cpt=cpt_name) %
                            Bax(conf='ins', bh3=1, a6=4, cpt=cpt_name) %
                            Bax(conf='ins', bh3=2, a6=3, cpt=cpt_name) %
                            Bax(conf='ins', bh3=2, a6=4, cpt=cpt_name))
            # Pore formation
            self.observable('pBax_%s' % cpt_name, Bax(pore='y', cpt=cpt_name))
            self.observable('pores_%s' % cpt_name, Pores(cpt=cpt_name))


    # -- VIRTUAL FUNCTIONS ----------------------------------------------------
    def within_compartment_rsf(self):
        raise NotImplementedError()

    def run_model(self):
        raise NotImplementedError()

    def get_module(self):
        return self.__module__.split('.')[-1]

    # -- METHODS FOR FITTING/CALIBRATION -----------------------------------
    def declare_nbd_scaling_parameters(self, nbd_sites):
        """
        Adds parameters and observables for fitting the NBD fluorescence data.

        Nominal values for the scaling parameters, as well as means and
        variances for their prior distributions, are set here.

        Note that the observables that are associated with the different
        NBD fluorescence states represent assumptions of the model, so they
        are subject to revision.

        Parameters
        ----------
        nbd_sites : list of strings
            A list of strings containing one or more of 'c3', 'c62', 'c120',
            'c122', 'c126'. Any entries other than these are ignored.

        .. todo:: Make core.Builder.declare_nbd_scaling_parameters flexible
            It needs to be able to take parameters that say which scaling
            parameters should be added. Also, it needs to have a way to
            specify which observables in the mechanistic model map to the
            c3, c62, etc. observables.
        """

        Bax = self['Bax']
        site_set = False

        if 'c3' in nbd_sites: 
            self.parameter('c3_scaling', 0.8, mean=np.log10(0.8), variance=0.04,
                           prior=Normal(0, 2))
            site_set = True
        if 'c62' in nbd_sites:
            self.parameter('c62_scaling', 0.9204, mean=np.log10(0.9204),
                           variance=0.04, prior=Normal(0, 2))
            site_set = True
        if 'c120' in nbd_sites:
            self.parameter('c120_scaling', 0.975, mean=np.log10(0.975),
                           variance=0.04, prior=Normal(0, 2))
            site_set = True
        if 'c122' in nbd_sites:
            self.parameter('c122_scaling', 0.952, mean=np.log10(0.952),
                           variance=0.04, prior=Normal(0, 2))
            site_set = True
        if 'c126' in nbd_sites:
            self.parameter('c126_scaling', 0.966, mean=np.log10(0.966),
                           variance=0.04, prior=Normal(0, 2))
            site_set = True

        if not site_set:
            raise Exception('Failed to set any NBD scaling parameters!')

    def prior(self, mcmc, position):
        """Calculate the prior probability of the set of parameters.

        The builder maintains a list of prior objects that corresponds in its
        order to the list of parameters in the parameter list (it includes only
        the parameters to be estimated).  To calculate the prior, the method
        iterates over the objects in self.priors and for each one gets the
        negative log probability of that position in parameter space. The
        negative log probabilities are then summed (equivalent to the product
        of the individual priors) and then returned by the function.
        """
        prior_prob = 0
        for i, prior in enumerate(self.priors):
            prior_prob += prior.pdf(position[i])
        return prior_prob

    def random_initial_values(self):
        """Return a set of random initial values for parameter estimation.

        Uses the prior distributions specified for the parameters to be
        estimated.

        Note that this must be called after the model has been built using one
        of the model building methods--otherwise the parameters (and the set of
        priors) won't have been created!

        Returns
        -------
        An array of initial parameter values in the same order as the
        parameters in self.estimate_params, and drawn from the prior
        distributions specified in self.priors. Note that the initial values
        are given in linear, rather than log, space.
        """
        initial_values_log = np.empty(len(self.priors))
        for i, prior in enumerate(self.priors):
            initial_values_log[i] = prior.random()
        return 10 ** initial_values_log

    # -- MECHANISTIC MOTIFS ------------------------------------------------
    def translocate_tBid_Bax(self):
        if LOGGING:
            print("core: translocate_tBid_Bax()")
        self.translocate_tBid()
        self.translocate_Bax()

    def translocate_tBid(self):
        if LOGGING:
            print("core: translocate_tBid()")

        assert len(self.cpt_list) != 0
        if len(self.cpt_list) > 1:
            # In a multi-compartment situation, we to rescale the forward
            # translocation rate by the stochastic scaling factor
            tBid_transloc_kf = self.parameter('tBid_transloc_kf', 1e-2,
                                          prior=Normal(-3, 1),
                                          factor=(1/float(self.scaling_factor)))
        else:
            # In a single-compartment situation, we need to multiply the
            # forward translocation rate by the concentration of Vesicles
            tBid_transloc_kf = self.parameter('tBid_transloc_kf', 1e-2,
                                              prior=Normal(-3, 1))

        tBid_transloc_kr = self.parameter('tBid_transloc_kr', 1e-1,
                                          prior=Normal(-1, 2))

        Vesicles = self['Vesicles']
        tBid = self['tBid']

        if len(self.cpt_list) == 1:
            cpt_name = self.cpt_list[0]
            self.rule('tBid_translocates_sol_to_%s' % cpt_name,
                 tBid(cpt='sol', conf='aq') + Vesicles() >>
                 tBid(cpt=cpt_name, conf='mem') + Vesicles(),
                 tBid_transloc_kf)
            self.rule('tBid_translocates_%s_to_sol' % cpt_name,
                 tBid(cpt=cpt_name, conf='mem') >> tBid(cpt='sol', conf='aq'),
                 tBid_transloc_kr)
        else:
            for cpt_name in self.cpt_list:
                self.rule('tBid_translocates_sol_to_%s' % cpt_name,
                     tBid(cpt='sol', conf='aq') >>
                     tBid(cpt=cpt_name, conf='mem'),
                     tBid_transloc_kf)
                self.rule('tBid_translocates_%s_to_sol' % cpt_name,
                     tBid(cpt=cpt_name, conf='mem') >>
                     tBid(cpt='sol', conf='aq'),
                     tBid_transloc_kr)


    def translocate_Bax(self):
        if LOGGING:
            print("core: translocate_Bax()")

        assert len(self.cpt_list) != 0
        if len(self.cpt_list) > 1:
            # In a multi-compartment situation, we need to rescale the forward
            # translocation rate by the stochastic scaling factor
            Bax_transloc_kf = self.parameter('Bax_transloc_kf', 1e-2,
                                          prior=Normal(-3, 1),
                                          factor=(1/float(self.scaling_factor)))
        else:
            # In a single-compartment situation, the scaling is provided
            # by the pseudo-first order reaction with Vesicles
            Bax_transloc_kf = self.parameter('Bax_transloc_kf', 1e-2,
                                          prior=Normal(-3, 1))

        Bax_transloc_kr = self.parameter('Bax_transloc_kr', 1e-1,
                            prior=Normal(-1, 1))

        Bax_mono = self['Bax'](bh3=None, a6=None)
        Vesicles = self['Vesicles']

        if len(self.cpt_list) == 1:
            cpt_name = self.cpt_list[0]
            self.rule('Bax_mono_translocates_sol_to_%s' % cpt_name,
                 Bax_mono(cpt='sol', conf='aq') + Vesicles() >>
                 Bax_mono(cpt=cpt_name, conf='mem') + Vesicles(),
                 Bax_transloc_kf)
            self.rule('Bax_mono_translocates_%s_to_sol' % cpt_name,
                 Bax_mono(cpt=cpt_name, conf='mem', pore='n') >>
                 Bax_mono(cpt='sol', conf='aq', pore='n'),
                 Bax_transloc_kr)
        else:
            for cpt_name in self.cpt_list:
                self.rule('Bax_mono_translocates_sol_to_%s' % cpt_name,
                     Bax_mono(cpt='sol', conf='aq') >>
                     Bax_mono(cpt=cpt_name, conf='mem'),
                     Bax_transloc_kf)
                self.rule('Bax_mono_translocates_%s_to_sol' % cpt_name,
                     Bax_mono(cpt=cpt_name, conf='mem', pore='n') >>
                     Bax_mono(cpt='sol', conf='aq', pore='n'),
                     Bax_transloc_kr)


    def tBid_binds_Bax(self, bax_site, bax_conf):
        """
        Can be used for binding to mem or ins state, at any site on Bax.

        Notes
        -----
        When iBax binds tBid at BH3: THERE IS A PROBLEM WITH THIS!!!

        When tBid is bound to Bax, it is prevented from recirculating back to
        the solution.  Therefore you get a sequence of events like:

        - tBid + Bax (full liposome) -> tBid + Bax(i) (full liposome)
        - Bax(p) (full) -> Bax(p) (empty) THIS HAPPENS VERY FAST
        - But also: Bax(i) + tBid (full liposome).

        When complexed in this way, tBid does not recycle back to the solution.
        Therefore tBid catalyses the creation of a species Bax(i) which over
        time shifts the equilibrium of tBid from c to full liposomes.

        However, if you are not interested in the permeabilization status of
        liposomes (e.g., if you are using the model only for looking at
        insertion rates) then this could be OK.

        Rates used in old version, iBax_binds_tBid_at_bh3: kf 1e-1, kr 2.5
        """

        if LOGGING:
            print('core: tBid_binds_Bax(bax_site=%s, bax_conf=%s)' %
                  (bax_site, bax_conf))

        # Forward rate of tBid binding to Bax
        kf = self.parameter('tBid_Bax_%s_%s_kf' % (bax_conf, bax_site), 1e-2,
                            factor=self.within_compartment_rsf(),
                            prior=Normal(-2, 2))
        # Reverse rate of tBid binding to Bax
        kr = self.parameter('tBid_Bax_%s_%s_kr' % (bax_conf, bax_site), 1.5,
                            prior=Normal(0, 2))

        # Create the dicts to parameterize the site that tBid binds to
        bax_site_bound = {bax_site:1}
        bax_site_unbound = {bax_site:None}

        tBid = self['tBid']
        Bax = self['Bax']

        for cpt_name in self.cpt_list:
            self.rule('tBid_Bax_%s_%s_bind_%s' %
                      (bax_conf, bax_site, cpt_name),
                      tBid(cpt=cpt_name, bh3=None) +
                      Bax(cpt=cpt_name, conf=bax_conf, **bax_site_unbound) >>
                      tBid(cpt=cpt_name, bh3=1) %
                      Bax(cpt=cpt_name, conf=bax_conf, **bax_site_bound),
                      kf)
            self.rule('tBid_Bax_%s_%s_unbind_%s' %
                      (bax_conf, bax_site, cpt_name),
                      tBid(cpt=cpt_name, bh3=1) %
                      Bax(cpt=cpt_name, conf=bax_conf, **bax_site_bound) >>
                      tBid(cpt=cpt_name, bh3=None) +
                      Bax(cpt=cpt_name, conf=bax_conf, **bax_site_unbound),
                      kr)

    def tBid_activates_Bax(self, bax_site='bh3', bax_bind_conf='mem',
                           bax_active_conf='ins'):
        """Default implementation of Bax activation by tBid.

        Takes arguments:
            - bax_site specifies the name of the site on Bax to which
              the bh3 site on tBid binds.

        Notes
        -----

        - Andrews suggests that tBid/Bax Kd should work out to 25nM (once in
          membrane, presumably)

        - Binding of tBid to Bax during the activation step should be transient

        - The timescale of 3c exposure has a half-life of about 50s (ln 2 /
          1.381e-2)--this could potentially correspond to the "iBax" form,
          though arguably this may be better indicated by the a5/6 insertion
          per the finding of Annis that Bax is multispanning prior to
          oligomerization

        - Binding of the BH3 (presumably by tBid?) occurs with an initial rate
          of ... (check fit)

        - When tBid is added, 50-80% of Bax binds to liposomes, though this
          goes down at high Bax/liposome ratios. Lovell fig S1 suggests that
          this reaches steady state by about 15 minutes.

        - The forward rate should be normalized according to protein/lipid
          ratio in some way

        Since the surface area of each vesicle is the same, the effect of
        changing vesicle concentration on the forward rate will be by altering
        the expected concentration of tBid and Bax per vesicle.  If there is
        one compartment, 100 Bax and 100 tBid, then the rate will be scaled to
        take into account 100 of each.  If there are 100 compartments, then the
        rate should be scaled to take into account 1 of each. The scaling
        should therefore most likely be simple linear scaling by dividing by
        the number of compartments, which is in the same units as the protein
        concentration (e.g., nanomolar).

        In the deterministic case, if the fundamental forward rate of binding
        between tBid and Bax is kf, this should really be normalized by the P/L
        ratio of both proteins. So for example, kf * tBid/ves * Bax/ves This
        because doubling of the vesicle concentration cuts the relative
        concentrations of both proteins by half, and hence scales the forward
        rate correspondingly.  In the SSA case, The forward rate should be
        kf*RSF * tBid * Bax (where the concentrations given are for that
        compartment). Since they represent concentrations On that individual
        compartment, the rate does not need to be normalized by the vesicle
        concentration.
        """

        if LOGGING:
            print('core: tBid_activates_Bax(bax_site=%s, bax_bind_conf=%s, '
                  'bax_active_conf=%s)' %
                  (bax_site, bax_bind_conf, bax_active_conf))

        self.tBid_binds_Bax(bax_site, bax_bind_conf)

        # Dissociation of tBid from iBax (EP -> E + P)
        kc = self.parameter('tBid_Bax_%s_%s_kc' % (bax_active_conf, bax_site),
                            1e-1, prior=Normal(-1, 2))

        # Create the dicts to parameterize the site that tBid binds to
        bax_site_bound = {bax_site:1}
        bax_site_unbound = {bax_site:None}

        tBid = self['tBid']
        Bax = self['Bax']

        # tBid dissociates from iBax after activation
        self.rule('tBid_unbinds_Bax_%s_%s' % (bax_active_conf, bax_site),
             tBid(bh3=1) % Bax(conf=bax_bind_conf, **bax_site_bound) >>
             tBid(bh3=None) + Bax(conf=bax_active_conf, **bax_site_unbound),
             kc)

    def tBid_activates_Bax_3step(self, bax_site='bh3', bax_bind_conf='mem',
                                 bax_active_conf='ins'):
        """tBid + Bax <-> tBid:Bax -> tBid:Bax* <-> Bax* + tBid"""

        if LOGGING:
            print('core: tBid_activates_Bax_3step(bax_site=%s, bax_bind_conf=%s, '
                  'bax_active_conf=%s)' %
                  (bax_site, bax_bind_conf, bax_active_conf))

        self.tBid_binds_Bax(bax_site, bax_bind_conf)
        self.tBid_binds_Bax(bax_site, bax_active_conf)

        kc = self.parameter('tBidBax_%s_%s_kc' % (bax_bind_conf, bax_active_conf),
                            1e-1, prior=Normal(-1, 2))

        # Create the dicts to parameterize the site that tBid binds to
        bax_site_bound = {bax_site:1}
        bax_site_unbound = {bax_site:None}

        tBid = self['tBid']
        Bax = self['Bax']

        # tBid:Bax undergoes Bax conformational change
        self.rule('tBidBax_%s_to_%s' % (bax_bind_conf, bax_active_conf),
             tBid(bh3=1) % Bax(conf=bax_bind_conf, **bax_site_bound) >>
             tBid(bh3=1) % Bax(conf=bax_active_conf, **bax_site_bound),
             kc)

    def basal_Bax_activation(self):
        if LOGGING:
            print "core: basal_Bax_activation"
        # Spontaneous rate of transition of Bax from the mitochondrial to the
        # inserted state
        Bax = self['Bax']

        basal_Bax_kf = self.parameter('basal_Bax_kf', 2e-3, prior=Normal(-3,1))

        for cpt_name in self.cpt_list:
            self.rule('basal_Bax_activation_%s' % cpt_name,
                      Bax(cpt=cpt_name, conf='mem', bh3=None, a6=None) >>
                      Bax(cpt=cpt_name, conf='ins', bh3=None, a6=None),
                      basal_Bax_kf)

    def Bax_reverses(self):
        """Reversion of inserted Bax to the peripherally bound form."""

        if LOGGING:
            print('core: Bax_reverses')

        Bax = self['Bax']

        # Reversion of active Bax (P -> S)
        # Estimated at 2e-4 +/- 1e-4 in Shamas-Din doi:10.1038/cddis.2014.234
        krev = self.parameter('iBax_reverse_k', 2e-4, prior=Normal(-4, 2))

        # iBax reverses back to mBax
        self.rule('iBax_reverses',
             Bax(conf='ins', bh3=None, a6=None) >>
             Bax(conf='mem', bh3=None, a6=None),
             krev)

    def tBid_reverses_Bax():
        # Rate of the EP->ES transition # FIXME
        Parameter('iBaxtBid_to_mBaxtBid_k', 1e-3)

        # REVERSIBILITY OF BAX ACTIVATION BY TBID (EP -> ES)
        Rule('iBaxtBid_to_mBaxtBid',
             tBid(loc='m', bh3=1) % Bax(loc='i', bh3=1) >>
             tBid(loc='m', bh3=1) % Bax(loc='m', bh3=1),
             iBaxtBid_to_mBaxtBid_k)

    def Bax_dimerizes(self, bax_conf='ins', reversible=True):
        """
        Notes
        -----
        - If the late phase of the 62c signal is an indication of dimerization,
          it starts to manifest around 500s.
        - In SATSOURA, Fig 4. appears to indicate that the Bax-Bax FRET reaches
          steady-state at around 12 minutes.
        """

        if LOGGING:
            print("core: Bax_dimerizes(bax_conf=%s, reversible=%s)" %
              (bax_conf, reversible))

        Bax = self['Bax']

        # Rate of dimerization formation/oligomerization of activated Bax (s^-1)
        Bax_dimerization_kf = self.parameter(
               'Bax_%s_dimerization_kf' % bax_conf, 1e-2,
               factor=self.within_compartment_rsf(),
               prior=Normal(-3, 2))

        self.rule('Bax_%s_forms_dimers_fwd' % bax_conf,
             Bax(cpt='ves', conf=bax_conf, bh3=None, a6=None) +
             Bax(cpt='ves', conf=bax_conf, bh3=None, a6=None) >>
             Bax(cpt='ves', conf=bax_conf, bh3=1, a6=None) %
             Bax(cpt='ves', conf=bax_conf, bh3=1, a6=None),
             Bax_dimerization_kf)

        if reversible:
            Bax_dimerization_kr = self.parameter(
                   'Bax_%s_dimerization_kr' % bax_conf, 1e-3,
                   prior=Normal(-3, 2))

            self.rule('Bax_%s_forms_dimers_rev' % bax_conf,
                 Bax(cpt='ves', conf=bax_conf, bh3=1, a6=None) %
                 Bax(cpt='ves', conf=bax_conf, bh3=1, a6=None) >>
                 Bax(cpt='ves', conf=bax_conf, bh3=None, a6=None) +
                 Bax(cpt='ves', conf=bax_conf, bh3=None, a6=None),
                 Bax_dimerization_kr)

    def Bax_tetramerizes(self, bax_conf='ins'):
        """
        This function depends on Bax_dimerization to be called as well.
        It forms dimers of dimers.

        Notes
        -----
        In Lovell Fig S1, about 80% of the Bax is at membranes (inserted,
        non-extractable, resistant to gel filtration etc.) after
        """

        if LOGGING:
            print("core: Bax_tetramerizes()")

        # Rate of dimerization formation/oligomerization of activated Bax (s^-1)
        Bax_tetramerization_kf = \
                self.parameter('Bax_%s_tetramerization_kf' % bax_conf, 1e-2,
                               factor=self.within_compartment_rsf(),
                               prior=Normal(-2, 2))
        Bax_tetramerization_kr = \
                self.parameter('Bax_%s_tetramerization_kr' % bax_conf, 1e-2,
                               prior=Normal(-2, 2))

        Bax = self['Bax']

        self.rule('Bax_Forms_Tetramers',
             MatchOnce(Bax(cpt='ves', conf=bax_conf, bh3=1, a6=None) %
                       Bax(cpt='ves', conf=bax_conf, bh3=1, a6=None)) +
             MatchOnce(Bax(cpt='ves', conf=bax_conf, bh3=2, a6=None) %
                       Bax(cpt='ves', conf=bax_conf, bh3=2, a6=None)) <>
             MatchOnce(Bax(cpt='ves', conf=bax_conf, bh3=1, a6=3) %
                       Bax(cpt='ves', conf=bax_conf, bh3=1, a6=4) %
                       Bax(cpt='ves', conf=bax_conf, bh3=2, a6=3) %
                       Bax(cpt='ves', conf=bax_conf, bh3=2, a6=4)),
             Bax_tetramerization_kf, Bax_tetramerization_kr)

    def Bax_nmerizes(self, n, bax_conf_state='ins'):
        Bax_nmerization_kf = self.parameter(
                'Bax_%s_%dmerization_kf' % (bax_conf_state, n),
                1e-2, factor=self.within_compartment_rsf(),
                prior=Normal(-2, 2))
        Bax_nmerization_kr = self.parameter(
                'Bax_%s_%dmerization_kr' % (bax_conf_state, n),
                1e-2, factor=self.within_compartment_rsf(),
                prior=Normal(-2, 2))
        Bax = self['Bax']

        for cpt_name in self.cpt_list:
            free_Bax = Bax(cpt=cpt_name, conf=bax_conf_state,
                           bh3=None, a6=None)
            lhs = free_Bax
            rhs = Bax(cpt=cpt_name, conf=bax_conf_state, bh3=0, a6=n-1)
            for i in range(n-1):
                lhs += free_Bax
                rhs = rhs % Bax(cpt=cpt_name, conf=bax_conf_state,
                                bh3=i+1, a6=i)

            self.rule('Bax_%s_forms_%dmers_fwd_%s' %
                      (bax_conf_state, n, cpt_name),
                      lhs >> rhs, Bax_nmerization_kf)

            self.observable('Bax%d_%s_%s' % (n, bax_conf_state, cpt_name), rhs)

        self.rule('Bax_%s_forms_%dmers_rev' % (bax_conf_state, n),
                  rhs >> lhs, Bax_nmerization_kr)

    def pores_from_Bax_schwarz(self, bax_conf='mem'):
        """The Schwarz model has pores form directly from the membrane-bound
        form without the pore former itself changing state. This appears
        a simpler approach that is worth using for comparing methods of
        permeabilization."""
        if LOGGING:
            print("core: pores_from_Bax_schwarz()")

        pore_formation_rate_k = self.parameter('pore_formation_rate_k', 1e-3,
                                               prior=Normal(-4, 1))

        Bax_mono = self['Bax'](bh3=None, a6=None)
        Pores = self['Pores']

        for cpt_name in self.cpt_list:
            self.rule('pores_from_Bax_schwarz_%s' % cpt_name,
                      Bax_mono(cpt=cpt_name, conf=bax_conf) >>
                      Bax_mono(cpt=cpt_name, conf=bax_conf) +
                      Pores(cpt=cpt_name),
                      pore_formation_rate_k)

    def pores_from_Bax_monomers(self, bax_conf='ins'):
        """Basically a way of counting the time-integrated amount of
        forward pore formation.

        The problem is that you need to be able to register that a pore has
        formed in the vesicle, but you don't want the same Bax to form more
        than one pore on the vesicle. If the rule simply said that inserted
        Bax produces pores at some rate, then there would be an accumulation
        of pores even inserted Bax was not increasing.

        So when a pore is formed, we set pore='y' do indicate that a particular
        Bax has formed a pore. That works well until we consider reversibility.
        If the pore is reversible, that is, if the Bax is to be able to revert
        from the inserted state back to the solution, then we must have
        it be able to return to a peripherally bound state without allowing it
        to create another pore again.

        If we keep pore='y' even when Bax is in the peripherally bound state,
        and interpret the semantics of this state as saying that "this Bax
        has formed a pore" rather than "this Bax is in a pore", then we could
        only reset the pore back to pore='n' once the Bax had returned to
        solution, which is addressed by the translocate_Bax macro.

        With that approach, there's no such thing as "pore reversibility"--
        only reversibility of Bax insertion, which is already encoded by
        the macro Bax_reverses().
        """
        if LOGGING:
            print("core: pores_from_Bax_monomers()")

        pore_formation_rate_k = self.parameter('pore_formation_rate_k', 1e-3,
                                               prior=Normal(-4, 1))

        Bax_mono = self['Bax'](bh3=None, a6=None)
        Pores = self['Pores']

        for cpt_name in self.cpt_list:
            self.rule('pores_from_Bax_monomers_%s' % cpt_name,
                      Bax_mono(cpt=cpt_name, conf=bax_conf, pore='n') >>
                      Bax_mono(cpt=cpt_name, conf=bax_conf, pore='y') +
                      Pores(cpt=cpt_name),
                      pore_formation_rate_k)

    def pores_from_Bax_dimers(self, bax_conf='ins'):
        """Basically a way of counting the time-integrated amount of
           forward pore formation.

        See notes for pores_from_Bax_monomers, above.
        """
        if LOGGING:
            print("core: pores_from_Bax_dimers(bax_conf=%s)" % bax_conf)

        Bax = self['Bax']
        Pores = self['Pores']

        pore_formation_rate_k = self.parameter('pore_formation_rate_k', 1e-3,
                                               prior=Normal(-3, 2))

        self.rule('pores_from_Bax_dimers',
             Bax(cpt='ves', conf=bax_conf, bh3=1, a6=None, pore='n') %
             Bax(cpt='ves', conf=bax_conf, bh3=1, a6=None, pore='n') >>
             Bax(cpt='ves', conf=bax_conf, bh3=1, a6=None, pore='y') %
             Bax(cpt='ves', conf=bax_conf, bh3=1, a6=None, pore='y') +
             Pores(cpt='ves'),
             pore_formation_rate_k)

    def pores_from_Bax_nmers(self, n, bax_loc_state='i'):
        pore_formation_rate_k = self.parameter('pore_formation_rate_k', 1e-4,
                prior=Normal(-4,1))

        Bax = self['Bax']
        Pores = self['Pores']

        lhs = Bax(loc=bax_loc_state, bh3=0, a6=n-1)
        rhs = Bax(loc='p', bh3=0, a6=n-1)
        for i in range(n-1):
            lhs = lhs % Bax(loc=bax_loc_state, bh3=i+1, a6=i)
            rhs = rhs % Bax(loc='p', bh3=i+1, a6=i)

        self.rule('pores_from_Bax_%s_%dmers' % (bax_loc_state, n),
             lhs >> rhs + Pores(), pore_formation_rate_k)

    def pores_aggregate(self, bax_conf='mem'):
        if LOGGING:
            print("core: pores_aggregate()")

        pore_aggregation_rate_k = self.parameter('pore_aggregation_rate_k',
                                                 5e-3, prior=Normal(-3, 1))

        Bax_mono = self['Bax'](bh3=None, a6=None)
        Pores = self['Pores']

        for cpt_name in self.cpt_list:
            self.rule('pores_aggregate_Bax_%s_%s' % (bax_conf, cpt_name),
                      Pores(cpt=cpt_name) +
                      Bax_mono(cpt=cpt_name, conf=bax_conf, pore='n') >>
                      Pores(cpt=cpt_name) +
                      Bax_mono(cpt=cpt_name, conf=bax_conf, pore='y'),
                      pore_aggregation_rate_k)

    # -- OTHER MOTIFS --------------------------------------------------------
    def Bax_auto_activates(self, activator_site='bh3', target_site='bh3'):
        """If this rule is ultimately going to be used in cases where
        polymerization is actually allowed, the fact that the product
        (activated) Bax is specified as fully unbound may be a problem.
        """

        if LOGGING:
            print "core: Bax_auto_activates"

        iBax_mBax_kf = self.parameter('iBax_mBax_kf', 1e-3)
        iBax_mBax_kr = self.parameter('iBax_mBax_kr', 10)
        iBaxmBax_to_iBax_iBax_k = self.parameter('iBaxmBax_to_iBax_iBax_k', 100)

        # Some useful aliases
        Bax = self['Bax']
        Bax_act_free = Bax(cpt='ves', conf='ins', **{activator_site:None})
        Bax_act_bound = Bax(cpt='ves', conf='ins', **{activator_site:1})
        Bax_tgt_mem_free = Bax(cpt='ves', conf='mem', **{target_site:None})
        Bax_tgt_mem_bound = Bax(cpt='ves', conf='mem', **{target_site:1})
        Bax_tgt_ins_free = Bax(cpt='ves', conf='ins', **{target_site:None})

        self.rule('iBax_binds_mBax',
                  Bax_act_free + Bax_tgt_mem_free <>
                  Bax_act_bound % Bax_tgt_mem_bound,
                  iBax_mBax_kf, iBax_mBax_kr)
        self.rule('iBax_activates_mBax',
                  Bax_act_bound % Bax_tgt_mem_bound >>
                  Bax_act_free + Bax_tgt_ins_free,
                  iBaxmBax_to_iBax_iBax_k)

    def Bax_auto_activates_one_step(self):
        """In the simplest case, Bax can only auto-activate another Bax if it
        is monomeric--this implies that it is not only bound to tBid, but also
        implies that it cannot be in a pore."""

        if LOGGING:
            print "core: Bax_auto_activates_one_step"

        iBax_activates_mBax_k = self.parameter('iBax_activates_mBax_k', 1e-4,
                                               prior=Normal(-4, 2))
        Bax = self['Bax']

        for cpt_name in self.cpt_list:
            self.rule('iBax_activates_mBax_%s' % cpt_name,
                      Bax(cpt=cpt_name, conf='ins', bh3=None, a6=None) +
                      Bax(cpt=cpt_name, conf='mem', bh3=None, a6=None) >>
                      Bax(cpt=cpt_name, conf='ins', bh3=None, a6=None) +
                      Bax(cpt=cpt_name, conf='ins', bh3=None, a6=None),
                      iBax_activates_mBax_k)

    def basal_bh3_exposure_auto2(self):
        if LOGGING:
            print "core: basal_bh3_exposure_auto2"

        Bax = self['Bax']

        #basal_Bax_kf = self.parameter('basal_Bax_kf', 2e-3)
        #basal_Bax_kr = self.parameter('basal_Bax_kr', 2e-3)

        #self.rule('basal_Bax_activation',
        #          Bax(bh3=None, loc='m') >> Bax(bh3=None, loc='bh3expos'),
        #          basal_Bax_kf, basal_Bax_kr)

        # Andrews suggests that tBid/Bax Kd should work out to 25nM
        # Forward rate of iBax binding to Bax (E + S -> ES)
        bh3Bax_mBax_kf = self.parameter('bh3Bax_mBax_kf', 1e-4)
        # Reverse rate of iBax binding to Bax (ES -> E + S)
        bh3Bax_mBax_kr = self.parameter('bh3Bax_mBax_kr', 50e-4)
        # Dissociation of iBax from iBax (EP -> E + P)
        bh3BaxmBax_to_bh3Bax_iBax_k = \
                        self.parameter('bh3BaxmBax_to_bh3Bax_iBax_k', 50e-4)
        Bax = self['Bax']

        self.rule('bh3Bax_binds_mBax',
                Bax(loc='m', a6=None, bh3=None) +
                Bax(loc='m', a6=None, bh3=None) <>
                Bax(loc='m', a6=None, bh3=1) %
                Bax(loc='m', a6=1, bh3=None),
                bh3Bax_mBax_kf, bh3Bax_mBax_kr)
        self.rule('bh3BaxmBax_to_bh3Bax_iBax',
                Bax(loc='m', a6=None, bh3=1) %
                Bax(loc='m', a6=1, bh3=None) >>
                Bax(loc='m', a6=None, bh3=None) +
                Bax(loc='i', a6=None, bh3=None),
                bh3BaxmBax_to_bh3Bax_iBax_k)

    ## -- MODEL BUILDING FUNCTIONS -------------------------------------------
    """
    Alternatives for activation: initial activation could occu
    - a6 of bh3, leading to
    - Fully activated state or simply c3 exposed, a9 inserted.
    - tBid can dissociate at the activation step, or not.
    - From there, fully activated Bax can bind tBid again at a6 or bh3;
    - Fully activated Bax could be reversible, or not.
    - Fully activated Bax could dimerize, or not
    - Partially activated (loosely inserted) Bax can bind tBid at a6 or bh3,
      and:
       - do nothing (i.e., just fluoresce) or become fully activated
    - Fully activated Bax can either remain bound to tBid (at a6, or bh3)
      (which can subsequently dissociate) - or dissociate immediately
      leading to
    """

    def build_model_from_dict(self, model_dict):
        model_attribute_sort_order = {
            'builder': 0,
            'bidtranslocation': 1,
            'baxtranslocation': 2,
            'activation': 3,
            'reversal': 4,
            'autoactivation': 5,
            'dimerization': 6,
            'nbd': 7,
            'bidfret': 8,
            'baxfret': 9,
            'bleach': 10,
            'timeoffset': 11,
        }
        model_string = ''

        def unrecognized_implementation(feature, implementation):
            raise ValueError("Don't know how to implement %s for feature %s" %
                             (implementation, feature))

        # Canonicalize the order of the model features based on the sort order
        # defined in the dict
        for feature in sorted(model_dict.keys(),
                              key=lambda x: model_attribute_sort_order[x]):
            implementation = model_dict[feature]
            # Regardless of what the feature is, if the implementation flag
            # is 0, that means don't implement it; skip to the next one
            if implementation == 0:
                continue
            # Build up a useful string notation for the model features
            if feature == 'builder':
                if implementation == 'one_cpt':
                    model_string += '1c_'
                elif implementation == 'lipo_sites':
                    model_string += 'ls_'
            else:
                model_string += feature[0].upper() + feature[1:5] + \
                                str(implementation)
            # Call the appropriate model macros based on the dict entries
            # Builder feature
            if feature == 'builder':
                # We don't have to do anything here because this builder to
                # use should have already been figured out earlier. However,
                # here we can make sure that the builder used goes into the
                # model name
                pass
            # Bax translocation
            elif feature == 'baxtranslocation':
                if implementation == 1:
                    self.translocate_Bax()
                else:
                    unrecognized_implementation(feature, implementation)
            # Bid translocation
            elif feature == 'bidtranslocation':
                if implementation == 1:
                    self.translocate_tBid()
                else:
                    unrecognized_implementation(feature, implementation)
            # Bax activation
            elif feature == 'activation':
                if implementation == 1:
                    self.basal_Bax_activation()
                elif implementation == 2:
                    self.tBid_activates_Bax(bax_site='bh3', bax_bind_conf='mem',
                                            bax_active_conf='ins')
                elif implementation == 3:
                    self.tBid_activates_Bax_3step(bax_site='bh3',
                            bax_bind_conf='mem', bax_active_conf='ins')
                else:
                    unrecognized_implementation(feature, implementation)
            # Bax activation reversal
            elif feature == 'reversal':
                if implementation == 1:
                    self.Bax_reverses()
                else:
                    unrecognized_implementation(feature, implementation)
            # Bax autoactivation
            elif feature == 'autoactivation':
                if implementation == 1:
                    self.Bax_auto_activates_one_step()
                elif implementation == 2:
                    self.Bax_auto_activates()
                else:
                    unrecognized_implementation(feature, implementation)
            # Bax dimerization
            elif feature == 'dimerization':
                if implementation == 1:
                    self.Bax_dimerizes(bax_conf='ins', reversible=False)
                elif implementation == 2:
                    self.Bax_dimerizes(bax_conf='ins', reversible=True)
                else:
                    unrecognized_implementation(feature, implementation)
            # NBD fluorescence
            elif feature == 'nbd':
                c0 = self.parameter('c0_scaling', 1., prior=None)
                c1 = self.parameter('c1_scaling', 5.,
                                    prior=UniformLinear(-1, 1))
                # 1: two confs: all iBax
                if implementation == 1:
                    self.expression('NBD',
                            (c0 * self['cBax_NBD'] +
                            c0 * self['mBax_NBD'] +
                            c1 * self['iBax_nopore_NBD']) / self['Bax_NBD_0'])
                # 2: two confs: dimer only
                elif implementation == 2:
                    self.expression('NBD',
                            (c0 * self['cBax_NBD'] +
                             c0 * self['mBax_NBD'] +
                             c0 * self['iBax_mono_NBD'] +
                             c1 * self['Bax2_NBD'])  / self['Bax_NBD_0'])
                # 3: three confs: 1st iBax (free or tBid bound), 2nd iBax dimer
                elif implementation == 3:
                    c2 = self.parameter('c2_scaling', 5.,
                                        prior=UniformLinear(-1, 1))
                    self.expression('NBD',
                            (c0 * self['cBax_NBD'] +
                             c0 * self['mBax_NBD'] +
                             c1 * self['tBidiBax_NBD'] +
                             c1 * self['iBax_mono_NBD'] +
                             c2 * self['Bax2_NBD']) / self['Bax_NBD_0'])
                # 4: four confs: tBid/Bax, iBax_mono, and dimer
                elif implementation == 4:
                    c2 = self.parameter('c2_scaling', 5.,
                                        prior=UniformLinear(-1, 1))
                    c3 = self.parameter('c3_scaling', 5.,
                                        prior=UniformLinear(-1, 1))
                    self.expression('NBD',
                            (c0 * self['cBax_NBD'] +
                             c0 * self['mBax_NBD'] +
                             c1 * self['tBidBax_NBD'] +
                             c2 * self['iBax_mono_NBD'] +
                             c3 * self['Bax2_NBD']) / self['Bax_NBD_0'])
                # 5: three confs: c/mBax -> tBid:Bax -> Bax*
                elif implementation == 5:
                    c2 = self.parameter('c2_scaling', 5.,
                                        prior=UniformLinear(-1, 1))
                    self.expression('NBD',
                            (c0 * self['cBax_NBD'] +
                             c0 * self['mBax_NBD_mono'] +
                             c1 * self['tBidBax_NBD'] +
                             c2 * self['iBax_NBD']) / self['Bax_NBD_0'])
                else:
                    unrecognized_implementation(feature, implementation)
            # Bid/Bax FRET
            elif feature == 'bidfret':
                bid_fret1 = self.parameter('bid_fret1', 30.,
                                           prior=UniformLinear(0, 2))
                # two states: unbound and bound during activation
                if implementation == 1:
                    self.expression('BidFRET',
                            (bid_fret1 * self['tBidBax_NBD']) / self['tBid_0'])
                # three states: free, tBid:mBax, and tBid:iBax
                elif implementation == 2:
                    bid_fret2 = self.parameter('bid_fret2', 30.,
                                               prior=UniformLinear(0, 2))
                    self.expression('BidFRET',
                            (bid_fret1 * self['tBidmBax_NBD'] +
                             bid_fret2 * self['tBidiBax_NBD']) /
                            self['tBid_0'])
                else:
                    unrecognized_implementation(feature, implementation)
            # Bax/Bax FRET
            elif feature == 'baxfret':
                bax_fret_scaling = self.parameter('bax_fret_scaling', 50.,
                                                  prior=UniformLinear(0, 2))
                if implementation == 1: # all Bax dimers
                    self.expression('BaxFRET',
                            (self['Bax2FRET'] / self['Bax_DAC_0']) *
                            bax_fret_scaling)
                else:
                    unrecognized_implementation(feature, implementation)
            elif feature == 'bleach': # Requires NBD to be present
                if implementation == 1:
                    self.monomer('Bleach', [])
                    self.parameter('Bleach_0', 1.0)
                    self.parameter('Bleach_k', 1.17e-5)
                    self.initial(self['Bleach'](), self['Bleach_0'])
                    self.rule('Bleaching', self['Bleach']() >> None,
                              self['Bleach_k'])
                    self.observable('Bleach_', self['Bleach']())
                    self.expression('NBD_bleach',
                                    self['NBD'] * self['Bleach_'])
                else:
                    unrecognized_implementation(feature, implementation)
            elif feature == 'timeoffset':
                if implementation == 'fit':
                    self.parameter('timeoffset', 100, prior=Uniform(-1, 3))
                else:
                    unrecognized_implementation(feature, implementation)
            else:
                raise ValueError("Don't know how to implement feature %s" %
                                 feature)

        self.model.name = model_string

    def build_model_t(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model t:"
        self.translocate_tBid_Bax()
        self.model.name = 't'

    def build_model_ta(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model ta:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3', bax_bind_conf='mem',
                                bax_active_conf='ins')
        self.model.name = 'ta'

    def build_model_tai(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model tai:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.tBid_binds_Bax(bax_site='bh3', bax_conf='ins')
        self.model.name = 'tai'

    def build_model_taid(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model taid:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.tBid_binds_Bax(bax_site='bh3', bax_conf='ins')
        self.Bax_dimerizes(bax_conf='ins')
        self.model.name = 'taid'

    def build_model_taidt(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model taidt:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.tBid_binds_Bax(bax_site='bh3', bax_conf='ins')
        self.Bax_dimerizes(bax_conf='ins')
        self.Bax_tetramerizes()
        self.model.name = 'taidt'

    def build_model_tair(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model tair:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.tBid_binds_Bax(bax_site='bh3', bax_conf='ins')
        self.Bax_reverses()
        self.model.name = 'tair'

    def build_model_taird(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model taird:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.tBid_binds_Bax(bax_site='bh3', bax_conf='ins')
        self.Bax_reverses()
        self.Bax_dimerizes(bax_conf='ins')
        self.model.name = 'taird'

    def build_model_tairdt(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model tairdt:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.tBid_binds_Bax(bax_site='bh3', bax_conf='ins')
        self.Bax_reverses()
        self.Bax_dimerizes(bax_conf='ins')
        self.Bax_tetramerizes()
        self.model.name = 'tairdt'

    def build_model_tad(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model tad:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.Bax_dimerizes(bax_conf='ins')
        self.model.name = 'tad'

    def build_model_tadt(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model tadt:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.Bax_dimerizes(bax_conf='ins')
        self.Bax_tetramerizes()
        self.model.name = 'tadt'

    def build_model_tar(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model tar:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.Bax_reverses()
        self.model.name = 'tar'

    def build_model_tard(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model tard:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.Bax_reverses()
        self.Bax_dimerizes(bax_conf='ins')
        self.model.name = 'tard'

    def build_model_tardt(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model tardt:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.Bax_reverses()
        self.Bax_dimerizes(bax_conf='ins')
        self.Bax_tetramerizes()
        self.model.name = 'tardt'

    # Model for NBD/Terbium data
    def build_model_nbd_terbium(self):
        #self.build_model_bax_heat_reversible()
        self.build_model_bax_heat()
        # Additional parameters
        pore_scaling = self.parameter('pore_scaling', 100.,
                                prior=Normal(2, 1))
        c0_scaling = self.parameter('c0_scaling', 1,
                                prior=None)
        c1_scaling = self.parameter('c1_scaling', 1.5,
                                prior=Normal(0, 0.5))
        c2_scaling = self.parameter('c2_scaling', 2,
                                prior=Normal(0, 0.5))
        # Get observables
        self.expression('ScaledPores', self['pores'] / pore_scaling)
        self.expression('NBD', (c0_scaling * (self['cBax'] + self['mBax']) +
                               c1_scaling * self['iBax'] +
                               c2_scaling * self['pBax']) / self['Bax_0'])
        self.model.name = 'nbd_terbium'

    # Model for NBD/Terbium data
    def build_model_nbd_terbium_c2_pores(self):
        #self.build_model_bax_heat_reversible()
        self.build_model_bax_heat_reversible()
        # Additional parameters
        pore_scaling = self.parameter('pore_scaling', 100.,
                                prior=Normal(2, 2))
        c0_scaling = self.parameter('c0_scaling', 1,
                                prior=None)
        c1_scaling = self.parameter('c1_scaling', 1.5,
                                prior=Normal(0, 0.5))
        c2_scaling = self.parameter('c2_scaling', 2,
                                prior=Normal(0, 0.5))
        # Get observables
        self.expression('Tb_sympy', 1 - sympy.exp(-self['pBax'] / pore_scaling))
        self.expression('Tb_math', 1 - math.e ** (-self['pBax'] / pore_scaling))
        #self.expression('Tb', self['pBax'] / pore_scaling)
        self.expression('NBD', (c0_scaling * (self['cBax'] + self['mBax']) +
                               c1_scaling * self['iBax'] +
                               c2_scaling * self['pBax']) / self['Bax_0'])
        self.model.name = 'nbd_terbium_c2_reversible_pores_exp'

    # Model for NBD + Bid/Bax FRET data
    def build_model_bid_bax_fret(self):
        self.translocate_tBid_Bax()
        self.tBid_binds_Bax('bh3', 'mem')

    def build_model_nbd_3_conf(self):
        self.build_model_bax_heat()
        c0 = self.parameter('c0_scaling', 1., prior=None)
        c1 = self.parameter('c1_scaling', 2., prior=Uniform(0, 1))
        c2 = self.parameter('c2_scaling', 5., prior=Uniform(0, 1))
        self.expression('NBD',
                (c0 * self['cBax'] +
                c0 * self['mBax'] +
                c1 * self['iBax_nopore'] +
                c2 * self['pBax']) / self['Bax_0'])

    def build_model_nbd_2_conf(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        c0 = self.parameter('c0_scaling', 1., prior=None)
        c1 = self.parameter('c1_scaling', 5., prior=Uniform(0, 1))
        self.expression('NBD',
                (c0 * self['cBax'] +
                c0 * self['mBax'] +
                c1 * self['iBax_nopore']) / self['Bax_0'])
        self.model.name = 'nbd_2_conf'

    def build_model_nbd_2_conf_dimer(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        #self.Bax_reverses()
        self.Bax_dimerizes(bax_conf='ins', reversible=True)
        #self.Bax_dimerizes(bax_conf='ins', reversible=False)
        c0 = self.parameter('c0_scaling', 1., prior=None)
        c1 = self.parameter('c1_scaling', 5., prior=Uniform(0, 1))
        # For 140320, special case
        self.monomer('Bleach', [])
        self.parameter('Bleach_0', 1.0)
        self.parameter('Bleach_k', 1.17e-5)
        self.initial(self['Bleach'](), self['Bleach_0'])
        self.rule('Bleaching', self['Bleach']() >> None, self['Bleach_k'])
        self.observable('Bleach_', self['Bleach']())
        # Expression for NBD fluorescence
        #self.expression('NBD',
        #        ((c0 * self['cBax'] +
        #         c0 * self['mBax'] +
        #         c1 * self['iBax_nopore']) / self['Bax_0']) *
        #         self['Bleach_'])
        #self.expression('NBD',
        #        (c0 * self['cBax'] +
        #        c0 * self['mBax_mono'] +
        #        c0 * self['iBax_mono'] +
        #        c1 * self['Bax2']) / self['Bax_0'])
        self.expression('NBD',
                ((c0 * self['cBax'] +
                 c0 * self['mBax'] +
                 c0 * self['iBax_mono'] +
                 c1 * self['Bax2'])  / self['Bax_0']) *
                 self['Bleach_'])
        self.model.name = 'nbd_2_conf_dimer'

    def build_model_nbd_2_conf_rev(self):
        self.build_model_nbd_2_conf()
        self.Bax_reverses()
        self.model.name = 'nbd_2_conf_rev'

    def build_model_nbd_2_conf_auto(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates_one_step()
        c0 = self.parameter('c0_scaling', 1., prior=None)
        c1 = self.parameter('c1_scaling', 5., prior=Uniform(0, 1))
        self.expression('NBD',
                (c0 * self['cBax'] +
                c0 * self['mBax'] +
                c1 * self['iBax_nopore']) / self['Bax_0'])
        self.model.name = 'nbd_2_conf_auto'

    def build_model_nbd_2_conf_auto_rev(self):
        self.build_model_nbd_2_conf_auto()
        self.Bax_reverses()
        self.model.name = 'nbd_2_conf_auto_rev'

    # Models incorporating dye release
    def build_model_bax_heat(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.pores_from_Bax_monomers(bax_conf='ins')
        self.model.name = 'bax_heat'

    def build_model_bax_heat_reversible(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_reverses()
        self.pores_from_Bax_monomers(bax_conf='ins')
        self.model.name = 'bax_heat_reversible'

    def build_model_bax_heat_dimer(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_dimerizes(bax_conf='ins')
        self.pores_from_Bax_dimers(bax_conf='ins')
        self.model.name = 'bax_heat_dimer'

    def build_model_bax_heat_dimer_reversible(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_reverses()
        self.Bax_dimerizes(bax_conf='ins')
        self.pores_from_Bax_dimers(bax_conf='ins')
        self.model.name = 'bax_heat_dimer_reversible'

    def build_model_bax_heat_auto1(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates_one_step()
        self.pores_from_Bax_monomers(bax_conf='ins')
        self.model.name = 'bax_heat_auto1'

    def build_model_bax_heat_auto2(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates()
        self.pores_from_Bax_monomers(bax_conf='ins')
        self.model.name = 'bax_heat_auto2'

    def build_model_bax_heat_bh3_auto2(self):
        self.translocate_Bax()
        self.bh3_exposure_auto2()
        self.pores_from_Bax_monomers(bax_conf='ins')
        self.model.name = 'bax_heat_auto2'

    def build_model_bax_heat_auto1_reversible_activation(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_reverses()
        self.Bax_auto_activates_one_step()
        self.pores_from_Bax_monomers(bax_conf='ins')
        self.model.name = 'bax_heat_auto1_reversible_activation'

    def build_model_bax_heat_auto2_reversible_activation(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_reverses()
        self.Bax_auto_activates()
        self.pores_from_Bax_monomers(bax_conf='ins')
        self.model.name = 'bax_heat_auto2_reversible_activation'

    def build_model_bax_heat_auto1_reversible(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates_one_step()
        self.Bax_reverses()
        self.pores_from_Bax_monomers(bax_conf='ins')
        self.model.name = 'bax_heat_auto1_reversible'

    def build_model_bax_heat_auto2_reversible(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates()
        self.Bax_reverses()
        self.pores_from_Bax_monomers(bax_conf='ins')
        self.model.name = 'bax_heat_auto2_reversible'

    def build_model_bax_heat_auto1_dimer(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates_one_step()
        self.Bax_dimerizes()
        self.pores_from_Bax_dimers(bax_conf='ins')
        self.model.name = 'bax_heat_auto1_dimer'

    def build_model_bax_heat_auto2_dimer(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates()
        self.Bax_dimerizes()
        self.pores_from_Bax_dimers(bax_conf='ins')
        self.model.name = 'bax_heat_auto2_dimer'

    def build_model_bax_heat_auto1_dimer_reversible(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates_one_step()
        self.Bax_dimerizes()
        self.Bax_reverses()
        self.pores_from_Bax_dimers(bax_conf='ins')
        self.model.name = 'bax_heat_auto1_dimer_reversible'

    def build_model_bax_heat_auto2_dimer_reversible(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates()
        self.Bax_dimerizes()
        self.Bax_reverses()
        self.pores_from_Bax_dimers(bax_conf='ins')
        self.model.name = 'bax_heat_auto2_dimer_reversible'

    # Heat models in which the heating doesn't make Bax pore-competent
    def build_model_bax_heat_bh3_exposure_auto2(self):
        self.translocate_Bax()
        self.basal_bh3_exposure_auto2()
        self.pores_from_Bax_monomers(bax_conf='ins')
        self.model.name = 'bax_heat_bh3_exposure_auto2'

    # Models incorporating dye release
    def build_model_bax_schwarz(self):
        self.translocate_Bax()
        self.pores_from_Bax_schwarz(bax_conf='mem')
        self.model.name = 'bax_schwarz'

    def build_model_bax_schwarz_irreversible(self):
        self.translocate_Bax()
        self.pores_from_Bax_monomers(bax_conf='mem')
        self.model.name = 'bax_schwarz_irreversible'

    def build_model_bax_schwarz_irreversible_aggregation(self):
        self.translocate_Bax()
        self.pores_from_Bax_monomers(bax_conf='mem')
        self.pores_aggregate(bax_conf='mem')
        self.model.name = 'bax_schwarz_irreversible_aggregation'

    def build_model_bax_schwarz_dimer(self):
        self.translocate_Bax()
        self.Bax_dimerizes(bax_loc_state='m')
        self.pores_from_Bax_dimers(bax_loc_state='m', reversible=False)
        self.model.name = 'bax_schwarz_dimer'

    def build_model_bax_schwarz_dimer_reversible(self):
        self.translocate_Bax()
        self.Bax_dimerizes(bax_loc_state='m')
        self.pores_from_Bax_dimers(bax_loc_state='m', reversible=True)
        self.model.name = 'bax_schwarz_dimer_reversible'

    def build_model_bax_schwarz_tetramer_reversible(self):
        self.translocate_Bax()
        self.Bax_dimerizes(bax_loc_state='m')
        self.Bax_tetramerizes(bax_loc_state='m')
        self.pores_from_Bax_tetramers(bax_loc_state='m', reversible=True)
        self.model.name = 'bax_schwarz_tetramer_reversible'

    def build_model_tap1(self):
        if LOGGING:
            print "---------------------------"
            print "core: Building model tap1:"
        self.build_model_ta()
        self.pores_from_Bax_monomers()
        self.model.name = 'tap1'

    def build_model_peptide_solution_dimer(self):
        # First, declare the model topology
        self.translocate_Bax()
        self.Bax_dimerizes(bax_loc_state='c')
        self.translocate_Bax_dimers()
        self.Bax_dimerizes(bax_loc_state='m')
        self.pores_from_Bax_dimers(bax_loc_state='m', reversible=False)

        # Create the expressions for the equilibrated concentrations
        # of solution monomer and dimer
        # First, some components we'll need
        Bax_0 = self['Bax_0']
        kf = self['Bax_c_dimerization_kf']
        kr = self['Bax_c_dimerization_kr']
        ka = self.expression('Bax_transloc_KA', (0.5*kf) / kr)
        Bax = self['Bax']
        solution = self['solution']

        # Concentration of solution Bax dimer
        Bax2_0 = self.expression('Bax2_0',
                             (1/ka)*(0.5*Bax_0*ka -
                                     0.125*(8*Bax_0*ka + 1)**0.5 + 0.125))
        # Concentration of solution Bax monomer
        Bax1_0 = self.expression('Bax1_0', Bax_0 - 2*Bax2_0)

        # Reset initial conditions accordingly
        self.model.initial_conditions = []
        bax_args = {'a6':None, 'lipo':None, 'c3':'s', 'c62':'s', 'c120':'s',
                'c122':'s', 'c126':'s', 'c184':'s'}
        self.initial(Bax(bh3=1, loc='c', **bax_args) %
                  Bax(bh3=1, loc='c', **bax_args), Bax2_0)
        self.initial(Bax(bh3=None, loc='c', **bax_args)**solution, Bax1_0)
        self.initial(self['Vesicles'](bax=None, tbid=None) **solution,
                     self['Vesicles_norm_0'])



