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
in child classes (e.g.,
:py:meth:`tbidbaxlipo.models.site_cpt.Builder.tBid_activates_Bax`).

The model builder classes contained in this file and also in
:py:mod:`tbidbaxlipo.models.one_cpt`, :py:mod:`tbidbaxlipo.models.n_cpt`, and
:py:mod:`tbidbaxlipo.models.site_cpt`, are used to manage the process of
building up alternative models. Each of these modules contains a class
`Builder` that contains an instance of a PySB model. The classes also contain a
series of methods implementing small sub-pieces of mechanism that can be termed
"motifs". These motifs can be recombined in different ways to create different
models.

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
  the BH3 exposed but inhibits further binding by tBid at that site. References:
  [Gavathiotis2008]_, [Kim2009]_.

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
              K1               K2                 K3

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
from bayessb.priors import Normal

class Builder(pysb.builder.Builder):

    # -- VIRTUAL FUNCTIONS ----------------------------------------------------
    def within_compartment_rsf(self):
        raise NotImplementedError()

    def translocate_tBid_Bax(self):
        raise NotImplementedError()

    def run_model(self):
        raise NotImplementedError()

    def declare_monomers(self):
        """Declares signatures for tBid and Bax."""
        self.monomer('tBid', ['bh3', 'loc'],
                {'loc': ['c', 'm']})
        self.monomer('Bax',
                        ['bh3', 'a6', 'loc', 'lipo',
                         'c3', 'c62', 'c120', 'c122', 'c126', 'c184'],
                        {'loc':  ['c', 'm', 'i', 'a', 'p', 'bh3expos'],
                         'c3':   ['s', 'm'],
                         'c62':  ['s', 'm'],
                         'c120': ['s', 'm'],
                         'c122': ['s', 'm'],
                         'c126': ['s', 'm'],
                         'c184': ['s', 'm']})
        self.monomer('Vesicles', ['bax'])
        self.monomer('Pores', [])

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
            self.parameter('c3_scaling', 0.8,
                           mean=np.log10(0.8), variance=0.04)
            site_set = True
        if 'c62' in nbd_sites:
            self.parameter('c62_scaling', 0.9204,
                   mean=np.log10(0.9204), variance=0.04)
            site_set = True
        if 'c120' in nbd_sites:
            self.parameter('c120_scaling', 0.975,
                   mean=np.log10(0.975), variance=0.04)
            site_set = True
        if 'c122' in nbd_sites:
            self.parameter('c122_scaling', 0.952,
                   mean=np.log10(0.952), variance=0.04)
            site_set = True
        if 'c126' in nbd_sites:
            self.parameter('c126_scaling', 0.966,
                   mean=np.log10(0.966), variance=0.04)
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
    def tBid_activates_Bax(self, bax_site='bh3'):
        """Default implementation of Bax activation by tBid.

        Takes two arguments:
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

        print('core: tBid_activates_Bax(bax_site=%s)' % bax_site)

        # Forward rate of tBid binding to Bax (E + S -> ES)
        kf = self.parameter('tBid_mBax_kf', 1e-2,
                            factor=self.within_compartment_rsf())
        # Reverse rate of tBid binding to Bax (ES -> E + S)
        kr = self.parameter('tBid_mBax_kr', 1.5)
        # Dissociation of tBid from iBax (EP -> E + P)
        kc = self.parameter('tBid_iBax_kc', 1e-1)

        # Create the dicts to parameterize the site that tBid binds to
        bax_site_bound = {bax_site:1}
        bax_site_unbound = {bax_site:None}

        tBid = self['tBid']
        Bax = self['Bax']

        self.rule('tBid_Bax_bind',
             tBid(loc='m', bh3=None) + Bax(loc='m', **bax_site_unbound) >>
             tBid(loc='m', bh3=1) % Bax(loc='m', **bax_site_bound),
             kf)
        self.rule('tBid_Bax_unbind',
             tBid(loc='m', bh3=1) % Bax(loc='m', **bax_site_bound) >>
             tBid(loc='m', bh3=None) + Bax(loc='m', **bax_site_unbound),
             kr)

        # tBid dissociates from iBax after activation
        self.rule('tBid_unbinds_iBax',
             tBid(bh3=1) % Bax(loc='m', **bax_site_bound) >>
             tBid(bh3=None) + Bax(loc='i', **bax_site_unbound),
             kc)

    def Bax_reverses(self):
        """Reversion of the inserted form of Bax to the solution."""

        print('core: Bax_reverses')

        Bax = self['Bax']
        sol = self['solution']
        ves = self['ves']

        # Reversion of active Bax (P -> S)
        krev = self.parameter('iBax_reverse_k', 1e-3)

        # iBax reverses back to mBax
        self.rule('iBax_reverses',
             Bax(loc='i', bh3=None, a6=None, lipo=ANY) ** ves >>
             Bax(loc='c', bh3=None, a6=None, lipo=None) ** sol,
             krev)

    def Bax_nmerizes(self, n, bax_loc_state='i'):
        Bax_nmerization_kf = self.parameter(
                'Bax_%s_%dmerization_kf' % (bax_loc_state, n),
                1e-2, factor=self.within_compartment_rsf())
        Bax_nmerization_kr = self.parameter(
                'Bax_%s_%dmerization_kr' % (bax_loc_state, n),
                1e-2, factor=self.within_compartment_rsf())
        Bax = self['Bax']

        free_Bax = Bax(loc=bax_loc_state, bh3=None, a6=None)
        lhs = free_Bax
        rhs = Bax(loc=bax_loc_state, bh3=0, a6=n-1)
        for i in range(n-1):
            lhs += free_Bax
            rhs = rhs % Bax(loc=bax_loc_state, bh3=i+1, a6=i)

        self.rule('Bax_%s_forms_%dmers' % (bax_loc_state, n),
             lhs <> rhs, Bax_nmerization_kf, Bax_nmerization_kr)
        self.observable('Bax%d_%s' % (n, bax_loc_state), rhs)

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

    def Bax_dimerizes(self, bax_loc_state='i'):
        """
        Notes
        -----
        - If the late phase of the 62c signal is an indication of dimerization,
          it starts to manifest around 500s.
        - In SATSOURA, Fig 4. appears to indicate that the Bax-Bax FRET reaches
          steady-state at around 12 minutes.
        """

        print("core: Bax_dimerizes(bax_loc_state=%s)" % bax_loc_state)

        # Rate of dimerization formation/oligomerization of activated Bax (s^-1)
        Bax_dimerization_kf = self.parameter(
               'Bax_%s_dimerization_kf' % bax_loc_state, 1e-2,
               factor=self.within_compartment_rsf())
        Bax_dimerization_kr = self.parameter(
               'Bax_%s_dimerization_kr' % bax_loc_state, 1e-2)

        Bax = self['Bax']

        self.rule('Bax_%s_forms_dimers' % bax_loc_state,
             Bax(loc=bax_loc_state, bh3=None, a6=None) +
             Bax(loc=bax_loc_state, bh3=None, a6=None) <>
             Bax(loc=bax_loc_state, bh3=1, a6=None) %
             Bax(loc=bax_loc_state, bh3=1, a6=None),
             Bax_dimerization_kf, Bax_dimerization_kr)
        self.observable('Bax2_%s' % bax_loc_state,
             Bax(loc=bax_loc_state, bh3=1, a6=None) %
             Bax(loc=bax_loc_state, bh3=1, a6=None))

    def Bax_tetramerizes(self, bax_loc_state='i'):
        """
        This function depends on Bax_dimerization to be called as well.

        Notes
        -----
        In Lovell Fig S1, about 80% of the Bax is at membranes (inserted,
        non-extractable, resistant to gel filtration etc.) after
        """

        print("core: Bax_tetramerizes()")

        # Rate of dimerization formation/oligomerization of activated Bax (s^-1)
        Bax_tetramerization_kf = self.parameter('Bax_tetramerization_kf', 1e-2,
               factor=self.within_compartment_rsf())
        Bax_tetramerization_kr = self.parameter('Bax_tetramerization_kr', 1e-2)

        Bax = self['Bax']

        self.rule('Bax_Forms_Tetramers',
             MatchOnce(Bax(loc=bax_loc_state, bh3=1, a6=None) %
                       Bax(loc=bax_loc_state, bh3=1, a6=None)) +
             MatchOnce(Bax(loc=bax_loc_state, bh3=2, a6=None) %
                       Bax(loc=bax_loc_state, bh3=2, a6=None)) <>
             Bax(loc=bax_loc_state, bh3=1, a6=3) %
             Bax(loc=bax_loc_state, bh3=1, a6=4) %
             Bax(loc=bax_loc_state, bh3=2, a6=3) %
             Bax(loc=bax_loc_state, bh3=2, a6=4), 
             Bax_tetramerization_kf, Bax_tetramerization_kr)

    def iBax_binds_tBid_at_bh3(self):
        """
        Notes
        -----
        THERE IS A PROBLEM WITH THIS!!!

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
        """
        print("core: iBax_binds_tBid_at_bh3()")

        # INHIBITION OF TBID BY BAX
        # Rate of tBid binding to iBax (E + P -> EP)
        kf = self.parameter('tBid_iBax_kf', 1e-1)
        # Rate of tBid binding to iBax (E + P -> EP)
        kr = self.parameter('tBid_iBax_kr', 2.5)

        tBid = self['tBid']
        Bax = self['Bax']

        # Binding between mtBid and iBax
        self.rule('tBid_iBax_bind_at_bh3',
             tBid(loc='m', bh3=None) + Bax(loc='i', bh3=None, a6=None) <>
             tBid(loc='m', bh3=1) % Bax(loc='i', bh3=1, a6=None),
             kf, kr)

    # -- OTHER MOTIFS --------------------------------------------------------
    def Bax_auto_activates(self, target_bax_site='a6'):
        print "core: Bax_auto_activates"
        # Andrews suggests that tBid/Bax Kd should work out to 25nM
        # Forward rate of iBax binding to Bax (E + S -> ES)
        iBax_mBax_kf = self.parameter('iBax_mBax_kf', 1e-3)
        # Reverse rate of iBax binding to Bax (ES -> E + S)
        iBax_mBax_kr = self.parameter('iBax_mBax_kr', 10)
        # Dissociation of iBax from iBax (EP -> E + P)
        iBaxmBax_to_iBax_iBax_k = self.parameter('iBaxmBax_to_iBax_iBax_k', 100) 

        Bax = self['Bax']

        self.rule('iBax_binds_mBax',
                Bax(loc='i', a6=None, bh3=None) +
                Bax(loc='m', a6=None, bh3=None) <>
                Bax(loc='i', a6=None, bh3=1) % Bax(loc='m', a6=1, bh3=None),
                iBax_mBax_kf, iBax_mBax_kr)
        self.rule('iBaxmBax_to_iBax_iBax',
                Bax(loc='i', a6=None, bh3=1) %
                Bax(loc='m', a6=1, bh3=None) >>
                Bax(loc='i', a6=None, bh3=None) +
                Bax(loc='i', a6=None, bh3=None),
                iBaxmBax_to_iBax_iBax_k)

    def Bax_auto_activates_one_step(self):
        print "core: Bax_auto_activates_one_step"
        # Andrews suggests that tBid/Bax Kd should work out to 25nM
        # Forward rate of iBax binding to Bax (E + S -> ES)
        iBax_activates_mBax_k = self.parameter('iBax_activates_mBax_k', 1e-4)
        # Reverse rate of iBax binding to Bax (ES -> E + S)
        #self.parameter('iBax_mBax_kr', 1e-2)
        #self.parameter('mBaxiBax_to_iBaxiBax_k', 1) 
        # Dissociation of iBax from iBax (EP -> E + P)
        #Parameter('iBax_iBax_kr', 2.5e-3)
        Bax = self['Bax']

        # Conformational change of Bax (ES -> EP)
        self.rule('iBax_activates_mBax',
            Bax(loc='i') + Bax(loc='m') >> Bax(loc='i') + Bax(loc='i'),
            iBax_activates_mBax_k)

    def Bax_aggregates_at_pores(self):
        Bax = self['Bax']
        aggregation_rate_k = self.parameter('aggregation_rate_k', 1e-4)
        self.rule('Bax_aggregates_at_pores',
             Bax(loc='p') + Bax(loc='m') >> Bax(loc='p') + Bax(loc='p'),
             aggregation_rate_k)

    def tBid_reverses_Bax():
        # Rate of the EP->ES transition # FIXME
        Parameter('iBaxtBid_to_mBaxtBid_k', 1e-3)

        # REVERSIBILITY OF BAX ACTIVATION BY TBID (EP -> ES)
        Rule('iBaxtBid_to_mBaxtBid',
             tBid(loc='m', bh3=1) % Bax(loc='i', bh3=1) >>
             tBid(loc='m', bh3=1) % Bax(loc='m', bh3=1),
             iBaxtBid_to_mBaxtBid_k)

    def basal_bh3_exposure_auto2(self):
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

    def basal_Bax_activation(self, reversible=False):
        print "core: basal_Bax_activation, reversible=%s" % reversible
        # Spontaneous rate of transition of Bax from the mitochondrial to the
        # inserted state
        # Implies average time is 10000 seconds???
        Bax = self['Bax']

        basal_Bax_kf = self.parameter('basal_Bax_kf', 2e-3,
                    prior=Normal(-3,1))
        self.rule('basal_Bax_activation',
                  Bax(bh3=None, loc='m') >> Bax(bh3=None, loc='i'),
                  basal_Bax_kf)

        if reversible:
            basal_Bax_kr = self.parameter('basal_Bax_kr', 1e-4,
                                          prior=Normal(-5, 1))
            self.rule('basal_Bax_activation_rev',
                      Bax(bh3=None, loc='i') >> Bax(bh3=None, loc='m'),
                      basal_Bax_kr)

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
    def build_model_t(self):
        print "---------------------------"
        print "core: Building model t:"
        self.translocate_tBid_Bax()
        self.model.name = 't'

    def build_model_ta(self):
        print "---------------------------"
        print "core: Building model ta:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.model.name = 'ta'

    def build_model_tai(self):
        print "---------------------------"
        print "core: Building model tai:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.iBax_binds_tBid_at_bh3()
        self.model.name = 'tai'

    def build_model_taid(self):
        print "---------------------------"
        print "core: Building model taid:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.iBax_binds_tBid_at_bh3()
        self.Bax_dimerizes()
        self.model.name = 'taid'

    def build_model_taidt(self):
        print "---------------------------"
        print "core: Building model taidt:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.iBax_binds_tBid_at_bh3()
        self.Bax_dimerizes()
        self.Bax_tetramerizes()
        self.model.name = 'taidt'

    def build_model_tair(self):
        print "---------------------------"
        print "core: Building model tair:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.iBax_binds_tBid_at_bh3()
        self.Bax_reverses()
        self.model.name = 'tair'

    def build_model_taird(self):
        print "---------------------------"
        print "core: Building model taird:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.iBax_binds_tBid_at_bh3()
        self.Bax_reverses()
        self.Bax_dimerizes()
        self.model.name = 'taird'

    def build_model_tairdt(self):
        print "---------------------------"
        print "core: Building model tairdt:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.iBax_binds_tBid_at_bh3()
        self.Bax_reverses()
        self.Bax_dimerizes()
        self.Bax_tetramerizes()
        self.model.name = 'tairdt'

    def build_model_tad(self):
        print "---------------------------"
        print "core: Building model tad:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.Bax_dimerizes()
        self.model.name = 'tad'

    def build_model_tadt(self):
        print "---------------------------"
        print "core: Building model tadt:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.Bax_dimerizes()
        self.Bax_tetramerizes()
        self.model.name = 'tadt'

    def build_model_tar(self):
        print "---------------------------"
        print "core: Building model tar:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.Bax_reverses()
        self.model.name = 'tar'

    def build_model_tard(self):
        print "---------------------------"
        print "core: Building model tard:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.Bax_reverses()
        self.Bax_dimerizes()
        self.model.name = 'tard'

    def build_model_tardt(self):
        print "---------------------------"
        print "core: Building model tardt:"
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='bh3')
        self.Bax_reverses()
        self.Bax_dimerizes()
        self.Bax_tetramerizes()
        self.model.name = 'tardt'

    # Model for NBD/Terbium data
    def build_model_nbd_terbium(self):
        self.build_model_bax_heat_reversible()
        # Additional parameters
        pore_scaling = self.parameter('pore_scaling', 100., estimate=True,
                                prior=Normal(2, 1))
        c0_scaling = self.parameter('c0_scaling', 1, estimate=False)
        c1_scaling = self.parameter('c1_scaling', 1.5,
                                prior=Normal(0, 0.5))
        c2_scaling = self.parameter('c2_scaling', 2,
                                prior=Normal(0, 0.5))
        # Get observables
        self.expression('ScaledPores', self['pores'] / pore_scaling)
        self.expression('NBD', (c0_scaling * self['cBax'] +
                               c1_scaling * self['iBax'] +
                               c2_scaling * self['pBax']) / self['Bax_0'])
        self.model.name = 'nbd_terbium'

    # Models incorporating dye release
    def build_model_bax_heat(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.pores_from_Bax_monomers(bax_loc_state='i', reversible=False)
        self.model.name = 'bax_heat'

    def build_model_bax_heat_aggregation(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.pores_from_Bax_monomers(bax_loc_state='i', reversible=False)
        self.Bax_aggregates_at_pores()
        self.model.name = 'bax_heat_aggregation'

    def build_model_bax_heat_reversible(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.pores_from_Bax_monomers(bax_loc_state='i', reversible=True)
        self.model.name = 'bax_heat_reversible'

    def build_model_bax_heat_reversible_aggregation(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_aggregates_at_pores()
        self.pores_from_Bax_monomers(bax_loc_state='i', reversible=True)
        self.model.name = 'bax_heat_reversible_aggregation'

    def build_model_bax_heat_dimer(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_dimerizes(bax_loc_state='i')
        self.pores_from_Bax_dimers(bax_loc_state='i', reversible=False)
        self.model.name = 'bax_heat_dimer'

    def build_model_bax_heat_dimer_reversible(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_dimerizes(bax_loc_state='i')
        self.pores_from_Bax_dimers(bax_loc_state='i', reversible=True)
        self.model.name = 'bax_heat_dimer_reversible'

    def build_model_bax_heat_auto1(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates_one_step()
        self.pores_from_Bax_monomers(bax_loc_state='i', reversible=False)
        self.model.name = 'bax_heat_auto1'

    def build_model_bax_heat_auto2(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates()
        self.pores_from_Bax_monomers(bax_loc_state='i', reversible=False)
        self.model.name = 'bax_heat_auto2'

    def build_model_bax_heat_bh3_auto2(self):
        self.translocate_Bax()
        self.bh3_exposure_auto2()
        self.pores_from_Bax_monomers(bax_loc_state='i', reversible=False)
        self.model.name = 'bax_heat_auto2'

    def build_model_bax_heat_auto1_reversible_activation(self):
        self.translocate_Bax()
        self.basal_Bax_activation(reversible=True)
        self.Bax_auto_activates_one_step()
        self.pores_from_Bax_monomers(bax_loc_state='i', reversible=False)
        self.model.name = 'bax_heat_auto1_reversible_activation'

    def build_model_bax_heat_auto2_reversible_activation(self):
        self.translocate_Bax()
        self.basal_Bax_activation(reversible=True)
        self.Bax_auto_activates()
        self.pores_from_Bax_monomers(bax_loc_state='i', reversible=False)
        self.model.name = 'bax_heat_auto2_reversible_activation'

    def build_model_bax_heat_auto1_reversible(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates_one_step()
        self.pores_from_Bax_monomers(bax_loc_state='i', reversible=True)
        self.model.name = 'bax_heat_auto1_reversible'

    def build_model_bax_heat_auto2_reversible(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates()
        self.pores_from_Bax_monomers(bax_loc_state='i', reversible=True)
        self.model.name = 'bax_heat_auto2_reversible'

    def build_model_bax_heat_auto1_dimer(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates_one_step()
        self.Bax_dimerizes()
        self.pores_from_Bax_dimers(bax_loc_state='i', reversible=False)
        self.model.name = 'bax_heat_auto1_dimer'

    def build_model_bax_heat_auto2_dimer(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates()
        self.Bax_dimerizes()
        self.pores_from_Bax_dimers(bax_loc_state='i', reversible=False)
        self.model.name = 'bax_heat_auto2_dimer'

    def build_model_bax_heat_auto1_dimer_reversible(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates_one_step()
        self.Bax_dimerizes()
        self.pores_from_Bax_dimers(bax_loc_state='i', reversible=True)
        self.model.name = 'bax_heat_auto1_dimer_reversible'

    def build_model_bax_heat_auto2_dimer_reversible(self):
        self.translocate_Bax()
        self.basal_Bax_activation()
        self.Bax_auto_activates()
        self.Bax_dimerizes()
        self.pores_from_Bax_dimers(bax_loc_state='i', reversible=True)
        self.model.name = 'bax_heat_auto2_dimer_reversible'

    # Heat models in which the heating doesn't make Bax pore-competent
    def build_model_bax_heat_bh3_exposure_auto2(self):
        self.translocate_Bax()
        self.basal_bh3_exposure_auto2()
        self.pores_from_Bax_monomers(bax_loc_state='i')
        self.model.name = 'bax_heat_bh3_exposure_auto2'

    # Models incorporating dye release
    def build_model_bax_schwarz(self):
        self.translocate_Bax()
        self.pores_from_Bax_monomers(bax_loc_state='m', reversible=False)
        self.model.name = 'bax_schwarz'

    def build_model_bax_schwarz_reversible(self):
        self.translocate_Bax()
        self.pores_from_Bax_monomers(bax_loc_state='m', reversible=True)
        self.model.name = 'bax_schwarz_reversible'

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
        self.initial(Bax(bh3=1, loc='c', **bax_args)**solution %
                  Bax(bh3=1, loc='c', **bax_args)**solution, Bax2_0)
        self.initial(Bax(bh3=None, loc='c', **bax_args)**solution, Bax1_0)
        self.initial(self['Vesicles'](bax=None) **solution, self['Vesicles_0'])

