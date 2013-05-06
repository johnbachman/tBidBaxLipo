"""
Basic implementation of interactions and mechanisms for tBid/Bax interactions
with each other and with membranes. Contains macros to build up a variety
of alternative models of tBid/Bax interactions.

Code organization
=================

Models of tBid/Bax interactions draw on the same pool of mechanistic facts, but
may differ in their implementation of how interactions with membranes are
modeled. As a result, some rules (e.g., those that specify how proteins
translocate to membranes) have to be written in a compartmen-specific way,
whereas others can be written in a generic way. To manage this degree of
complexity, macros implementing tBid/Bax interactions are implemented here in a
generic way (e.g.,
:py:meth:`tbidbaxlipo.models.core.Builder.tBid_activates_Bax`); if
compartment-specific modifications are necessary, these functions are overriden
in child classes (e.g.,
:py:meth:`tbidbaxlipo.models.site_cpt.Builder.tBid_activates_Bax`).

The model builder classes contained in this file and also in
:py:mod:`tbidbaxlipo.models.one_cpt`, :py:mod:`tbidbaxlipo.models.n_cpt`,
and :py:mod:`tbidbaxlipo.models.site_cpt`, are used to manage the process of
building up alternative models. Each of these modules contains a class
`Builder` that contains an instance of a PySB model. The classes also contain a
series of methods implementing small sub-pieces of mechanism that can be termed
"motifs". These motifs can be recombined in different ways to create different
models.

The pattern for model construction used here does not rely on the SelfExporter
class of PySB. Instead, the ``Builder`` class contains an instance of a PySB
model object. Monomers, Parameters, Rules, etc. are added to this model object
by invoking the wrapper functions included in
:py:mod:`tbidbaxlipo.model.core`. These include

- :py:meth:`tbidbaxlipo.model.core.Builder.monomer`
- :py:meth:`tbidbaxlipo.model.core.Builder.parameter`
- :py:meth:`tbidbaxlipo.model.core.Builder.rule`
- :py:meth:`tbidbaxlipo.model.core.Builder.compartment`
- :py:meth:`tbidbaxlipo.model.core.Builder.initial`
- :py:meth:`tbidbaxlipo.model.core.Builder.observable`

Each of these functions invokes the appropriate PySB component constructor
while setting the argument ``_export=False`` so that the SelfExporter is not
invoked. The created components are then added to the instance of the model
contained in the builder class.

In addition, the model builder class implements the ``__getitem__`` method so
that invoking ``self['component_name']`` from within any method returns the
component with the given name.

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

To-do list
==========

.. todo:: List outline of the mechanistic alternatives before implementing them.

.. todo:: Think about how models could be built from rules or macros using tables.

.. todo:: Make it so macros can be invoked from within the motifs in this

paradigm.

"""

__author__ = 'johnbachman'

from pysb import *
from pysb.macros import *
import numpy as np

class Builder(object):

    # -- VIRTUAL FUNCTIONS ----------------------------------------------------
    def within_compartment_rsf(self):
        raise NotImplementedError()

    def translocate_tBid_Bax(self):
        raise NotImplementedError()

    def run_model(self):
        raise NotImplementedError()

    # -- CONSTRUCTOR AND MONOMER DECLARATIONS --------------------------------
    def __init__(self, params_dict=None):
        """Base constructor for all model builder classes.

        Initializes collections of the parameters to estimate, as well
        as the means and variances of their priors.

        Parameters
        ----------
        params_dict : dict
            The params_dict allows any parameter value to be overriden
            by name; any parameters not included in the dict will be set
            to default values. For example, if params_dict contains::

                {'tBid_Bax_kf': 1e-2}

            then the parameter tBid_Bax_kf will be assigned a value of 1e-2;
            all other parameters will take on default values. However,
            note that the parameter value given will be multiplied by any
            scaling factor passed in when the parameter is declared.
        """

        self.model = Model('tBid_Bax', _export=False)
        self.estimate_params = []
        self.parameter_means = np.array([])
        self.parameter_variances = np.array([])
        self.params_dict = params_dict

    def declare_monomers(self):
        """Declares signatures for tBid and Bax."""
        self.monomer('tBid', ['bh3', 'loc'],
                {'loc': ['c', 'm']})
        self.monomer('Bax',
                        ['bh3', 'a6', 'loc',
                         'c3', 'c62', 'c120', 'c122', 'c126', 'c184'],
                        {'loc':  ['c', 'm', 'i', 'a', 'p'],
                         'c3':   ['s', 'm'],
                         'c62':  ['s', 'm'],
                         'c120': ['s', 'm'],
                         'c122': ['s', 'm'],
                         'c126': ['s', 'm'],
                         'c184': ['s', 'm']})
        self.monomer('Vesicles', [])

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
        if (position[0] < -0.1 or position[0] > 0.7):
            return np.inf
        if any(position[1:] > 0) or any(position[1:] < -5):
            return np.inf
        else:
            return np.sum((position - self.parameter_means)**2 / \
                          (2 * self.parameter_variances))

    def random_initial_values(self):
        """Return a set of random initial values for parameter estimation.

        Uses the lognormal prior distributions specified for the parameters to
        be estimated. Currently it uses the specified variance, though it may
        be worth considering using a larger variance to ensure that the
        starting positions are "overdispersed."

        Note that this must be called after the model has been built using one
        of the model building methods--otherwise the parameters won't have been
        created!

        Returns
        -------
        An array of initial parameter values in the same order as the
        parameters in self.estimate_params, and drawn from lognormal (base 10)
        distributions with the means given by self.parameter_means and the
        variances given by self.parameter_variances.
        """
        while (True):
            test_vals = (self.parameter_means +
                            (np.random.randn(len(self.estimate_params)) * \
                             self.parameter_variances))

            if any(test_vals[1:] > 0) or any(test_vals[1:] < -5):
                continue
            elif test_vals[0] < -0.1 or test_vals[0] > 0.7:
                continue
            else:
                return 10 ** test_vals

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

        #bind(tBid(loc='m'), 'bh3', Bax(loc='m'), bax_site,
        #  [tBid_mBax_kf, tBid_mBax_kr])
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
        """Reversion of the inserted form of Bax to the loosely associated
        state, at which point it can return to the solution."""

        print('core: Bax_reverses')

        Bax = self['Bax']

        # Reversion of active Bax (P -> S)
        krev = self.parameter('iBax_reverse_k', 1e-2)

        # iBax reverses back to mBax
        self.rule('iBax_reverses',
             Bax(loc='i', bh3=None, a6=None) >> Bax(loc='m', bh3=None, a6=None),
             krev)

    def Bax_dimerizes(self):
        """
        Notes
        -----
        - If the late phase of the 62c signal is an indication of dimerization,
          it starts to manifest around 500s.
        - In SATSOURA, Fig 4. appears to indicate that the Bax-Bax FRET reaches
          steady-state at around 12 minutes.
        """

        print("core: Bax_dimerizes")
     
        # Rate of dimerization formation/oligomerization of activated Bax (s^-1)
        Bax_dimerization_kf = self.parameter('Bax_dimerization_kf', 1e-2)# was 1
        Bax_dimerization_kr = self.parameter('Bax_dimerization_kr', 1e-2)

        Bax = self['Bax']

        self.rule('Bax_Forms_Dimers',
             Bax(loc='i', bh3=None, a6=None) +
             Bax(loc='i', bh3=None, a6=None) <>
             Bax(loc='i', bh3=1, a6=None) % Bax(loc='i', bh3=1, a6=None),
             Bax_dimerization_kf, Bax_dimerization_kr)

    def Bax_tetramerizes(self):
        """
        This function depends on Bax_dimerization to be called as well.

        Notes
        -----
        In Lovell Fig S1, about 80% of the Bax is at membranes (inserted,
        non-extractable, resistant to gel filtration etc.) after       
        """

        print("core: Bax_tetramerizes()")

        # Rate of dimerization formation/oligomerization of activated Bax (s^-1)
        Bax_tetramerization_kf = self.parameter('Bax_tetramerization_kf', 1e-2)
        Bax_tetramerization_kr = self.parameter('Bax_tetramerization_kr', 1e-2) 

        Bax = self['Bax']

        self.rule('Bax_Forms_Tetramers',
             MatchOnce(Bax(loc='i', bh3=1, a6=None) %
                       Bax(loc='i', bh3=1, a6=None)) +
             MatchOnce(Bax(loc='i', bh3=2, a6=None) %
                       Bax(loc='i', bh3=2, a6=None)) <>
             Bax(loc='i', bh3=1, a6=3) % Bax(loc='i', bh3=1, a6=4) % 
             Bax(loc='i', bh3=2, a6=3) % Bax(loc='i', bh3=2, a6=4), 
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
    def Bax_auto_activates(target_bax_site='a6'):
        # Andrews suggests that tBid/Bax Kd should work out to 25nM
        # Forward rate of iBax binding to Bax (E + S -> ES)
        Parameter('iBax_mBax_kf', 1e-4)
        # Reverse rate of iBax binding to Bax (ES -> E + S)
        Parameter('iBax_mBax_kr', 1e-2)
        Parameter('mBaxiBax_to_iBaxiBax_k', 1) 
        # Dissociation of iBax from iBax (EP -> E + P)
        #Parameter('iBax_iBax_kr', 2.5e-3)

        bind(Bax(loc='i'), 'bh3', Bax(loc='m'), target_bax_site, # FULL
          [iBax_mBax_kf, iBax_mBax_kr])
        bind(Bax(loc='p'), 'bh3', Bax(loc='m'), target_bax_site, # FULL
          [iBax_mBax_kf, iBax_mBax_kr])

        # Create the dicts to parameterize the site that iBax binds to
        target_bax_site_bound = {target_bax_site:1}
        target_bax_site_unbound = {target_bax_site:None}

        # Conformational change of Bax (ES -> EP)
        Rule('mBaxiBax_to_iBaxiBax',
             Bax(loc='i', bh3=1) % Bax(loc='m', **target_bax_site_bound) >>
             Bax(loc='i', bh3=None) + Bax(loc='i', **target_bax_site_unbound),
             mBaxiBax_to_iBaxiBax_k)

    def Bax_aggregates_at_pores():
        Parameter('aggregation_rate_k', 1e-4)
        Rule('Bax_aggregates_at_pores',
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

    def basal_Bax_activation():
        # Spontaneous rate of transition of Bax from the mitochondrial to the
        # inserted state
        # Implies average time is 10000 seconds???
        Parameter('basal_Bax_kf', 1e-3)
        Parameter('basal_Bax_kr', 10)

        two_state_equilibrium(Bax(bh3=None, dye='f'), 'm', 'i',
                [basal_Bax_kf, basal_Bax_kr], sitename='loc')
        two_state_equilibrium(Bax(bh3=None, dye='e'), 'm', 'i',
                [basal_Bax_kf, basal_Bax_kr], sitename='loc')

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

    # -- CONSTRUCTOR WRAPPER FUNCTIONS ---------------------------------------
    def monomer(self, *args, **kwargs):
        """Adds a parameter to the Builder's model instance."""
        m = Monomer(*args, _export=False, **kwargs)
        self.model.add_component(m)
        return m

    def parameter(self, name, value, factor=1, estimate=True,
                  mean=-3.0, variance=2.0):
        """Adds a parameter to the Builder's model instance.

        Examines the params_dict attribute of the Builder instance (which is
        set in the constructor,
        :py:meth:`tbidbaxlipo.models.core.Builder.__init__`).  If the
        parameter with the given name is in the ``params_dict``, then the value
        in the ``params_dict`` is used to construct the parameter, and the
        argument ``value`` is ignored. If the parameter is not in
        ``params_dict``, then the parameter is assigned ``value``.

        Furthermore, in all cases the parameter value is multiplied by a
        scaling factor specified by the argument ``factor``. This allows rate
        scaling factor that are dependent on lipid concentrations or on units
        (deterministic vs. stochastic) to be handled by keeping the same
        parameter value but passing in the appropriate value for ``factor``.

        Parameters
        ----------
        name : string
            The name of the parameter to add
        value : number
            The value of the parameter
        factor : number
            A scaling factor to be applied to the parameter value.
        estimate : boolean
            Specifies whether the parameter should be included among the
            parameters to estimate, contained in the set
            ``Builder.estimate_params``. Defaults to True.
        mean : float
            The mean value of the log10 of the parameter's prior (lognormal)
            distribution. Defaults to -3. Ignored if ``estimate`` is False.
        variance : float
            The variance (in log10 space) of the parameter's prior (lognormal)
            distribution. Defaults to 2.0. Ignored if ``estimate`` is False.
        """

        if self.params_dict is None:
            param_val = value * factor
        else:
            if name in self.params_dict:
                param_val = self.params_dict[name] * factor
            else:
                param_val = value * factor

        p = Parameter(name, param_val, _export=False)
        self.model.add_component(p)

        if estimate:
            self.estimate_params.append(p)
            self.parameter_means = np.append(self.parameter_means, mean)
            self.parameter_variances = np.append(self.parameter_variances,
                                                 variance)

        return p

    def rule(self, *args, **kwargs):
        """Adds a rule to the Builder's model instance."""
        r = Rule(*args, _export=False, **kwargs)
        self.model.add_component(r)
        return r 

    def compartment(self, *args, **kwargs):
        """Adds a compartment to the Builder's model instance."""
        c = Compartment(*args, _export=False, **kwargs)
        self.model.add_component(c)
        return c 

    def observable(self, *args, **kwargs):
        """Adds an observable to the Builder's model instance."""
        o = Observable(*args, _export=False, **kwargs)
        self.model.add_component(o)
        return o

    def initial(self, *args):
        """Adds an initial condition to the Builder's model instance."""
        self.model.initial(*args)

    def __getitem__(self, index):
        """Returns the component with the given string index
        from the instance of the model contained by the Builder."""
        return self.model.all_components()[index]

