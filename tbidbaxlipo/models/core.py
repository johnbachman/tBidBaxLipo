"""
TODOs:

- Since the models that I'm most interested in are the two_cpt and the site_cpt
  ones, this should ideally be written in such a way that the models don't have
  to be re-implemented--perhaps flag the rules that need to be given
  compartment specific implementations by adding a parameter?

- I should list, in this document, an outline of the alternative to mechanism
  that I want to explore, before implementing them.

- Since many of the mechanistic alternatives involve differences in binding
  site, or alternative conformational changes, it might be more useful to
  approach this by building up a set of rules as objects and then
  programmatically composing models from these rules using tables.

tBid Binding to Bax
-------------------

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

The Models
----------

Model 0
    - tBid and Bax translocate to membranes.

Model 1. 
    - tBid and Bax translocate to membranes
    - tBid causes Bax insertion into membranes

Model 2-inh.
    - tBid and Bax translocate to membranes
    - tBid causes Bax insertion into membranes
    - tBid 

"""
__author__ = 'johnbachman'

from pysb import *
#from pysb.core import SelfExporter
from pysb.macros import *
#import inspect

#SelfExporter.do_export = False

class Builder(object):

    # VIRTUAL FUNCTIONS
    # =================

    def within_compartment_rsf(self):
        raise NotImplementedError()

    def translocate_tBid_Bax():
        raise NotImplementedError()

    def run_model():
        raise NotImplementedError()

    # DEFAULT IMPLEMENTATIONS
    # =======================

    def __init__(self, params_dict=None):
        """Base constructor for all model builder classes.

        Parameters
        ----------
        params_dict : dict
            The params_dict allows any parameter value to be overriden
            by name; any parameters not included in the dict will be set
            to default values. For example, if params_dict contains::

                {'tBid_Bax_kf': 1e-2}

            then the parameter tBid_Bax_kf will be assigned a value of 1e-2;
            all other parameters will take on default values.
        """

        self.model = Model('tBid_Bax', _export=False)
        self.params_dict = params_dict

    def declare_monomers(self):
        """Declares signatures for tBid and Bax."""
        self.monomer('tBid', ['bh3', 'loc'],
                {'loc': ['c', 'm']})
        self.monomer('Bax', ['bh3', 'a6', 'loc'],
                {'loc': ['c','m', 'i', 'p']})

    def tBid_activates_Bax(self, bax_site='a6'):
        """Default implementation of Bax activation by tBid.

        Takes two arguments:
            - bax_site specifies the name of the site on Bax to which
              the bh3 site on tBid binds.
            - vesicles_conc is a Parameter object containing the concentration
              (in the same units as the proteins) of liposomes in the system.
              This value is used to scale the association rate of Bax and tBid
              on the membrane.
        """
        # CONSTRAINTS:
        # - Andrews suggests that tBid/Bax Kd should work out to 25nM (once in membrane, presumably)
        # - Binding of tBid to Bax during the activation step should be transient
        # - The timescale of 3c exposure has a half-life of about 50s (ln 2 / 1.381e-2)
        #   --this could potentially correspond to the "iBax" form, though arguably
        #   this may be better indicated by the a5/6 insertion per the finding of Annis
        #   that Bax is multispanning prior to oligomerization
        # - Binding of the BH3 (presumably by tBid?) occurs with an initial rate of ... (check fit)
        # - When tBid is added, 50-80% of Bax binds to liposomes, though this goes down at high
        #   Bax/liposome ratios. Lovell fig S1 suggests that this reaches steady state by about
        #   15 minutes.
        # FIXME The forward rate should be normalized according to protein/lipid ratio in some way

        # Since the surface area of each vesicle is the same, the effect of
        # changing vesicle concentration on the forward rate will be by
        # altering the expected concentration of tBid and Bax per vesicle.  If
        # there is one compartment, 100 Bax and 100 tBid, then the rate will be
        # scaled to take into account 100 of each.  If there are 100
        # compartments, then the rate should be scaled to take into account 1
        # of each. The scaling should therefore most likely be simple linear
        # scaling by dividing by the number of compartments, which is in the
        # same units as the protein concentration (e.g., nanomolar).

        # In the deterministic case, if the fundamental forward rate of binding
        # between tBid and Bax is kf, this should really be normalized by the
        # P/L ratio of both proteins. So for example,
        # kf * tBid/ves * Bax/ves
        # This because doubling of the vesicle concentration cuts the relative
        # concentrations of both proteins by half, and hence scales the forward
        # rate correspondingly. 
        # In the SSA case, 
        # The forward rate should be kf*RSF * tBid * Bax (where the concentrations
        # given are for that compartment). Since they represent concentrations
        # On that individual compartment, the rate does not need to be normalized
        # by the vesicle concentration.
        print('core: tBid_activates_Bax(bax_site=%s)' % bax_site)

        # Forward rate of tBid binding to Bax (E + S -> ES)
        kf = self.parameter('tBid_mBax_kf', 1e-2, factor=self.within_compartment_rsf())
        # Reverse rate of tBid binding to Bax (ES -> E + S)
        kr = self.parameter('tBid_mBax_kr', 1.5)
        # Dissociation of tBid from iBax (EP -> E + P)
        kc = self.parameter('tBid_iBax_kc', 1e-1)
        # Reversion of active Bax (P -> S)
        krev = self.parameter('iBax_reverse_k', 1e-2)

        # Create the dicts to parameterize the site that tBid binds to
        bax_site_bound = {bax_site:1}
        bax_site_unbound = {bax_site:None}

        #bind(tBid(loc='m'), 'bh3', Bax(loc='m'), bax_site,
        #  [tBid_mBax_kf, tBid_mBax_kr])
        tBid = self['tBid']
        Bax = self['Bax']
        self.rule('tBid_Bax_bind',
             tBid(loc='m', bh3=None) + Bax(loc='m', bh3=None, a6=None) >>
             tBid(loc='m', bh3=1) % Bax(loc='m', bh3=None, a6=1),
             kf)
        self.rule('tBid_Bax_unbind',
             tBid(loc='m', bh3=1) % Bax(loc='m', bh3=None, a6=1) >>
             tBid(loc='m', bh3=None) + Bax(loc='m', bh3=None, a6=None),
             kr)

        # tBid dissociates from iBax
        self.rule('tBid_unbinds_iBax',
             tBid(bh3=1) % Bax(loc='m', **bax_site_bound) >>
             tBid(bh3=None) + Bax(loc='i', **bax_site_unbound),
             kc)

        # iBax reverses back to mBax
        self.rule('iBax_reverses',
             Bax(loc='i', bh3=None, a6=None) >> Bax(loc='m', bh3=None, a6=None),
             krev)

        # Activation
        self.observable('tBidBax', tBid(loc='m', bh3=1) % Bax(loc='m', a6=1))
        self.observable('iBax', Bax(loc='i', bh3=None, a6=None))

    def Bax_dimerizes(self):
        """CONSTRAINTS:
         - If the late phase of the 62c signal is an indication of dimerization, it
           starts to manifest around 500s.
         - In SATSOURA, Fig 4. appears to indicate that the Bax-Bax FRET reaches
           steady-state at around 12 minutes.
        """
        print("tBid_Bax: Bax_dimerizes()")
     
        # Rate of dimerization formation/oligomerization of activated Bax (s^-1). 
        Bax_dimerization_kf = self.parameter('Bax_dimerization_kf', 1e-2) # was 1
        Bax_dimerization_kr = self.parameter('Bax_dimerization_kr', 1e0)

        Bax = self['Bax']

        self.rule('Bax_Forms_Dimers',
             Bax(loc='i', bh3=None, a6=None) + Bax(loc='i', bh3=None, a6=None) <>
             Bax(loc='i', bh3=1, a6=None) % Bax(loc='i', bh3=1, a6=None),
             Bax_dimerization_kf, Bax_dimerization_kr)

        self.observable('Bax2', 
             MatchOnce(Bax(loc='i', bh3=1, a6=None) % Bax(loc='i', bh3=1, a6=None)))

    def Bax_tetramerizes(self):
        """ CONSTRAINTS:
        In Lovell Fig S1, about 80% of the Bax is at membranes (inserted,
        non-extractable, resistant to gel filtration etc.) after       
        """
        print("tBid_Bax: Bax_tetramerizes()")

        """ This function depends on Bax_dimerization to be called as well."""
        # Rate of dimerization formation/oligomerization of activated Bax (s^-1). 
        Bax_tetramerization_kf = self.parameter('Bax_tetramerization_kf', 1e-2) # was 1
        Bax_tetramerization_kr = self.parameter('Bax_tetramerization_kr', 1e-2) 

        Bax = self['Bax']

        self.rule('Bax_Forms_Tetramers',
             MatchOnce(Bax(loc='i', bh3=1, a6=None) % Bax(loc='i', bh3=1, a6=None)) +
             MatchOnce(Bax(loc='i', bh3=2, a6=None) % Bax(loc='i', bh3=2, a6=None)) <>
             Bax(loc='i', bh3=1, a6=3) % Bax(loc='i', bh3=1, a6=4) % 
             Bax(loc='i', bh3=2, a6=3) % Bax(loc='i', bh3=2, a6=4), 
             Bax_tetramerization_kf, Bax_tetramerization_kr)

        self.observable('Bax4',
             MatchOnce(Bax(loc='i', bh3=1, a6=3) % Bax(loc='i', bh3=1, a6=4) % 
             Bax(loc='i', bh3=2, a6=3) % Bax(loc='i', bh3=2, a6=4)))

    # FIXME
    """
    THERE IS A PROBLEM WITH THIS!!!
    When tBid is bound to Bax, it is prevented from recirculating back to the solution.
    Therefore you get a sequence of events like:
    tBid + Bax (full liposome) -> tBid + Bax(i) (full liposome)
    Bax(p) (full) -> Bax(p) (empty) THIS HAPPENS VERY FAST
    But also: Bax(i) + tBid (full liposome). When complexed in this way,
    tBid does not recycle back to the solution. Therefore tBid catalyses the creation
    of a species Bax(i) which over time shifts the equilibrium of tBid from c to full liposomes.
    """
    def iBax_binds_tBid():
        print("tBid_Bax: iBax_binds_tBid()")

        # INHIBITION OF TBID BY BAX
        kf = self.parameter('tBid_iBax_kf', 1e-1) # Rate of tBid binding to iBax (E + P -> EP)
        kr = self.parameter('tBid_iBax_kr', 2.5) # Rate of tBid binding to iBax (E + P -> EP)

        tBid = self['tBid']
        Bax = self['Bax']

        # Binding between mtBid and iBax (activation back-reaction--should be slow)
        self.rule('tBid_iBax_bind',
             tBid(loc='m', bh3=None) + Bax(loc='i', bh3=None, a6=None) <>
             tBid(loc='m', bh3=1) % Bax(loc='i', bh3=1, a6=None),
             kf, kr)

## MODEL BUILDING FUNCTIONS
    def build_model0(self):
        print "---------------------------"
        print "core: Building model 0:"

        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='a6')

    def build_model1(self):
        """Activation, dimerization."""
        print "---------------------------"
        print "core: Building model 1:"

        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='a6')
        self.Bax_dimerizes()

    def build_model2(self):
        """Bax tetramerization."""
        print "---------------------------"
        print "core: Building model 2:"

        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='a6')
        self.Bax_dimerizes()
        self.Bax_tetramerizes()

    def build_model3(self):
        """iBax-tBid binding, tetramerization."""
        print "---------------------------"
        print "core: Building model 3:"

        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='a6')
        self.iBax_binds_tBid()
        self.Bax_dimerizes()
        self.Bax_tetramerizes()

## OTHER FUNCS ###################################################

    def Bax_auto_activates(target_bax_site='a6'):
        # Andrews suggests that tBid/Bax Kd should work out to 25nM
        Parameter('iBax_mBax_kf', 1e-4) # Forward rate of iBax binding to Bax (E + S -> ES)
        Parameter('iBax_mBax_kr', 1e-2)   # Reverse rate of iBax binding to Bax (ES -> E + S)
        Parameter('mBaxiBax_to_iBaxiBax_k', 1) 
        #Parameter('iBax_iBax_kr', 2.5e-3)   # Dissociation of iBax from iBax (EP -> E + P)

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

    def Bax_reverses():
        Parameter('Bax_i_to_c', 1e-4)
        Rule('Bax_reverses',
             Bax(loc='i', bh3=None, a6=None) >> Bax(loc='c', bh3=None, a6=None),
             Bax_i_to_c)
        Rule('Bax_reverses_p',
             Bax(loc='p', bh3=None, a6=None) >> Bax(loc='c', bh3=None, a6=None),
             Bax_i_to_c)

    def Bax_aggregates_at_pores():
        Parameter('aggregation_rate_k', 1e-4)
        Rule('Bax_aggregates_at_pores',
             Bax(loc='p') + Bax(loc='m') >> Bax(loc='p') + Bax(loc='p'),
             aggregation_rate_k)

    def tBid_reverses_Bax():
        Parameter('iBaxtBid_to_mBaxtBid_k', 1e-3) # Rate of the EP->ES transition # FIXME

        # REVERSIBILITY OF BAX ACTIVATION BY TBID (EP -> ES)
        Rule('iBaxtBid_to_mBaxtBid',
             tBid(loc='m', bh3=1) % Bax(loc='i', bh3=1) >>
             tBid(loc='m', bh3=1) % Bax(loc='m', bh3=1),
             iBaxtBid_to_mBaxtBid_k)

    def basal_Bax_activation():
        # Spontaneous rate of transition of Bax from the mitochondrial to the inserted state
        Parameter('basal_Bax_kf', 1e-3) # Implies average time is 10000 seconds???
        Parameter('basal_Bax_kr', 10)

        two_state_equilibrium(Bax(bh3=None, dye='f'), 'm', 'i',
                [basal_Bax_kf, basal_Bax_kr], sitename='loc')
        two_state_equilibrium(Bax(bh3=None, dye='e'), 'm', 'i',
                [basal_Bax_kf, basal_Bax_kr], sitename='loc')

    def monomer(self, *args, **kwargs):
        m = Monomer(*args, _export=False, **kwargs)
        self.model.add_component(m)
        return m

    def parameter(self, name, value, factor=1):
        if self.params_dict is None:
            param_val = value * factor
        else:
            if name in self.params_dict:
                param_val = self.params_dict[name] * factor
            else:
                param_val = value * factor

        p = Parameter(name, param_val, _export=False)
        self.model.add_component(p)
        return p

    def rule(self, *args, **kwargs):
        r = Rule(*args, _export=False, **kwargs)
        self.model.add_component(r)
        return r 

    def compartment(self, *args, **kwargs):
        c = Compartment(*args, _export=False, **kwargs)
        self.model.add_component(c)
        return c 

    def observable(self, *args, **kwargs):
        o = Observable(*args, _export=False, **kwargs)
        self.model.add_component(o)
        return o

    # Note: initial is not a component
    def initial(self, *args):
        self.model.initial(*args)

    def __getitem__(self, index):
        return self.model.all_components()[index]

