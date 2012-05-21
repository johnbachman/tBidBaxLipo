__author__ = 'johnbachman'

from pysb import *
from pysb.core import SelfExporter
from pysb.macros import *
#import inspect

class tBid_Bax(object):

## VIRTUAL FUNCTIONS

    #{{{# within_compartment_rsf()
    def within_compartment_rsf(self):
        raise NotImplementedError()
    #}}}

    #{{{# translocate_tBid_Bax()
    def translocate_tBid_Bax():
        raise NotImplementedError()
    #}}}

    #{{{# run_model()
    def run_model():
        raise NotImplementedError()
    #}}} 

## DEFAULT IMPLEMENTATIONS
    #{{{# __init__    
    def __init__(self, params_dict=None):
        self.model = Model('tBid_Bax', _export=False)

        # The params_dict allows any parameter value to be overriden
        # by name; any parameters not included in the dict will be set
        # to default values
        self.params_dict = params_dict
    #}}}

    #{{{# declare_monomers()
    def declare_monomers(self):
        self.monomer('tBid', ['bh3', 'loc'],
                {'loc': ['c', 'm']})
        self.monomer('Bax', ['bh3', 'a6', 'loc'],
                {'loc': ['c','m', 'i', 'p']})
        #self.monomer('Vesicles', [], {})
        #self.monomer('Pores')

        #Monomer('Pore', [], {})
        #Monomer('Vesicles', ['dye'],
        #        {'dye': ['f', 'e']})
    #}}}

    #{{{# tBid_activates_Bax()
    def tBid_activates_Bax(self, bax_site='a6'):
        """Takes two arguments:
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

        # Since the surface area of each vesicle is the same,
        # the effect of changing vesicle concentration on the forward rate will be
        # by altering the expected concentration of tBid and Bax per vesicle.
        # If there is one compartment, 100 Bax and 100 tBid, then the rate will be scaled
        # to take into account 100 of each.
        # If there are 100 compartments, then the rate should be scaled to take into account
        # 1 of each. The scaling should therefore most likely be simple linear scaling
        # by dividing by the number of compartments, which is in the same units as the
        # protein concentration (e.g., nanomolar).

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
        print("tBid_Bax: tBid_activates_Bax(bax_site=" + bax_site + ")")

        # Forward rate of tBid binding to Bax (E + S -> ES)
        kf = self.parameter('tBid_mBax_kf', 1, factor=self.within_compartment_rsf())
        #tBid_mBax_kf = Parameter('tBid_mBax_kf', 0)
        #Parameter('tBid_mBax_kf', 0) # Forward rate of tBid binding to Bax (E + S -> ES)
        # Reverse rate of tBid binding to Bax (ES -> E + S)
        kr = self.parameter('tBid_mBax_kr', 1)
        # Dissociation of tBid from iBax (EP -> E + P)
        #tBid_iBax_kc = Parameter('tBid_iBax_kc', 10)

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
        #Rule('tBid_unbinds_iBax',
        #     tBid(bh3=1) % Bax(loc='m', **bax_site_bound) >>
        #     tBid(bh3=None) + Bax(loc='i', **bax_site_unbound),
        #     tBid_iBax_kc)

    #}}}

    #{{{# build_model0
    def build_model0(self):
        print "---------------------------"
        print "tBid_Bax: Building model 0:"

        self.translocate_tBid_Bax()
        self.tBid_activates_Bax(bax_site='a6')

        #dye_release(Bax(loc='i', bh3=None))
        #pores_from_Bax_monomers()
        #Bax_dimerizes()
        #dye_release(Pore())
        #dye_release(Bax(loc='i', bh3=1, a6=None) % Bax(loc='i', bh3=1, a6=None))
    #}}}


## OTHER FUNCS ###################################################
#{{{
    #{{{# Bax_auto_activates()
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
    #}}}

    #{{{# Bax_reverses()
    def Bax_reverses():
        Parameter('Bax_i_to_c', 1e-4)
        Rule('Bax_reverses',
             Bax(loc='i', bh3=None, a6=None) >> Bax(loc='c', bh3=None, a6=None),
             Bax_i_to_c)
        Rule('Bax_reverses_p',
             Bax(loc='p', bh3=None, a6=None) >> Bax(loc='c', bh3=None, a6=None),
             Bax_i_to_c)
    #}}}

    #{{{# Bax_dimerizes()
    def Bax_dimerizes():
        print("Bax_dimerizes()")
        """CONSTRAINTS:
         - If the late phase of the 62c signal is an indication of dimerization, it
           starts to manifest around 500s.
         - In SATSOURA, Fig 4. appears to indicate that the Bax-Bax FRET reaches
           steady-state at around 12 minutes.
        """
     
        # Rate of dimerization formation/oligomerization of activated Bax (s^-1). 
        Bax_dimerization_kf = Parameter('Bax_dimerization_kf', 1e-2) # was 1
        Bax_dimerization_kr = Parameter('Bax_dimerization_kr', 1e0)

    #    Rule('Bax_Forms_Dimers',
    #         Bax(loc='i', bh3=None, a6=None) + Bax(loc='i', bh3=None, a6=None) <>
    #         Bax(loc='i', bh3=1, a6=None) % Bax(loc='i', bh3=1, a6=None),
    #         Bax_dimerization_kf, Bax_dimerization_kr)
        Rule('Bax_Forms_Dimers',
             Bax(loc='i', bh3=None, a6=None) + Bax(loc='i', bh3=None, a6=None) <>
             Pore(),
             Bax_dimerization_kf, Bax_dimerization_kr)
    #}}}

    #{{{# Bax_tetramerizes
    def Bax_tetramerizes():
        """ CONSTRAINTS:
        In Lovell Fig S1, about 80% of the Bax is at membranes (inserted,
        non-extractable, resistant to gel filtration etc.) after       
        """

        """ This function depends on Bax_dimerization to be called as well."""
        # Rate of dimerization formation/oligomerization of activated Bax (s^-1). 
        Parameter('Bax_tetramerization_kf', 1e-2) # was 1
        Parameter('Bax_tetramerization_kr', 1e-2) 

        Rule('Bax_Forms_Tetramers',
             MatchOnce(Bax(loc='i', bh3=1, a6=None) % Bax(loc='i', bh3=1, a6=None)) +
             MatchOnce(Bax(loc='i', bh3=2, a6=None) % Bax(loc='i', bh3=2, a6=None)) <>
             Bax(loc='i', bh3=1, a6=3) % Bax(loc='i', bh3=1, a6=4) % 
             Bax(loc='i', bh3=2, a6=3) % Bax(loc='i', bh3=2, a6=4), 
             Bax_tetramerization_kf, Bax_tetramerization_kr)
        Rule('Bax_Forms_Tetramers_p',
             MatchOnce(Bax(loc='p', bh3=1, a6=None) % Bax(loc='p', bh3=1, a6=None)) +
             MatchOnce(Bax(loc='p', bh3=2, a6=None) % Bax(loc='p', bh3=2, a6=None)) <>
             Bax(loc='p', bh3=1, a6=3) % Bax(loc='p', bh3=1, a6=4) % 
             Bax(loc='p', bh3=2, a6=3) % Bax(loc='p', bh3=2, a6=4), 
             Bax_tetramerization_kf, Bax_tetramerization_kr)
    #}}}

    #{{{# Bax_aggregates_at_pores()  
    def Bax_aggregates_at_pores():
        Parameter('aggregation_rate_k', 1e-4)
        Rule('Bax_aggregates_at_pores',
             Bax(loc='p') + Bax(loc='m') >> Bax(loc='p') + Bax(loc='p'),
             aggregation_rate_k)
    #}}}#

    #{{{# Bax_inhibits_tBid()
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
    def Bax_inhibits_tBid():
        # INHIBITION OF TBID BY BAX
        Parameter('tBid_iBax_kf', 1e-1) # Rate of tBid binding to iBax (E + P -> EP)
        Parameter('tBid_iBax_kr', 1) # Rate of tBid binding to iBax (E + P -> EP)

        # Binding between mtBid and iBax (activation back-reaction--should be slow)
        Rule('tBid_binds_iBax_f', tBid(loc='m', bh3=None) + Bax(loc='i', bh3=None) <>
             tBid(loc='m', bh3=1) % Bax(loc='i', bh3=1),
             tBid_iBax_kf, tBid_iBax_kr)
        Rule('tBid_binds_pBax_f', tBid(loc='m', bh3=None) + Bax(loc='p', bh3=None) <>
             tBid(loc='m', bh3=1) % Bax(loc='p', bh3=1),
             tBid_iBax_kf, tBid_iBax_kr)
        #Rule('tBid_binds_iBax_e', tBid(loc='m', dye='e', bh3=None) + Bax(loc='i', dye='e', bh3=None) >>
        #     tBid(loc='m', dye='e', bh3=1) % Bax(loc='i', dye='e', bh3=1),
        #     tBid_iBax_kf)
    #}}}

    #{{{# tBid_reverses_Bax()
    def tBid_reverses_Bax():
        Parameter('iBaxtBid_to_mBaxtBid_k', 1e-3) # Rate of the EP->ES transition # FIXME

        # REVERSIBILITY OF BAX ACTIVATION BY TBID (EP -> ES)
        Rule('iBaxtBid_to_mBaxtBid',
             tBid(loc='m', bh3=1) % Bax(loc='i', bh3=1) >>
             tBid(loc='m', bh3=1) % Bax(loc='m', bh3=1),
             iBaxtBid_to_mBaxtBid_k)
    #}}}

    #{{{# basal_Bax_activation()
    def basal_Bax_activation():
        # Spontaneous rate of transition of Bax from the mitochondrial to the inserted state
        Parameter('basal_Bax_kf', 1e-3) # Implies average time is 10000 seconds???
        Parameter('basal_Bax_kr', 10)

        two_state_equilibrium(Bax(bh3=None, dye='f'), 'm', 'i',
                [basal_Bax_kf, basal_Bax_kr], sitename='loc')
        two_state_equilibrium(Bax(bh3=None, dye='e'), 'm', 'i',
                [basal_Bax_kf, basal_Bax_kr], sitename='loc')
    #}}}
#}}}

#{{{ self.set...
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


#}}}

