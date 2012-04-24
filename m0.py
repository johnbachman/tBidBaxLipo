__author__ = 'johnbachman'

#{{{# IMPORTS
from pysb import *
from pysb.macros import *
from pylab import *
from pysb.integrate import odesolve
from util.fitting import fit, fit_initial, mse
import util.fitting as fitting
#}}}

#{{{#  COMMENTS
"""
Simple model of tBid-Bax interaction, consisting only of tBid's activation of
Bax. The basic enzyme-substrate model.

Requirements of the model (these should be incorporated into tests):
- partitioning of tBid into membranes should occur rapidly and lead to nearly all
  (> 90-95%) of tBid in the membrane at equilibrium
- 

Things that are unknown:
- Does the partitioning of tBid and/or Bax to the membrane saturate,
  or is the leveling off due to non-ideal partitioning or aggregation?

The way this is written,

- Bax equilibrates equally between the cytosolic and peripherally bound state to both full and
empty liposomes.

- But, once in the "inserted" state, it can't dissociate from the vesicle!!
  (i.e., there's no reversal either to the cytosolic or peripherally bound states).

- Further, there's no reversal from the porated back to the inserted state--instead Bax
  goes from porated to the peripherally bound state on an empty vesicle.

There is a tradeoff between k_pore_rev and k_eflx. If k_pore_rev is slow, this means that pores,
once formed, are very stable, and so the trend is for almost all Bax to end up as oligomerized. Once
pore-state Bax reaches steady state, dye release occurs basically linearly afterwards.

Also note that in this model, once Bax is in the inserted state, it cannot revert back to the peripherally
bound state without first passing through the pore state. Also note that it cannot dissociate from the 
membrane unless it is in the 'm' (peripherally bound) state.

There is a cycle in that Bax can go from m to m via
mtBid + mBax <-> mtBid:mBax (<)-> mtBid + iBax (<)-> pBax
              K1               K2                 K3

Therefore we need to have
tBid_Bax_kf * tBid_Bax_kc * k_pore = tBid_Bax_kr * (tBid_iBax_kf) * k_pore_rev

The product of the rate constants in both directions should be the same--even though dye is released in this process,
one would imagine that the process would proceed exactly the same even without dye.

If the whole chain is at equilibrium, then each individual link must be at equilibrium.

Behavior of the model:---------
- Why does mftBid not come down as the number of empty vesicles goes up? Theoretically,
as tBid cycles back and forth to the cytoplasm, it should unbind, return to dye='none',
and then be dye='e' after rebinding an empty vesicle.


In the Almeida model, the pore is presumed to have a particular average lifetime,
and the transition from P_lf to P* is determined by some first order rate,
and then 
The rate of the efflux reaction is proportional to the amount of pBax(full),
the lifetime of the transition between pBax(full) -> pBax(empty), along with the
efflux rate constant.


So imagine that there are precisely 100 (V) vesicles, and 100 Bax molecules. Now
suppose that at a given time there is precisely 1 Bax molecule in the pore state
on a full vesicle, and this is sufficient to trigger dye release. Further
suppose that pore formation is irreversible. Within that timestep (say, 1 second),
you would expect that the pBax(full) concentration would go to 0, and the CF
(efflux) function would change by precisely 1%.

"""
#}}}

model = Model('tBidBax')

#{{{# MONOMERS
Monomer('tBid', ['bh3', 'loc', 'dye'],
        {'loc': ['c', 'm'],
         'dye': ['none', 'f', 'e']})
Monomer('Bax', ['bh3', 'a6', 'loc', 'dye'],
        {'loc': ['c','m', 'i', 'p'],
         'dye': ['none', 'f', 'e']})
Monomer('Vesicles', ['dye'],
        {'dye': ['f', 'e']})
#}}}

#{{{# INITIAL CONDITIONS
#{{{# CONSTANTS
v0 = 0.6/2.0
L_0 = 50000  # Initial condition for lipid concentration
fraction_dim = 0.02
#}}}
# in units of nM
Parameter('tBid_0', 20)
Parameter('Bax_0', 100)
Parameter('Vesicles_0', L_0)
#Parameter('Bax2_0', fraction_dim * Bax_0.value)

Initial(tBid(loc='c', bh3=None, dye='none'), tBid_0)
Initial(Bax(loc='c', bh3=None, a6=None, dye='none'), Bax_0)
Initial(Vesicles(dye='f'), Vesicles_0)
#Initial(Bax(loc='i', bh3=1, dye='f', a6=None) % Bax(loc='i', bh3=1, dye='f', a6=None), Bax2_0)
#Initial(Vesicles(dye='e'), Parameter('eVesicles_0', Vesicles_0.value*0.1))
#}}}

#{{{# OBSERVABLES
Observe('eVes', Vesicles(dye='e'))
Observe('iBax', Bax(loc='i', bh3=None))
Observe('cBax', Bax(loc='c'))
#Observe('pBax', Bax(loc='p'))
#Observe('Bax2', MatchOnce(Bax(loc='i', bh3=1) % Bax(loc='i', bh3=1)))
Observe('mfBax', Bax(dye='f'))
#Observe('meBax', Bax(dye='e'))
#Observe('fVes', Vesicles(dye='f'))
Observe('mftBid', tBid(loc='m', dye='f'))
#Observe('metBid', tBid(loc='m', dye='e'))
Observe('mtBid', tBid(loc='m'))
Observe('mBax', Bax(loc='m'))
#Observe('tBidBax', tBid(bh3=1) % Bax(bh3=1))
#Observe('tBidiBax', tBid(bh3=1) % Bax(bh3=1, loc='i'))

#}}}

#{{{# MODEL MACROS
#{{{# translocate_Bax()
def translocate_Bax():
    # Translocation of proteins to vesicles (i.e., from c to m state)
    # STARTING
    Parameter('Bax_transloc_kf', 1e-4) # nM^-1 s^-1 (Implies 1e5 M^-1 s^-1 forward rate)
    Parameter('Bax_transloc_kr', 7e-2) # Reverse rate (s^-1); implies 700nM binding, per recent Andrews paper

    #Parameter('Bax_transloc_kr', 70) # May be the rate-limiting step for returning Bax to cytosol

    direct_catalysis_reversible(Bax(bh3=None, loc='c', dye='none'), Vesicles(dye='e'),
                                Bax(bh3=None, loc='m', dye='e'), [Bax_transloc_kf, Bax_transloc_kr])
    direct_catalysis_reversible(Bax(bh3=None, loc='c', dye='none'), Vesicles(dye='f'),
                                Bax(bh3=None, loc='m', dye='f'), [Bax_transloc_kf, Bax_transloc_kr])
#}}}

#{{{# translocate_tBid()
def translocate_tBid():
    # Translocation of proteins to vesicles (i.e., from c to m state)
    #Parameter('tBid_transloc_kf', 1e-3) # nM^-1 s^-1 (Implies 1e6 M^-1 s^-1 forward rate)
    #Parameter('tBid_transloc_kr', 7e-1) # Reverse rate (s^-1); implies 700nM binding, per recent Andrews paper
                                        # for Bax
    Parameter('tBid_transloc_kf', 1e-3) # nM^-1 s^-1 (Implies 1e6 M^-1 s^-1 forward rate)
    Parameter('tBid_transloc_kr', 7e-1) # Reverse rate (s^-1); implies 700nM binding, per recent Andrews paper

    ## RULES #####
    # tBid
    direct_catalysis_reversible(tBid(bh3=None, loc='c', dye='none'), Vesicles(dye='f'),
                                tBid(bh3=None, loc='m', dye='f'), [tBid_transloc_kf, tBid_transloc_kr])
    direct_catalysis_reversible(tBid(bh3=None, loc='c', dye='none'), Vesicles(dye='e'),
                                tBid(bh3=None, loc='m', dye='e'), [tBid_transloc_kf, tBid_transloc_kr])
#}}}

#{{{# basal_Bax_activation()
def basal_Bax_activation():
    # Spontaneous rate of transition of Bax from the mitochondrial to the inserted state
    Parameter('basal_Bax_kf', 1e-3) # Implies average time is 10000 seconds???
    Parameter('basal_Bax_kr', 10)

    two_state_equilibrium(Bax(bh3=None, dye='f'), 'm', 'i', [basal_Bax_kf, basal_Bax_kr], sitename='loc')
    two_state_equilibrium(Bax(bh3=None, dye='e'), 'm', 'i', [basal_Bax_kf, basal_Bax_kr], sitename='loc')
#}}}

#{{{# tBid_activates_Bax()
def tBid_activates_Bax(bax_site='bh3'):
    # Andrews suggests that tBid/Bax Kd should work out to 25nM
    Parameter('tBid_mBax_kf', 1e-5) # Forward rate of tBid binding to Bax (E + S -> ES)
    Parameter('tBid_mBax_kr', 2.5e-3)   # Reverse rate of tBid binding to Bax (ES -> E + S)
    Parameter('mBaxtBid_to_iBaxtBid_k', 1e-3) 
    Parameter('tBid_iBax_kr', 2.5e-3)   # Dissociation of tBid from iBax (EP -> E + P)

    # Should add detailed, reversible catalysis macro to make this easier!
    bind(tBid(loc='m', dye='f'), 'bh3', Bax(loc='m', dye='f'), bax_site, # FULL
      [tBid_mBax_kf, tBid_mBax_kr])
    bind(tBid(loc='m', dye='e'), 'bh3', Bax(loc='m', dye='e'), bax_site, # EMPTY
      [tBid_mBax_kf, tBid_mBax_kr])

    # Create the dicts to parameterize the site that tBid binds to
    bax_site_bound = {bax_site:1}
    bax_site_unbound = {bax_site:None}

    # Conformational change of Bax (ES -> EP)
    Rule('mBaxtBid_to_iBaxtBid',
         tBid(loc='m', bh3=1) % Bax(loc='m', **bax_site_bound) >>
         tBid(loc='m', bh3=1) % Bax(loc='i', **bax_site_bound),
         mBaxtBid_to_iBaxtBid_k)

    # tBid dissociates from iBax
    Rule('tBid_unbinds_iBax_f',
         tBid(loc='m', bh3=1) % Bax(loc='i', **bax_site_bound) >>
         tBid(loc='m', bh3=None) + Bax(loc='i', **bax_site_unbound),
         tBid_iBax_kr)
    #Rule('tBid_unbinds_iBax_e', tBid(loc='m', dye='e', bh3=1) % Bax(loc='i', dye='e', bh3=1) >>
    #     tBid(loc='m', dye='e', bh3=None) + Bax(loc='i', dye='e', bh3=None),
    #     tBid_iBax_kr)
#}}}

#{{{# Bax_auto_activates()
def Bax_auto_activates(target_bax_site='a6'):
    # Andrews suggests that tBid/Bax Kd should work out to 25nM
    Parameter('iBax_mBax_kf', 1e-5) # Forward rate of iBax binding to Bax (E + S -> ES)
    Parameter('iBax_mBax_kr', 2.5e-3)   # Reverse rate of iBax binding to Bax (ES -> E + S)
    Parameter('mBaxiBax_to_iBaxiBax_k', 1e-3) 
    #Parameter('iBax_iBax_kr', 2.5e-3)   # Dissociation of iBax from iBax (EP -> E + P)

    # Should add detailed, reversible catalysis macro to make this easier!
    bind(Bax(loc='i', dye='f'), 'bh3', Bax(loc='m', dye='f'), target_bax_site, # FULL
      [iBax_mBax_kf, iBax_mBax_kr])
    bind(Bax(loc='i', dye='e'), 'bh3', Bax(loc='m', dye='e'), target_bax_site, # EMPTY
      [iBax_mBax_kf, iBax_mBax_kr])

    # Create the dicts to parameterize the site that iBax binds to
    target_bax_site_bound = {target_bax_site:1}
    target_bax_site_unbound = {target_bax_site:None}

    # Conformational change of Bax (ES -> EP)
    Rule('mBaxiBax_to_iBaxiBax',
         Bax(loc='i', bh3=1) % Bax(loc='m', **target_bax_site_bound) >>
         Bax(loc='i', bh3=None) + Bax(loc='i', **target_bax_site_unbound),
         mBaxiBax_to_iBaxiBax_k)

    # iBax dissociates from iBax
    #Rule('iBax_unbinds_iBax_f',
    #     Bax(loc='m', bh3=1) % Bax(loc='i', bh3=1) >>
    #     tBid(loc='m', bh3=None) + Bax(loc='i', bh3=None),
    #     tBid_iBax_kr)
    #Rule('tBid_unbinds_iBax_e', tBid(loc='m', dye='e', bh3=1) % Bax(loc='i', dye='e', bh3=1) >>
    #     tBid(loc='m', dye='e', bh3=None) + Bax(loc='i', dye='e', bh3=None),
    #     tBid_iBax_kr)
#}}}

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
    Parameter('tBid_iBax_kf', 1e-3) # Rate of tBid binding to iBax (E + P -> EP)

    # Binding between mtBid and iBax (activation back-reaction--should be slow)
    Rule('tBid_binds_iBax_f', tBid(loc='m', dye='f', bh3=None) + Bax(loc='i', dye='f', bh3=None) >>
         tBid(loc='m', dye='f', bh3=1) % Bax(loc='i', dye='f', bh3=1),
         tBid_iBax_kf)
    Rule('tBid_binds_iBax_e', tBid(loc='m', dye='e', bh3=None) + Bax(loc='i', dye='e', bh3=None) >>
         tBid(loc='m', dye='e', bh3=1) % Bax(loc='i', dye='e', bh3=1),
         tBid_iBax_kf)
#}}}

#{{{# tBid_reverses_Bax()
def tBid_reverses_Bax():
    Parameter('iBaxtBid_to_mBaxtBid_k', 1e-3) # Rate of the EP->ES transition # FIXME

    # REVERSIBILITY OF BAX ACTIVATION BY TBID (EP -> ES)
    Rule('iBaxtBid_to_mBaxtBid',
         tBid(loc='m', bh3=1) % Bax(loc='i', bh3=1) >> tBid(loc='m', bh3=1) % Bax(loc='m', bh3=1),
         iBaxtBid_to_mBaxtBid_k)
#}}}

#{{{# dye_release()
def dye_release(reversible_pore=False):
    # Rate of pore formation/oligomerization of activated Bax (s^-1). 
    # Slowing this down makes Bax linger in the pre-oligomerized (inserted) state.
    # THIS IS CLEARLY WRONG, BECAUSE THIS REACTION SHOULD BE AT LEAST 2ND ORDER
    Parameter('k_pore', 0.1) # was 1


    # Rate of dye release while pores persist. Also a scaling factor for how much pore-state
    # Bax is required to trigger efflux
    # Rate was 100 M^-1 s^-1 in Cecropin A paper, implying 1e-7 in nM^-1a; but that was for smaller liposomes
    # Time taken for dye to be released from a liposome--used to change state of the Bax protein.
    Parameter('k_eflx', 1000)
    # Rate of dye release from liposomes--scaling factor for how much dye is released when a liposome
    # is permeabilized.
    Parameter('lipo_eflx', 100)

    # SECTION RELATED TO DYE RELEASE
    Rule('Free_Bax_Forms_Pores', Bax(loc='i', bh3=None) >> Bax(loc='p', bh3=None), k_pore)

    if (reversible_pore):
        # Rate of closure of pores; making this slower makes pores more stable
        #Parameter('k_pore_rev', 1e-3)
        Parameter('k_pore_rev', 10)
        Rule('Bax_Forms_Pores_rev', Bax(loc='p') >> Bax(loc='i'), k_pore_rev)

    Rule('Bax_Full_to_Empty', Bax(loc='p', dye='f') >> Bax(loc='p', dye='e'), k_eflx)
    #Rule('BaxPore_Releases_Dye', Bax(loc='p') >> Bax(loc='m'), k_pore_rev)
    Rule('Dye_Release', Vesicles(dye='f') + Bax(loc='p', dye='f') >>
          Vesicles(dye='e') + Bax(loc='p', dye='f'), lipo_eflx)
#}}}

#{{{# Bax_dimerizes()
def Bax_dimerizes(dimer_diss_rate=0):
    # Rate of dimerization formation/oligomerization of activated Bax (s^-1). 
    Parameter('Bax_dimerization_kf', 1e-1) # was 1
    Parameter('Bax_dimerization_kr', dimer_diss_rate) 

    Rule('Bax_Forms_Dimers',
         Bax(loc='i', bh3=None, a6=None) + Bax(loc='i', bh3=None, a6=None) <>
         Bax(loc='i', bh3=1, a6=None) % Bax(loc='i', bh3=1, a6=None),
         Bax_dimerization_kf, Bax_dimerization_kr)
#}}}

#{{{# Bax_tetramerizes
def Bax_tetramerizes(tetramer_diss_rate=0):
    """ This function depends on Bax_dimerization to be called as well."""
    # Rate of dimerization formation/oligomerization of activated Bax (s^-1). 
    Parameter('Bax_tetramerization_kf', 1e-1) # was 1
    Parameter('Bax_tetramerization_kr', tetramer_diss_rate) 

    Rule('Bax_Forms_Tetramers',
         Bax(loc='i', bh3=1, a6=None) % Bax(loc='i', bh3=1, a6=None) +
         Bax(loc='i', bh3=2, a6=None) % Bax(loc='i', bh3=2, a6=None) <>
         Bax(loc='i', bh3=1, a6=3) % Bax(loc='i', bh3=1, a6=4) % 
         Bax(loc='i', bh3=2, a6=3) % Bax(loc='i', bh3=2, a6=4), 
         Bax_tetramerization_kf, Bax_tetramerization_kr)
#}}}

#{{{# dye_release_dimeric()
def dye_release_dimeric():
    # Rate of dye release while pores persist. Also a scaling factor for how much pore-state
    # Bax is required to trigger efflux
    # Rate was 100 M^-1 s^-1 in Cecropin A paper, implying 1e-7 in nM^-1a; but that was for smaller liposomes
    # Time taken for dye to be released from a liposome--used to change state of the Bax protein.
    Parameter('k_eflx', 1000)
    # Rate of dye release from liposomes--scaling factor for how much dye is released when a liposome
    # is permeabilized.
    Parameter('lipo_eflx', 100)

    Rule('Bax_Full_to_Empty',
         Bax(loc='i', bh3=1, dye='f') % Bax(loc='i', bh3=1, dye='f') >>
         Bax(loc='i', bh3=1, dye='e') % Bax(loc='i', bh3=1, dye='e'),
         k_eflx)
    Rule('Dye_Release',
         Vesicles(dye='f') + Bax(loc='i', bh3=1, dye='f') % Bax(loc='i', bh3=1, dye='f') >>
         Vesicles(dye='e') + Bax(loc='i', bh3=1, dye='f') % Bax(loc='i', bh3=1, dye='f'),
         lipo_eflx)
#}}}

#{{{# dye_release_tetrameric()
def dye_release_tetrameric():
    # Rate of dye release while pores persist. Also a scaling factor for how much pore-state
    # Bax is required to trigger efflux
    # Rate was 100 M^-1 s^-1 in Cecropin A paper, implying 1e-7 in nM^-1a; but that was for smaller liposomes
    # Time taken for dye to be released from a liposome--used to change state of the Bax protein.
    Parameter('k_eflx', 1000)
    # Rate of dye release from liposomes--scaling factor for how much dye is released when a liposome
    # is permeabilized.
    Parameter('lipo_eflx', 100)

    Rule('Bax_Full_to_Empty_Tetrameric',
         Bax(loc='i', bh3=1, a6=3, dye='f') % Bax(loc='i', bh3=1, a6=4, dye='f') % 
         Bax(loc='i', bh3=2, a6=3, dye='f') % Bax(loc='i', bh3=2, a6=4, dye='f') >>
         Bax(loc='i', bh3=1, a6=3, dye='e') % Bax(loc='i', bh3=1, a6=4, dye='e') % 
         Bax(loc='i', bh3=2, a6=3, dye='e') % Bax(loc='i', bh3=2, a6=4, dye='e'),
         k_eflx)
    Rule('Dye_Release_Tetrameric',
         Vesicles(dye='f') +
         Bax(loc='i', bh3=1, a6=3, dye='f') % Bax(loc='i', bh3=1, a6=4, dye='f') %  # Bax tetramer
         Bax(loc='i', bh3=2, a6=3, dye='f') % Bax(loc='i', bh3=2, a6=4, dye='f') >>
         Vesicles(dye='e') +
         Bax(loc='i', bh3=1, a6=3, dye='f') % Bax(loc='i', bh3=1, a6=4, dye='f') %  # Bax tetramer
         Bax(loc='i', bh3=2, a6=3, dye='f') % Bax(loc='i', bh3=2, a6=4, dye='f'),
         lipo_eflx)
#}}}
#}}}

#{{{# MODEL BUILDING FUNCTIONS
def build_model0():
    translocate_Bax()
    basal_Bax_activation()
    tBid_activates_Bax()
    dye_release()
    #dye_release_dimeric()

def build_model1():
    print "Building model 1: Activation, dimerization"
    translocate_Bax()
    translocate_tBid()
    tBid_activates_Bax(bax_site='a6')
    Bax_dimerizes(dimer_diss_rate=2.5e-1)
    dye_release_dimeric()

#def build_model1r():
#    print "Building model 1r"
#    translocate_Bax()
#    translocate_tBid()
#    tBid_activates_Bax(bax_site='a6')
#    dye_release_dimeric(

def build_model2():
    print "Building model 2: Bax-tBid inhibition"
    translocate_Bax()
    translocate_tBid()
    tBid_activates_Bax()
    Bax_inhibits_tBid()
    Bax_dimerizes(dimer_diss_rate=2.5e-1)
    dye_release_dimeric()

def build_model2r():
    print "Building model 2r"
    translocate_Bax()
    translocate_tBid()
    tBid_activates_Bax()
    Bax_inhibits_tBid()
    dye_release_dimeric(reversible_pore=True)

def build_model3():
    print "Building model 3"
    translocate_Bax()
    translocate_tBid()
    tBid_activates_Bax()
    Bax_inhibits_tBid()
    tBid_reverses_Bax()
    dye_release_dimeric(reversible_pore=False)

def build_model3r():
    print "Building model 3r"
    translocate_Bax()
    translocate_tBid()
    tBid_activates_Bax(bax_site='a6')
    Bax_inhibits_tBid()
    #tBid_reverses_Bax()
    Bax_dimerizes(dimer_diss_rate=0)
    dye_release_dimeric()

def build_model4():
    print "Building model 4"
    translocate_Bax()
    translocate_tBid()
    tBid_activates_Bax(bax_site='a6')
    Bax_inhibits_tBid()
    #tBid_reverses_Bax()
    Bax_auto_activates(target_bax_site='a6') 
    Bax_dimerizes(dimer_diss_rate=0)
    Bax_tetramerizes(tetramer_diss_rate=0)
    #dye_release_dimeric(reversible_pore=False)
    dye_release_tetrameric()
#}}} 

build_model2()

def load_params(param_set):
    for i, param_name in enumerate(param_set.param_dict):
        model_pnames = [param.name for param in model.parameters]
        if (param_name in model_pnames):
            model.parameters[param_name].value = param_set.param_dict[param_name]
        else:
            print(("WARNING: parameter %s in the given parameter set does not " +
                   "exist in this model, it is being ignored.") % param_name)

#{{{# RUNNING THE MODEL
def run_model(tmax=4000, fittype='explin'):
    t = linspace(0, tmax, 100)
    x = odesolve(model, t)
    figure()
    #plot(t, x['metBid']/tBid_0.value, label='metBid')
    plot(t, x['mftBid']/tBid_0.value, label='mftBid')
    plot(t, x['mfBax']/Bax_0.value, label='mfBax')
    #plot(t, x['meBax']/Bax_0.value, label='meBax')
    plot(t, x['tBidBax']/Bax_0.value, label='tBidBax')
    plot(t, x['tBidiBax']/Bax_0.value, label='tBidiBax')
    #plot(t, (x['iBax']+x['pBax'])/Bax_0.value, label='ipBax')
    plot(t, (x['iBax'])/Bax_0.value, label='iBax')
    plot(t, (2*x['Bax2'])/Bax_0.value, label='Bax2')
    #plot(t, x['pBax']/Bax_0.value, label='pBax')
    plot(t, x['eVes']/Vesicles_0.value, label='eVes')
    legend()

    # Fit eVes trajectory
    k = fitting.Parameter(0.0025)
    k2 = fitting.Parameter(0.00025)
    fmax = fitting.Parameter(4)
    fmax2 = fitting.Parameter(0.4)
    m = fitting.Parameter(0.01)

    def linear(t):      return (m()*t)
    def single_exp (t): return ((fmax()*(1 - exp(-k()*t))))
    def exp_lin(t):     return ((fmax()*(1 - exp(-k()*t))) + (m()*t))
    def double_exp(t):  return ((fmax()*(1 - exp(-k()*t)))  + (fmax2()*(1 - exp(-k2()*t))))
    def exp_exp(t):     return ((fmax()*(1 - exp((1- exp(-k()*t))   ))))

    if (True):
        timecourse = []
        for val in x['eVes']/Vesicles_0.value:
            p  = 1 - val
            timecourse.append(-math.log(float(p)))

        time = t
        if (fittype == 'linear'):
            fit(linear, [m], array(timecourse), array(time))
            #fit_initial(single_exp, [k, fmax], array(timecourse), array(time))
            fitfunc = linear
        elif (fittype == 'single_exp'):
            fit(single_exp, [k, fmax], array(timecourse), array(time))
            #fit_initial(single_exp, [k, fmax], array(timecourse), array(time))
            fitfunc = single_exp
        elif (fittype == 'exp_lin'):
            fit(exp_lin, [k, fmax, m], array(timecourse), array(time))
            #fit_initial(exp_lin, [k, fmax, m], array(timecourse), array(time))
            fitfunc = exp_lin
        elif (fittype == 'double_exp'):
            fit(double_exp, [k, fmax, k2, fmax2], array(timecourse), array(time))
            #fit_initial(double_exp, [k, fmax, k2, fmax2], array(timecourse), array(time))
            fitfunc = double_exp
        elif (fittype == 'exp_exp'):
            fit(exp_exp, [k, fmax], array(timecourse), array(time))
            #fit_initial(exp_exp, [k, fmax], array(timecourse), array(time))
            fitfunc = exp_exp
        else:
            raise Exception('unknown fit type')

        mse_val = mse(fitfunc, array(timecourse), array(time))

        # Perform the fit
        #fit(biphasic, [ki_parm, k0_parm, kt_parm], array(timecourse), array(time))
        #fit(biphasic, [ki_parm, kt_parm], array(timecourse), array(time))
        #fit(biphasic, [ki_parm, k0_parm], array(timecourse), array(time))
        #print("k0=" + str(k0_parm()) + ", ki=" + str(ki_parm()) + ", kt=" + str(kt_parm()) )

        # Plot original values along with fitted values
        figure()
        fit_vals = map(fitfunc, time) 
        plot(time, timecourse, label='eVes Model')
        plot(time, fit_vals, label='eVes Fit')
        legend(loc='lower right')
        #return x
#}}}


