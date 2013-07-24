__author__ = 'johnbachman'

# IMPORTS
from pysb import *
from pysb.macros import *
from pylab import *
from pysb.integrate import odesolve
from util.fitting import fit, fit_initial, mse
import util.fitting as fitting

model = Model('tBidBax_2c_model')

from model_core import *

# COMPARTMENTS
Compartment('solution', dimension=3, parent=None)
Compartment('full_ves', dimension=2, parent=solution)
Compartment('empty_ves', dimension=2, parent=solution)

# INITIAL CONDITIONS
# Concentration of vesicles, in nanomolar
#Parameter('NUM_VESICLES', 0.038) # 1000nm vesicles
#Parameter('NUM_VESICLES', 5) # 100nm vesicles
Parameter('Vesicles_0', 5)
Parameter('tBid_0', 20)
Parameter('Bax_0', 100)

Initial(tBid(loc='c', bh3=None) ** solution, tBid_0)
Initial(Bax(loc='c', bh3=None, a6=None) ** solution, Bax_0)
Initial(Vesicles() ** full_ves, Vesicles_0)
#Initial(Bax(loc='i', bh3=1, dye='f', a6=None) % Bax(loc='i', bh3=1, dye='f', a6=None), Bax2_0)
#Initial(Vesicles(dye='e'), Parameter('eVesicles_0', Vesicles_0.value*0.1))

# OBSERVABLES
#Observable('eVes', Vesicles(dye='e'))
#Observable('pores', Pores())
Observable('ctBid', tBid(loc='c'))
Observable('mtBid', tBid(loc='m'))
Observable('cBax', Bax(loc='c'))
Observable('mBax', Bax(loc='m'))
#Observable('iBax', Bax(loc='i'))
#Observable('eiBax', Bax(loc='i') ** empty_ves)
#Observable('ePore', Pore() ** empty_ves)
#Observable('fPore', Pore() ** full_ves)
#Observable('pBax', Bax(loc='p'))
#Observable('pBax', Bax(loc='p'))
#Observable('Bax2', MatchOnce(Bax(bh3=1, a6=None) % Bax(bh3=1, a6=None)))
#Observable('Bax4', MatchOnce(
#    Bax(bh3=1, a6=3) % Bax(bh3=1, a6=4) %
#    Bax(bh3=2, a6=3) % Bax(bh3=2, a6=4)))
#Observable('mfBax', Bax(dye='f'))
#Observable('meBax', Bax(dye='e'))
#Observable('eVes', Vesicles() ** empty_ves)
#Observable('fVes', Vesicles() ** full_ves)
#Observable('mftBid', tBid(loc='m', dye='f'))
#Observable('metBid', tBid(loc='m', dye='e'))
#Observable('tBidBax', tBid(bh3=1) % Bax(a6=1))
#Observable('tBidiBax', tBid(bh3=1) % Bax(bh3=1, loc='i'))

# MODEL MACROS
def translocate_tBid_Bax_2c():
    print("translocate_tBid_Bax_2c()")

 
    # DEPENDENCIES: Expects the compartments full_ves and empty_ves to be defined

    # Translocation of proteins to vesicles (i.e., from c to m state)
    # Lovell measures translocation rate of tBid at 10^12 molecules per sec
    # = 10^12/10^23 = 10^-11 = 0.01nM/sec
    # k * 20nM tBid * ves = 0.01nm/sec
    # k = 0.01 / (20 * ves)
    # for 5nm vesicles: k = 0.01 / (100) = 1e-4
    # in nM:
    # diffusion limited: 10^6 M^-1 s^-1 = 10^-3
    # reverse rate:
    Parameter('tBid_transloc_kf', 1e-2) # nM^-1 s^-1 (Implies 1e6 M^-1 s^-1 forward rate)
    Parameter('tBid_transloc_kr', 0) # Reverse rate (s^-1)

    Parameter('Bax_transloc_kf', 1e-3) # nM^-1 s^-1 (Implies 1e5 M^-1 s^-1 forward rate)
    Parameter('Bax_transloc_kr', 1e-2) # Reverse rate (s^-1); implies 700nM binding,
                                       # per recent Andrews paper

    for cpt in [full_ves, empty_ves]:
        Rule('tBid_translocates_sol_to_%s' % cpt.name,
             Vesicles() ** cpt + tBid(loc='c') ** solution >>
             Vesicles() ** cpt + tBid(loc='m') ** cpt,
             tBid_transloc_kf)
        Rule('tBid_translocates_%s_to_sol' % cpt.name,
             tBid(loc='m', bh3=None) ** cpt >> tBid(loc='c', bh3=None) ** solution,
             tBid_transloc_kr)
        Rule('Bax_translocates_sol_to_%s' % cpt.name,
             Vesicles() ** cpt + Bax(loc='c') ** solution >>
             Vesicles() ** cpt + Bax(loc='m') ** cpt,
             Bax_transloc_kf)
        Rule('Bax_translocates_%s_to_sol' % cpt.name,
             Bax(loc='m', bh3=None, a6=None) ** cpt >>
             Bax(loc='c', bh3=None, a6=None) ** solution,
             Bax_transloc_kr)

def dye_release(pore_forming_species):
    print("dye_release(" + str(pore_forming_species) + ")")
    # Rate of dye release while pores persist. Also a scaling factor for how much pore-state
    # Bax is required to trigger efflux
    # Slowing this down makes Bax linger in the pre-oligomerized (inserted) state.
    # Rate was 100 M^-1 s^-1 in Cecropin A paper, implying 1e-7 in nM^-1; but that was for
    # smaller liposomes
    # Time taken for dye to be released from a liposome--used to change state of the Bax protein.
    Parameter('k_eflx', 1000)

    # Rate of dye release from liposomes--scaling factor for how much dye is
    # released when a liposome is permeabilized.
    #Parameter('lipo_eflx', k_eflx.value / Vesicles_0.value)
    Parameter('lipo_eflx', 10)

    # Note
    Rule('Vesicles_Full_to_Empty',
         Vesicles() ** full_ves  + pore_forming_species ** full_ves >>
         Vesicles() ** empty_ves + pore_forming_species ** full_ves,
         lipo_eflx)

    Rule('Bax_Full_to_Empty',
         pore_forming_species ** full_ves >> pore_forming_species ** empty_ves,
         k_eflx)


######################################


# MODEL BUILDING FUNCTIONS
def build_2c_model0():
    print "--------------------------------"
    print "Building 2c model 0:"
    translocate_tBid_Bax_2c()
    #tBid_activates_Bax(bax_site='a6', vesicles_conc=Vesicles_0)
    #tBid_activates_Bax(bax_site='a6')
    #Bax_dimerizes()
    #dye_release(Bax(loc='i', bh3=None))
    #dye_release(Pore())
    #dye_release(Bax(loc='i', bh3=1, a6=None) % Bax(loc='i', bh3=1, a6=None))
 
def build_model1():
    print """Building model 1: Translocation, Bax activation, and monomeric pore formation,
             and aggregation."""
    translocate_Bax()
    translocate_tBid()
    tBid_activates_Bax(bax_site='a6')
    pores_from_Bax_monomers()
    #Bax_aggregates_at_pores()

def build_model1r():
    translocate_Bax()
    translocate_tBid()
    tBid_activates_Bax(bax_site='a6')
    pores_from_Bax_monomers()
    Bax_reverses()
    #Bax_aggregates_at_pores()

def build_model2():
    print "Building model 2: Translocation, Bax activation, and dimeric pore formation."
    translocate_Bax()
    translocate_tBid()
    tBid_activates_Bax(bax_site='a6')
    Bax_dimerizes()
    pores_from_Bax_dimers()
    Bax_aggregates_at_pores()

def build_model3():
    print "Building model 3: Translocation, Bax activation, "
    print "tetrameric pore formation and aggregation."
    translocate_Bax()
    translocate_tBid()
    tBid_activates_Bax(bax_site='a6')
    Bax_dimerizes()
    Bax_tetramerizes()
    pores_from_Bax_tetramers()
    Bax_reverses()
    #Bax_auto_activates(target_bax_site='a6')
    #Bax_aggregates_at_pores()

def build_model4():
    print "Building model 4: Translocation, Bax activation/auto-activation, "
    print "tetrameric pore formation."

    """ In this model, the autoactivation has no effect, perhaps because
        dimerization is so favorable that there is not sufficient free iBax/pBax to
        kick in the positive feedback effect."""

    translocate_Bax()
    translocate_tBid()
    tBid_activates_Bax(bax_site='a6')
    Bax_reverses()
    Bax_dimerizes()
    Bax_tetramerizes()
    pores_from_Bax_tetramers()
    #Bax_auto_activates(target_bax_site='a6')
    #Bax_inhibits_tBid()

build_2c_model0()

def load_params(param_set):
    for i, param_name in enumerate(param_set.param_dict):
        model_pnames = [param.name for param in model.parameters]
        if (param_name in model_pnames):
            model.parameters[param_name].value = param_set.param_dict[param_name]
        else:
            print(("WARNING: parameter %s in the given parameter set does not " +
                   "exist in this model, it is being ignored.") % param_name)

# RUNNING THE MODEL
def run_model(tmax=12000, fittype='explin'):
    t = linspace(0, tmax, 1000)
    x = odesolve(model, t)
    figure()
    plot(t, (x['ctBid'])/tBid_0.value, label='ctBid')
    plot(t, (x['mtBid'])/tBid_0.value, label='mtBid')
    plot(t, (x['cBax'])/Bax_0.value, label='cBax')
    plot(t, (x['mBax'])/Bax_0.value, label='mBax')
    #plot(t, x['metBid']/tBid_0.value, label='metBid')
    #plot(t, x['mtBid']/tBid_0.value, label='mtBid')
    #plot(t, x['mftBid']/tBid_0.value, label='mftBid')
    #plot(t, x['mfBax']/Bax_0.value, label='mfBax')
    #plot(t, x['meBax']/Bax_0.value, label='meBax')
    #plot(t, x['tBidBax']/Bax_0.value, label='tBidBax')
    #plot(t, x['tBidiBax']/Bax_0.value, label='tBidiBax')
    #plot(t, (x['iBax']+x['pBax'])/Bax_0.value, label='ipBax')
    #plot(t, (x['iBax'])/Bax_0.value, label='iBax')
    #plot(t, (x['eiBax'])/Bax_0.value, label='eiBax')
    #plot(t, (2*x['ePore'])/Bax_0.value, label='ePore')
    #plot(t, (2*x['fPore'])/Bax_0.value, label='fPore')
    #plot(t, (x['pBax'])/Bax_0.value, label='pBax')
    #plot(t, (2*(x['Bax2']+(2*x['Bax4'])))/Bax_0.value, label='Bax2')
    #plot(t, (2*x['Bax2'])/Bax_0.value, label='Bax2')
    #plot(t, (4*x['Bax4'])/Bax_0.value, label='Bax4')
    #plot(t, x['pBax']/Bax_0.value, label='pBax')
    #plot(t, x['eVes']/Vesicles_0.value, label='eVes')
    legend()
    xlabel("Time (seconds)")
    ylabel("Normalized Concentration")

    #figure()
    #plot(t, (x['pores']/model.parameters['NUM_VESICLES'].value), label='pores')
    #xlabel("Time (seconds)")
    #ylabel("Pores per vesicle")


