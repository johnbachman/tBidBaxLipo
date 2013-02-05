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
"""Note on vesicle concentration SATSOURA:
It was found that the measured r was often higher than the
estimated r, with large variations observed from preparation-to-preparation.
This discrepancy was mainly the result of a lipid concentration that was
lower than expected, showing that during the resuspension and extrusion steps,
not unexpectedly a significant amount of lipid was not incorporated into
liposomes. This is one reason why different repeats of the same experiment
might give different values for bound Bax, since this value depends on r, and
the actual r can vary from experiment-to-experiment."""
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

    """SATSOURA: Immunoblotting suggests that in the absence of tBid, less than
    ~10% of Bax or EGFP-Bax is stably bound to liposomes (data not shown and [50]).
    However, this does not rule out the possibility that some Bax may transiently
    bind to membranes in the absence of tBid, an interac- tion that could be
    disrupted during gel filtration. Indeed, the presence of lipid membranes alone
    can induce a transient Bax conformational change [17], and in healthy cells a
    small fraction of Bax is found loosely associated with mitochondrial membranes
    (as shown by the fact that this fraction is carbonate extractable) [10,32]

    SATSOURA: However, when observing the diffusion of the EGFP-Bax in the presence of lipo-
    somes, but in the absence of tBid (Fig. 2F), it appeared that the total amount
    of EGFP-Bax associated with the liposomes is always less than 5% (for r>1), or
    that the interaction is too transient (lasting less than a few ms) to be
    detected by FCS.
    """
    """ SATSOURA: In contrast, im- munoblotting against tBid showed that for
        all proteins-to-liposome ra- tios explored, about 90% of tBid was associated
        with the membrane (data not shown).  """
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

    # Should also set tBid state?
    # How do I figure out how much/how fast to set tBid from f to e? 
    # Using same rate as for lipos, the amount of tBid that transitions from f to e
    # is proportional to the amount of lipo that transitions.
    # So the idea is that tBid ** full_ves is spread equally across these liposomes
    # and the decay rate of tBid f->e is the same as the decay rate of vesicles
    # from f->e.
    #Rule('tBid_Full_to_Empty',
    #     tBid() ** full_ves  + pore_forming_species ** full_ves >>
    #     tBid() ** empty_ves + pore_forming_species ** full_ves,
    #     lipo_eflx, move_connected=True)

    """

Interestingly, when pore formation is irreversible, this model of dye release
(with dimers) produces the same steady state fraction of dye release,
regardless of the amount of vesicles (given the same rates for all of the Bax
and dye release steps). This must be because the decay of the vesicles from
full to empty essentially tracks the Bax dimerization rate???

Part of the problem here though is that the rates of tBid/Bax interaction are not
controlled for the amount of lipid there is. With fewer liposomes, there is a greater
tBid/Bax to liposome ratio, and hence encounter frequencies will be higher.
So forward rates should be scaled by this concentration; they will then be
in terms of molar ratios of tBid/Bax to liposome/lipid.
   
In addition, lipid amounts will affect the translocation rates of tBid and Bax
to membranes. This is currently taken into account though in the translocation
steps.
 
"""
               
 


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


#  COMMENTS
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

# EXTRA
"""
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

    if (False):
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
    return x
"""
