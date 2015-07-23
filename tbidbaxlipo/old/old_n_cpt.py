from pysb import *
from pysb.macros import *
import numpy as np
from tbidbaxlipo.util.fitting import fit, mse
from tbidbaxlipo.util import color_iter
from scipy.stats import poisson
from pysb import kappa
from pysb import bng
from tbidbaxlipo.models import core
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties

class Builder(core.Builder):

    def within_compartment_rsf(self):
        return 1.0

    def __init__(self, scaling_factor=10, params_dict=None):
        # Sets self.model = Model(), and self.params_dict
        core.Builder.__init__(self, params_dict=params_dict)

        # This is the only variable that needs to be set to rescale
        # the simulation according to the nominal initial concentrations
        # given below (which match those of the deterministic simulations)
        self.scaling_factor = scaling_factor

        self.parameter('Vesicles_0', 5, factor=self.scaling_factor)
        self.parameter('tBid_0', 20, factor=self.scaling_factor)
        self.parameter('Bax_0', 100, factor=self.scaling_factor)

        self.declare_monomers()

        # COMPARTMENTS
        solution = self.compartment('solution', dimension=3, parent=None)
        for i in range(1, int(self['Vesicles_0'].value)+1):
            self.compartment('c' + str(i), dimension=2, parent=solution)

        tBid = self['tBid']
        Bax = self['Bax']
        self.initial(tBid(bh3=None, loc='c') ** solution, self['tBid_0'])
        self.initial(Bax(loc='c', bh3=None, a6=None, lipo=None,
                     c3='s', c62='s', c120='s', c122='s', c126='s', c184='s')
                     ** solution, self['Bax_0'])

        # OBSERVABLES
        self.observable('ctBid', tBid(loc='c') ** solution)
        self.observable('mtBid', tBid(loc='m'))
        self.observable('cBax', Bax(loc='c') ** solution)
        self.observable('mBax', Bax(loc='m'))

        # Activation
        self.observable('iBax', Bax(loc='i'))
        # FIXME this should be based on the Bax activation site!!!
        self.observable('tBidBax', tBid(loc='m', bh3=1) % Bax(loc='m', bh3=1))

        # pore formation
        for cpt in self.model.compartments:
            if (not cpt.name == 'solution'):
                self.observable('mtBid_%s' % cpt.name, tBid() ** cpt)
                self.observable('Bax_%s' % cpt.name, Bax() ** cpt)
                self.observable('mBax_%s' % cpt.name, Bax(loc='m') ** cpt)
                self.observable('iBax_%s' % cpt.name, Bax(loc='i') ** cpt)
                self.observable('tBidBax_%s' % cpt.name,
                            tBid(loc='m', bh3=1) % Bax(loc='m', a6=1) ** cpt)
                #self.observable('Bax2_%s' % cpt.name,
                #                MatchOnce(Bax(bh3=1) % Bax(bh3=1) ** cpt))

    # Translocation
    def translocate_Bax(self):
        print("n_cpt: translocate_Bax()")

        Bax_transloc_kf = self.parameter('Bax_transloc_kf', 1e-2,
                                          factor=(1/float(self.scaling_factor)))
        Bax_transloc_kr = self.parameter('Bax_transloc_kr', 1e-1)

        Bax = self['Bax']
        solution = self['solution']

        for cpt in self.model.compartments:
            if (not cpt.name == 'solution'):
                self.rule('Bax_translocates_sol_to_%s' % cpt.name,
                     Bax(loc='c', bh3=None, a6=None) ** solution >>
                     Bax(loc='m', bh3=None, a6=None) ** cpt,
                     Bax_transloc_kf)
                self.rule('Bax_translocates_%s_to_sol' % cpt.name,
                     Bax(loc='m', bh3=None, a6=None) ** cpt >>
                     Bax(loc='c', bh3=None, a6=None) ** solution,
                     Bax_transloc_kr)

    def translocate_tBid(self):
        print("n_cpt: translocate_tBid()")

        # Translocation of proteins to vesicles (i.e., from c to m state)
        tBid_transloc_kf = self.parameter('tBid_transloc_kf', 1e-1,
                                          factor=(1/float(self.scaling_factor)))
        tBid_transloc_kr = self.parameter('tBid_transloc_kr', 0)

        tBid = self['tBid']
        solution = self['solution']

        for cpt in self.model.compartments:
            if (not cpt.name == 'solution'):
                self.rule('tBid_translocates_sol_to_%s' % cpt.name,
                     tBid(loc='c', bh3=None) ** solution >>
                     tBid(loc='m', bh3=None) ** cpt,
                     tBid_transloc_kf)
                self.rule('tBid_translocates_%s_to_sol' % cpt.name,
                     tBid(loc='m', bh3=None) ** cpt >>
                     tBid(loc='c', bh3=None) ** solution,
                     tBid_transloc_kr)

    def translocate_tBid_Bax(self):
        print("n_cpt: translocate_tBid_Bax()")
        self.translocate_tBid()
        self.translocate_Bax()

    # Pore formation
    def pores_from_Bax_monomers(self, bax_loc_state='i', reversible=False):
        print("n_cpt: pores_from_Bax_monomers()")

        pore_formation_rate_k = self.parameter('pore_formation_rate_k', 1e-3)

        Bax = self['Bax']
        Pores = self['Pores']

        self.rule('Pores_From_Bax_Monomers', 
             Bax(loc=bax_loc_state) >> Bax(loc='p') + Pores(),
             pore_formation_rate_k)

        if reversible:
            pore_reverse_rate_k = self.parameter('pore_reverse_rate_k', 1e-3)

            self.rule('Pores_reverse',
                 Bax(loc='p') >> Bax(loc='c') ** solution,
                 pore_reverse_rate_k)

        # Pore observables
        self.observable('pBax', Bax(loc='p'))
        self.observable('pores', Pores())
        for i, cpt in enumerate(self.model.compartments):
            if (not cpt.name == 'solution'):
                self.observable('pores_%s' % cpt.name, Pores() ** cpt)

    def pores_from_Bax_dimers(self):
        """Basically a way of counting the time-integrated amount of
           forward pore formation."""
        raise NotImplementedError()
        Parameter('pore_formation_rate_k', 100)
        #Parameter('scaled_pore_formation_rate_k', 0.03)
        Parameter('pore_recycling_rate_k', 100)

        Rule('Pores_From_Bax_Dimers', 
             MatchOnce(Bax(loc='i', bh3=1, a6=None) %
                       Bax(loc='i', bh3=1, a6=None)) >>
             MatchOnce(Bax(loc='p', bh3=1, a6=None) %
                       Bax(loc='p', bh3=1, a6=None)) + Pores(),
             pore_formation_rate_k)

    # Running the model
    def get_compartment_observables(self, obs_basename):
        """Get the names for an observable across all compartments.

        Parameters
        ----------
        obs_basename : string
            The basename of the desired observable, e.g. 'mBax'.

        Returns
        -------
        list of strings
            A list of strings of the observable names for the desired
            observable in each compartment, e.g. ['mBax_c1', 'mBax_c2', ...]
        """

        return ['%s_%s' % (obs_basename, cpt.name)
                for cpt in self.model.compartments
                if not cpt.name == 'solution']

    def run_model(self, tmax=12000, num_sims=1):
        xrecs = []
        dr_all = []
        for i in range(0, num_sims):
            print "Running SSA simulation %d of %d..." % (i+1, num_sims)
            ssa_result = bng.run_ssa(self.model, t_end=tmax, n_steps=100,
                                     cleanup=True)
            xrecs.append(ssa_result)
            dr_all.append(get_dye_release(self.model, 'pores', ssa_result))

        xall = np.array([x.tolist() for x in xrecs])
        x_std = np.recarray(xrecs[0].shape, dtype=xrecs[0].dtype,
                            buf=np.std(xall, 0))
        x_avg = np.recarray(xrecs[0].shape, dtype=xrecs[0].dtype,
                            buf=np.mean(xall, 0))

        plot_simulation(x_avg, x_std, dr_all, self.model, num_sims)

def plot_simulation(x_avg, x_std, dr_all, model, num_sims, figure_ids=[0, 1]):
    ci = color_iter()
    marker = ','
    linestyle = ''

    tBid_0 = model.parameters['tBid_0'].value
    Bax_0 = model.parameters['Bax_0'].value
    Vesicles_0 = model.parameters['Vesicles_0'].value

    # Translocation
    plt.figure(figure_ids[0])
    plt.errorbar(x_avg['time'], x_avg['ctBid']/tBid_0,
             yerr=(x_std['ctBid']/tBid_0)/np.sqrt(num_sims),
             label='ctBid',
             color=ci.next(), marker=marker, linestyle=linestyle)
    plt.errorbar(x_avg['time'], x_avg['mtBid']/tBid_0,
             yerr=(x_std['mtBid']/tBid_0)/np.sqrt(num_sims),
             label='mtBid',
             color=ci.next(), marker=marker, linestyle=linestyle)
    plt.errorbar(x_avg['time'], x_avg['cBax']/Bax_0,
             yerr=(x_std['cBax']/Bax_0)/np.sqrt(num_sims),
             label='cBax',
             color=ci.next(), marker=marker, linestyle=linestyle)
    plt.errorbar(x_avg['time'], x_avg['mBax']/Bax_0,
             yerr=(x_std['mBax']/Bax_0)/np.sqrt(num_sims),
             label='mBax',
             color=ci.next(), marker=marker, linestyle=linestyle)

    # Activation
    plt.errorbar(x_avg['time'], x_avg['iBax']/Bax_0,
                 yerr=(x_std['iBax']/Bax_0)/np.sqrt(num_sims),
                 label='iBax', color=ci.next(), marker=marker,
                 linestyle=linestyle)
    plt.errorbar(x_avg['time'], x_avg['tBidBax']/tBid_0,
                 yerr=(x_std['tBidBax']/tBid_0)/np.sqrt(num_sims),
                 label='tBidBax', color=ci.next(), marker=marker,
                 linestyle=linestyle)

    # Pore Formation
    plt.errorbar(x_avg['time'], x_avg['pBax']/Bax_0,
                 yerr=(x_std['pBax']/Bax_0)/np.sqrt(num_sims),
                 color=ci.next(),
                 label='pBax', marker=marker, linestyle=linestyle)

    ci = color_iter()

    # Dye release calculated exactly ----------
    dr_avg = np.mean(dr_all, 0)
    dr_std = np.std(dr_all, 0)
    plt.errorbar(x_avg['time'], dr_avg,
                yerr=dr_std/np.sqrt(num_sims), label='dye_release',
                color=ci.next(), marker=marker, linestyle=linestyle)

    fontP = FontProperties()
    fontP.set_size = ('small')
    plt.legend(loc='upper center', prop=fontP, ncol=5,
               bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=True)
    plt.xlabel("Time (seconds)")
    plt.ylabel("Normalized Concentration")
    plt.show()

    # Plot pores/vesicle in a new figure ------
    ci = color_iter()
    plt.figure(figure_ids[1])
    plt.errorbar(x_avg['time'], x_avg['pores']/Vesicles_0,
                 yerr=(x_std['pores']/Vesicles_0)/np.sqrt(num_sims),
                 label='pores/ves', color=ci.next(), marker=marker,
                 linestyle=linestyle)
    plt.legend(loc='upper center', prop=fontP, ncol=1,
               bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=True)
    plt.show()
    #F_t = 1 - dr_avg
    #pores_poisson = -log(F_t)
    #plot(x_avg['time'], pores_poisson, color=ci.next(), label='-ln F(t),
    #     stoch', marker=marker, linestyle=linestyle)
    #xlabel("Time (seconds)")
    #ylabel("Pores/vesicle")
    #title("Pores/vesicle")
    #legend()

    #xlabel("Time (seconds)")
    #ylabel("Dye Release")
    #title("Dye release calculated via compartmental model")

    return

# COMMENTS
"""
Things I would want to know:
* Is pore formation actually independent?
  - To know this would want to know if pores per liposome actually follows a Poisson distribution
    under various assumptions
* What is the distribution of Bax molecules among liposomes?
  - At steady state, 
* 
SOME OBSERVATIONS
* When looking simply at binding and unbinding, one finds that, especially with large
numbers of compartments, the distribution of molecules localized to liposomes follows
a Poisson distribution. At stationarity, one cal average out the distribution over
time and get a smooth distribution that way.

* Notably (perhaps) when an irreversible step is in place, e.g., the irreversible
insertion of Bax, the uneven distribution among compartments can get "frozen" in
place. So there is no smoothing effect by averaging over many timepoints, since
the distribution of inserted Bax is fixed once all of the Bax has been re-localized.

* So, the next step is to consider pore formation. The first thing to do would be
to have a forward rate of pore formation that was derived from Bax activation.
In the trivial case of pores simple being synonymous with iBax, one could see 
average "pores per liposome" by plotting the compartment mean, and could see the
distribution at any given timepoint and evaluate it's Poisson-ness.

* In the case of reversible Bax, one would expect to see the pore distribution
become more clearly/smoothly Poisson-like as time went on, because the pore
number would increase and pore "events" would be more smoothly distributed
among liposomes.

* What about when Bax can auto-activate? There are at least two cases to
 consider:

1) when Bax can "grow" pores, depleting Bax from the solution while not triggering
new (initial) pore forming events

2) when Bax can auto-recruit, potentially triggering additional pore
formation events (either initial or secondary)

The question is, what is the distribution of Bax localization, and
what is the distribution of initial pore formation? Is the estimate
gained from f(0) reliable in each of these cases?

TODO:
* Need to update run_ssa code to output BNG output as it runs
* Need to commit and push changes to bng.py
* Need to figure out how to convert rates in det vs. stochastic model

Converting rates:
- In a deterministic model, there is a rate defining the propensity of
the transition c -> m. This is modeled as a pseudo-first order reaction
that is proportional to the amount of lipid/liposomes in the system.
In the stochastic model, there is a propensity of a molecule to jump
to a particular compartment; the summation of the propensity over all
of the compartments should be equal to the propensity of the transition
in the two compartment model, with a stochastic conversion. This is a property
of the law of probability that P(A or B or C) = P(A) + P(B) + P(C). 

So then we need a way of converting from the deterministic rate to the
stochastic rate.

Deterministic rate constant k (in (nmol/L)-1 s-1)

k/V = gamma (nanomoles-1 s-1)

The volume in the titration expt is 60uL, giving

gamma_mol = k / 60 10^-6  L = k / 60*10^-6 (nmol * s)

(This will be larger than K)

gamma_molec = gamma_mol / A = k / (60*10^-6 * 10^-9 moles * sec * 6.022 * 10^23)
gamma_molec = k / (60^10^-6 * 6.022 * 10^14 * molec * sec)
gamma_molec = k / 60 * 6.022 * 10^8 * molec * sec)
gamma_molec = k / 361.32 * 10^8 * molec * sec)
gamma_molec = k / 3.6132 * 10^10 * molec * sec)

"""

# OTHER
def rate_calcs(nm_rate):

    vol = 60e-6
    print "Assumed volume: %g L" % vol

    nm_conc = 100 # assume nM
    molar_conc = nm_conc * 1e-9
    print "Example concentration: %g Molar" % molar_conc

    A = 6.022e23 # Avogadro's number
    total_n = molar_conc * vol * A
    print "Total molecules at molar conc of %g in vol of %g: %g" % \
          (molar_conc, vol, total_n)

    # To convert nanomolar rate to molar rate, multiply the
    # nanomolar rate by 1e9
    print "deterministic (nm-1 s-1) rate: %g" % nm_rate # assume nm^-1 s^-1
    molar_rate = nm_rate * 1e9
    print "deterministic (M-1 s-1) rate: %g" % molar_rate
    molec_rate = molar_rate / (vol * A)
    print "molecular (molec^-1 s^-1) rate: %g" % molec_rate

    scaled_n = 1000
    print("scaled number of molecules meant to represent conc of %g Molar: %g "\
          % (molar_conc, scaled_n))
    actual_to_scaled = scaled_n / float(total_n)
    print "stochastic scaling factor of %g" % actual_to_scaled
    scaled_molec_rate = molec_rate / actual_to_scaled
    print "scaled molec rate: %g" % scaled_molec_rate

def run_model_ode(tmax=12000):
    from pysb.integrate import odesolve
    t = np.linspace(0, tmax, 1000)
    x = odesolve(model, t)

    ci = color_iter()
    plt.figure(1)
    # Translocation
    plt.plot(t, (x['ctBid'])/tBid_0.value, label='ctBid', color=ci.next())
    plt.plot(t, (x['mtBid'])/tBid_0.value, label='mtBid', color=ci.next())
    plt.plot(t, (x['cBax'])/Bax_0.value, label='cBax', color=ci.next())
    plt.plot(t, (x['mBax'])/Bax_0.value, label='mBax', color=ci.next())
    plt.legend()
