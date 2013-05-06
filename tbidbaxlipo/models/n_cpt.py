__author__ = 'johnbachman'

from pysb import *
from pysb.macros import *
import numpy as np
#from pylab import *
from tbidbaxlipo.util.fitting import fit, mse
from tbidbaxlipo.util import color_iter 
from scipy.stats import poisson
from pysb import kappa
from pysb import bng
from tbidbaxlipo.models import core
from matplotlib import pyplot as plt

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
        for i in range(1, self['Vesicles_0'].value+1):
            self.compartment('c' + str(i), dimension=2, parent=solution)

        tBid = self['tBid']
        Bax = self['Bax']

        self.initial(tBid(bh3=None, loc='c') ** solution, self['tBid_0'])
        self.initial(Bax(loc='c', bh3=None, a6=None,
                     c3='s', c62='s', c120='s', c122='s', c126='s', c184='s')
                     ** solution, self['Bax_0'])

        # OBSERVABLES
        self.observable('ctBid', tBid(loc='c') ** solution)
        self.observable('mtBid', tBid(loc='m'))
        self.observable('cBax', Bax(loc='c') ** solution)
        self.observable('mBax', Bax(loc='m'))

        # Activation
        self.observable('iBax', Bax(loc='i'))
        self.observable('tBidBax', tBid(loc='m', bh3=1) % Bax(loc='m', a6=1))

        # pore formation
        for cpt in self.model.compartments:
            if (not cpt.name == 'solution'):
                self.observable('tBid_%s' % cpt.name, tBid() ** cpt)
                self.observable('Bax_%s' % cpt.name, Bax() ** cpt)
                self.observable('iBax_%s' % cpt.name, Bax(loc='i') ** cpt)
                self.observable('tBidBax_%s' % cpt.name,
                            tBid(loc='m', bh3=1) % Bax(loc='m', a6=1) ** cpt)
                #self.observable('Bax2_%s' % cpt.name,
                #                MatchOnce(Bax(bh3=1) % Bax(bh3=1) ** cpt))

    # MODEL MACROS
    def translocate_tBid_Bax(self):
        print("n_cpt: translocate_tBid_Bax()")

        # Translocation of proteins to vesicles (i.e., from c to m state)
        tBid_transloc_kf = self.parameter('tBid_transloc_kf', 1e-1,
                                          factor=(1/float(self.scaling_factor)))
        tBid_transloc_kr = self.parameter('tBid_transloc_kr', 0)

        Bax_transloc_kf = self.parameter('Bax_transloc_kf', 1e-2,
                                          factor=(1/float(self.scaling_factor)))
        Bax_transloc_kr = self.parameter('Bax_transloc_kr', 1e-2)

        tBid = self['tBid']
        Bax = self['Bax']
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
                self.rule('Bax_translocates_sol_to_%s' % cpt.name,
                     Bax(loc='c', bh3=None, a6=None) ** solution >>
                     Bax(loc='m', bh3=None, a6=None) ** cpt,
                     Bax_transloc_kf)
                self.rule('Bax_translocates_%s_to_sol' % cpt.name,
                     Bax(loc='m', bh3=None, a6=None) ** cpt >>
                     Bax(loc='c', bh3=None, a6=None) ** solution,
                     Bax_transloc_kr)

    def pores_from_Bax_monomers(self):
        print("pores_from_Bax_monomers()")

        Parameter('pore_formation_rate_k', 1e-3)
        #Parameter('pore_recycling_rate_k', 0)

        Rule('Pores_From_Bax_Monomers', 
             Bax(loc='i') >> Bax(loc='p') + Pores(),
             pore_formation_rate_k)

        # Pore observables
        Observable('pBax', Bax(loc='p'))
        Observable('pores', Pores())    
        for i, cpt in enumerate(model.compartments):
            if (not cpt.name == 'solution'):
                Observable('pores_%s' % cpt.name, Pores() ** cpt)

    def pores_from_Bax_dimers(self):
        """Basically a way of counting the time-integrated amount of
           forward pore formation."""

        Parameter('pore_formation_rate_k', 100)
        #Parameter('scaled_pore_formation_rate_k', 0.03)
        Parameter('pore_recycling_rate_k', 100)

        Rule('Pores_From_Bax_Dimers', 
             MatchOnce(Bax(loc='i', bh3=1, a6=None) %
                       Bax(loc='i', bh3=1, a6=None)) >>
             MatchOnce(Bax(loc='p', bh3=1, a6=None) %
                       Bax(loc='p', bh3=1, a6=None)) + Pores(),
             pore_formation_rate_k)

    def run_model(self, tmax=12000, num_sims=1, figure_ids=[0, 1]):
        xrecs = []
        dr_all = []
        for i in range(0, num_sims):
            ssa_result = bng.run_ssa(self.model, t_end=tmax, n_steps=100,
                                     cleanup=True)
            xrecs.append(ssa_result)
            #dr_all.append(get_dye_release(model, 'pores', ssa_result))

        xall = np.array([x.tolist() for x in xrecs])
        x_std = np.recarray(xrecs[0].shape, dtype=xrecs[0].dtype,
                            buf=np.std(xall, 0))
        x_avg = np.recarray(xrecs[0].shape, dtype=xrecs[0].dtype,
                            buf=np.mean(xall, 0))

        ci = color_iter()
        marker = ','
        linestyle = ''
        plt.figure(1)

        tBid_0 = self['tBid_0']
        Bax_0 = self['Bax_0']

        # Translocation
        plt.figure(figure_ids[0])
        plt.errorbar(x_avg['time'], x_avg['ctBid']/tBid_0.value,
                 yerr=x_std['ctBid']/tBid_0.value,
                 color=ci.next(), marker=marker, linestyle=linestyle)
        plt.errorbar(x_avg['time'], x_avg['mtBid']/tBid_0.value,
                 yerr=x_std['mtBid']/tBid_0.value,
                 color=ci.next(), marker=marker, linestyle=linestyle)
        plt.errorbar(x_avg['time'], x_avg['cBax']/Bax_0.value,
                 yerr=x_std['cBax']/Bax_0.value,
                 color=ci.next(), marker=marker, linestyle=linestyle)
        plt.errorbar(x_avg['time'], x_avg['mBax']/Bax_0.value,
                 yerr=x_std['mBax']/Bax_0.value,
                 color=ci.next(), marker=marker, linestyle=linestyle)

        # Activation
        plt.errorbar(x_avg['time'], x_avg['iBax']/Bax_0.value,
                 yerr=x_std['iBax']/Bax_0.value, label='iBax',
                 color=ci.next(), marker=marker, linestyle=linestyle)
        plt.errorbar(x_avg['time'], x_avg['tBidBax']/tBid_0.value,
                 yerr=x_std['tBidBax']/tBid_0.value,
                 color=ci.next(), marker=marker, linestyle=linestyle)

        # Dye release calculated exactly ----------
        #dr_avg = mean(dr_all, 0)
        #dr_std = std(dr_all, 0)
        #errorbar(x_avg['time'], dr_avg,
        #         yerr=dr_std, label='dye_release',
        #         color=ci.next(), marker=marker, linestyle=linestyle)


        # Pore Formation
        #plot(x['time'], x['pBax']/Bax_0.value, label='pBax')
        #leg = legend()
        #ltext = leg.get_texts()
        #setp(ltext, fontsize='small')

        #xlabel("Time (seconds)")
        #ylabel("Normalized Concentration")

        #ci = color_iter()
        # Plot pores/vesicle in a new figure ------
        #figure(2)
        #errorbar(x_avg['time'], x_avg['pores'] / float(NUM_COMPARTMENTS),
        #         yerr=x_std['pores']/float(NUM_COMPARTMENTS), label='pores',
        #         color=ci.next(), marker=marker, linestyle=linestyle)

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

        return xrecs[0]

# SSA data analysis functions
def plot_compartment_mean(model, observable_name, output):
    total = 0
    for i, cpt in enumerate(model.compartments):
        if (not cpt.name == 'solution'):
            total += output['%s_%s' % (observable_name, cpt.name)]
    mean = total / len(model.compartments)
    plt.figure()
    plt.plot(output['time'], mean)

def get_dye_release(model, pore_observable_name, output):
    num_timepoints = len(output['time'])
    dye_release = zeros(num_timepoints)
    for t in range(0, num_timepoints):
        permeabilized_count = 0.
        for i, cpt in enumerate(model.compartments):
            if (not cpt.name == 'solution' and 
                output['%s_%s' % (pore_observable_name, cpt.name)][t] > 0):
                permeabilized_count += 1
        frac_permeabilized = permeabilized_count / float(NUM_COMPARTMENTS)
        dye_release[t] = frac_permeabilized

    return dye_release

def plot_compartment_all(model, observable_name, output):
    plt.figure()
    for i, cpt in enumerate(model.compartments):
        if (not cpt.name == 'solution'):
            plt.plot(output['time'], output['%s_%s' %
                                            (observable_name, cpt.name)])

def plot_predicted_p(model, observable_name, output):
    # Iterate through data, and at each time point, get the fraction
    # of compartments with 0 of the observable
    num_timepoints = len(output['time']) 
    f0 = []
    #for t in range(0, num_timepoints):
        #f0_t
        #for i, cptgg

def get_compartment_dist(model, observable_name, output):

    end_counts = []
    num_timepoints = len(output['time']) 
    for t in range(num_timepoints-1, num_timepoints):
        for i, cpt in enumerate(model.compartments):
            if (not cpt.name == 'solution'):
                end_counts.append(output['%s_%s' %
                                         (observable_name, cpt.name)][t])
    return np.array(end_counts)

def combine_dists(model, observable_name, output_array):
    totals = []
    for output in output_array:
        if (totals == []):
            totals = get_compartment_dist(model, observable_name, output)
        else:
            totals += get_compartment_dist(model, observable_name, output)
    return totals

def plot_combined_dist(model, observable_name, output_array):
    totals = combine_dists(model, observable_name, output_array)
    plot_dist(totals)

def plot_dist(data):
    plt.figure()

    #end_counts = get_compartment_dist(model, observable_name, output)

    bins = range(int(min(data)), int(max(data)+1))
    print "end_counts: "
    print data

    total_count = sum(data)
    print "total membrane-bound tBid: "
    print total_count

    h, bins = histogram(data, bins=bins)
    h = h / float(len(data))
    print "bins: "
    print bins
    print "h: "
    print h
    print "sum h: %f" % sum(h)
    #bar(bins[0:len(bins)-1], h)
    plt.plot(bins[0:len(bins)-1], h)

    mean_counts = np.mean(data)
    print "mean counts: "
    print mean_counts

    poisson_counts = poisson.pmf(bins, mean_counts)
    #poisson_counts = poisson.pmf(bins, ^)
    print "sum poisson:  %f" % sum(poisson_counts)
    print poisson_counts
    plt.plot(bins, poisson_counts)


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
