from pysb import *
from numpy import array, recarray, mean, std
from matplotlib import pyplot as plt
from tbidbaxlipo.util import color_iter 
from scipy.stats import poisson
from pysb import bng
from pysb import kappa
from tbidbaxlipo.models import core

class Builder(core.Builder):

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

        # The compartment list needs to be built before declaring monomers
        self.cpt_list = ['c%d' % cpt_num for cpt_num in
                         range(1, int(self['Vesicles_0'].value)+1)]
        self.cpt_list.append('solution')

        # Now, monomers can be declared
        self.declare_monomers(self.cpt_list)

        # Declare initial conditions
        self.initial(self['tBid'](bh3=None, loc='c', cpt='solution'),
                     self['tBid_0'])
        self.initial(self['Bax'](bh3=None, a6=None, loc='c', cpt='solution'),
                     self['Bax_0'])

    def declare_monomers(self, cpt_list):
        """The monomer declarations need to override the default from the base
        class because each monomer has to have a list of possible compartments
        associated with it."""

        self.monomer('tBid', ['bh3', 'loc', 'cpt'],
                {'loc': ['c', 'm'], 'cpt':cpt_list}),
        self.monomer('Bax', ['bh3', 'a6', 'loc', 'cpt'],
                {'loc': ['c','m', 'i', 'p'],
                 'cpt': cpt_list})
        self.monomer('Pores', ['cpt'], {'cpt': cpt_list})

    # Translocation
    def translocate_Bax(self):
        print("site_cpt: translocate_Bax()")
        Bax_transloc_kf = self.parameter('Bax_transloc_kf', 1e-2,
                                         factor=(1/float(self.scaling_factor)))
        Bax_transloc_kr = self.parameter('Bax_transloc_kr', 1e-1)

        Bax = self['Bax']

        for cpt_name in self.cpt_list:
            if (not cpt_name == 'solution'):
                self.rule('Bax_translocates_sol_to_%s' % cpt_name,
                     Bax(loc='c', cpt='solution', bh3=None, a6=None) >>
                     Bax(loc='m', cpt=cpt_name, bh3=None, a6=None),
                     Bax_transloc_kf)
                self.rule('Bax_translocates_%s_to_sol' % cpt_name,
                     Bax(loc='m', cpt=cpt_name, bh3=None, a6=None) >>
                     Bax(loc='c', cpt='solution', bh3=None, a6=None),
                     Bax_transloc_kr)

        # Observables
        self.observable('cBax', Bax(loc='c', cpt='solution'))
        self.observable('mBax', Bax(loc='m'))
        for cpt_name in self.cpt_list:
            if (not cpt_name == 'solution'):
                self.observable('Bax_%s' % cpt_name, Bax(cpt=cpt_name))

    def translocate_tBid(self):
        print("site_cpt: translocate_tBid()")

        tBid_transloc_kf = self.parameter('tBid_transloc_kf', 1e-1,
                                         factor=(1/float(self.scaling_factor)))
        tBid_transloc_kr = self.parameter('tBid_transloc_kr', 0)

        tBid = self['tBid']

        for cpt_name in self.cpt_list:
            if (not cpt_name == 'solution'):
                self.rule('tBid_translocates_sol_to_%s' % cpt_name,
                     tBid(loc='c', bh3=None, cpt='solution') >>
                     tBid(loc='m', bh3=None, cpt=cpt_name),
                     tBid_transloc_kf)
                self.rule('tBid_translocates_%s_to_sol' % cpt_name,
                     tBid(loc='m', bh3=None, cpt=cpt_name) >>
                     tBid(loc='c', bh3=None, cpt='solution'),
                     tBid_transloc_kr)

        # Observables
        self.observable('ctBid', tBid(loc='c', cpt='solution'))
        self.observable('mtBid', tBid(loc='m'))
        for cpt_name in self.cpt_list:
            if (not cpt_name == 'solution'):
                self.observable('mtBid_%s' % cpt_name, tBid(cpt=cpt_name))

    def translocate_tBid_Bax(self):
        print("site_cpt: translocate_tBid_Bax()")
        self.translocate_tBid()
        self.translocate_Bax()

    # Binding
    def tBid_activates_Bax(self, bax_site='bh3'):
        print("site_cpt: tBid_activates_Bax(bax_site=" + bax_site + ")")

        # Forward rate of tBid binding to Bax (E + S -> ES)
        tBid_mBax_kf = self.parameter('tBid_mBax_kf', 1e-2)
        # Reverse rate of tBid binding to Bax (ES -> E + S)
        tBid_mBax_kr = self.parameter('tBid_mBax_kr', 1.5)
        # Dissociation of tBid from iBax (EP -> E + P)
        tBid_iBax_kc = self.parameter('tBid_iBax_kc', 1e-1)

        # Create the dicts to parameterize the site that tBid binds to
        bax_site_bound = {bax_site:1}
        bax_site_unbound = {bax_site:None}

        tBid = self['tBid']
        Bax = self['Bax']

        for cpt_name in self.cpt_list:
            if (not cpt_name == 'solution'):
                self.rule('tBid_Bax_bind_%s' % cpt_name,
                     tBid(loc='m', bh3=None, cpt=cpt_name) +
                     Bax(loc='m', cpt=cpt_name, **bax_site_unbound) >>
                     tBid(loc='m', bh3=1, cpt=cpt_name) %
                     Bax(loc='m', cpt=cpt_name, **bax_site_bound),
                     tBid_mBax_kf)
                self.rule('tBid_Bax_unbind_%s' % cpt_name,
                     tBid(loc='m', bh3=1, cpt=cpt_name) %
                     Bax(loc='m', cpt=cpt_name, **bax_site_bound) >>
                     tBid(loc='m', bh3=None, cpt=cpt_name) +
                     Bax(loc='m', cpt=cpt_name, **bax_site_unbound),
                     tBid_mBax_kr)

        # tBid dissociates from iBax
        self.rule('tBid_unbinds_iBax',
             tBid(loc='m', bh3=1) % Bax(loc='m', **bax_site_bound) >>
             tBid(loc='m', bh3=None) + Bax(loc='i', **bax_site_unbound),
             tBid_iBax_kc)

        # Observables
        # -----------
        self.observable('iBax', Bax(loc='i'))
        self.observable('tBidBax', tBid(loc='m', bh3=1) %
                                   Bax(loc='m', **bax_site_bound))

        for cpt_name in self.cpt_list:
            if (not cpt_name == 'solution'):
                self.observable('iBax_%s' % cpt_name,
                                Bax(loc='i', cpt=cpt_name))
                self.observable('tBidBax_%s' % cpt_name,
                                tBid(loc='m', bh3=1, cpt=cpt_name) %
                                Bax(loc='m', cpt=cpt_name, **bax_site_bound))

    # Pore formation
    def pores_from_Bax_monomers(self, bax_loc_state='i', reversible=False):
        raise NotImplementedError()
        print("site_cpt: pores_from_Bax_monomers()")

        pore_formation_rate_k = self.parameter('pore_formation_rate_k', 1e-3)

        Bax = self['Bax']
        Pores = self['Pores']

        # Compartment dependent
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
        for cpt_name in self.cpt_list:
            if (not cpt.name == 'solution'):
                self.observable('pores_%s' % cpt_name, Pores(cpt=cpt_name))

    def pores_from_Bax_dimers():
        raise NotImplementedError()

    def rate_calcs(nm_rate):

        vol = 60e-6
        print "Assumed volume: %g L" % vol

        nm_conc = 100 # assume nM
        molar_conc = nm_conc * 1e-9
        print "Example concentration: %g Molar" % molar_conc

        A = 6.022e23 # Avogadro's number
        total_n = molar_conc * vol * A
        print "Total molecules at molar conc of %g in vol of %g: %g" % (molar_conc, vol, total_n)

        # To convert nanomolar rate to molar rate, multiply the
        # nanomolar rate by 1e9
        print "deterministic (nm-1 s-1) rate: %g" % nm_rate # assume nm^-1 s^-1
        molar_rate = nm_rate * 1e9
        print "deterministic (M-1 s-1) rate: %g" % molar_rate
        molec_rate = molar_rate / (vol * A)
        print "molecular (molec^-1 s^-1) rate: %g" % molec_rate

        scaled_n = 1000
        print("scaled number of molecules meant to represent conc of %g Molar: %g " \
                % (molar_conc, scaled_n))
        actual_to_scaled = scaled_n / float(total_n)
        print "stochastic scaling factor of %g" % actual_to_scaled
        scaled_molec_rate = molec_rate / actual_to_scaled
        print "scaled molec rate: %g" % scaled_molec_rate

    # RUNNING THE MODEL
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
                for cpt in self.cpt_list if not cpt == 'solution']

    def run_model(self, tmax=12000, num_sims=1, use_kappa=True,
                  figure_ids=[0, 1]):
        xrecs = []   # The array to store the simulation data
        dr_all = []  # TODO: Delete this

        # Run multiple simulations and collect data
        for i in range(0, num_sims):
            # Run simulation using Kappa:
            if use_kappa:
                ssa_result = kappa.run_simulation(self.model, time=tmax,
                                                  points=100,
                                                  output_dir='simdata')
                xrecs.append(ssa_result)
            # Run simulation using BNG SSA implementation:
            else:
                ssa_result = bng.run_ssa(self.model, t_end=tmax, n_steps=100,
                                         cleanup=True)
                xrecs.append(ssa_result)
                #dr_all.append(get_dye_release(model, 'pores', ssa_result))

        # Convert the multiple simulations in an array...
        xall = array([x.tolist() for x in xrecs])

        # ...and calculate the Mean and SD across the simulations
        x_std = recarray(xrecs[0].shape, dtype=xrecs[0].dtype, buf=std(xall, 0))
        x_avg = recarray(xrecs[0].shape, dtype=xrecs[0].dtype,
                         buf=mean(xall, 0))

        # Plotting parameters, aliases
        ci = color_iter()
        marker = 'x'
        linestyle = '-'
        tBid_0 = self['tBid_0']
        Bax_0 = self['Bax_0']

        # Translocation: plot cyto/mito tBid, and cyto/mito Bax
        plt.ion()
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

        # Activation: plot iBax and tBidBax
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
        #     stoch',
        #        marker=marker, linestyle=linestyle)
        #xlabel("Time (seconds)")
        #ylabel("Pores/vesicle")
        #title("Pores/vesicle")
        #legend()

        #xlabel("Time (seconds)")
        #ylabel("Dye Release")
        #title("Dye release calculated via compartmental model")

        return xrecs[0]

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
