Properties of unimolecular dye release models
=============================================

The simplest possible model of vesicle permeabilization, one first proposed
by Schwarz ([Schwarz1992b]_), consists of the following elementary reactions:

.. math::

    P_c + L \rightleftharpoons P_l

.. math::

    m \cdot P_l \rightarrow Pore

Where :math:`P_c` and :math:`P_l` represent the peptide or protein in its
aqueous (cytosolic) or lipid-associated state, :math:`L` represents the
concentration of vesicles, and :math:`m` is the number of monomers involved
in pore formation.

:ref:`As mentioned earlier <previously_published_models>`, Schwarz's analysis
indicates that when pore formation is independent, then pores by definition
have a Poisson distribution across liposomes, which allows the experimental
efflux data to be used to calculate the average number of pores per vesicle
(and hence also the total number of pores, given that the concentration of
vesicles is known).

Assuming that the peptide/protein :math:`P` comes rapidly to steady-state
in its association with vesicles, pore formation proceeds linearly with
a rate proportional to :math:`P_l`, that is

.. math::

    \frac{d\ Pore}{dt} = k_{pore} [P_l]^n


Where the exponent :math:`n \leq m`. As a function of time, this will be
observed simply as a straight line for pore formation over time:

.. math::

    Pore(t) = k_{app} t

With the Poisson assumption of the independence of pore formation, this means
that dye release will follow a single exponential,

.. math::

    \mathrm{Frac. permeabilized} = 1 - e^{-k_{app}t}

When the data are well-fit by a single exponential (e.g. as they are in the
case of detergent-activated, truncated Bax, [Saito2000]_), then the scaling of
the apparent rate of pore formation, :math:`k_{app}`, as a function of the
concentration of total :math:`P`, can be analyzed, yielding estimates of the
values of the elementary parameters :math:`k_{pore}` and :math:`n`.

However, it is important to note the assumptions that this type of simple model
depends on:

**First, in its simplest form, it suggests that liposomes do not saturate with
protein.** That is, there is no limit to the amount of protein that can bind to
them.  In his work Schwarz addressed this by formulating a thermodynamic
relationship determining the amount of protein at the membrane at equilibrium
([Rizzo1987]_, [Schwarz1992b]_). In his formulation, a dimensionless partition
coefficient is multiplied by a thermodynamic interaction parameter which is
meant to take into account possible unfavorable interactions between proteins
as they accumulate in the membrane (e.g., repulsion of charges); this is
multiplied by the amount of lipid. The equation is

.. math::

    r_1 = \frac{\Gamma_1 \cdot c_1}{\alpha_1}

where :math:`r_1` denotes the resulting protein to lipid molar ratio (which
determines :math:`P_l` in our notation), :math:`\Gamma_1` is the dimensionless
partition coefficient, :math:`\alpha_1` is the thermodynamic interaction
coefficient, and :math:`c_1` is the concentration of free (aqueous) monomeric
peptide. A value of ` for :math:`\alpha_1` indicates that there is no
unfavorable interactions among proteins on the membrane, and hence binding is
not impaired. Schwarz suggests that peptides used at sufficiently low
concentrations will have :math:`\alpha_1` values close to 1.

**Second, it requires that the amount of membrane-bound protein be small
relative to the total protein, that is,** :math:`P_l << P_{total}`. If this is
not the case, then the kinetics will slow down as a substantial amount of the
aqueous protein available for partitioning to and permeabilizing and membranes
is depleted. This would be particularly problematic for permeabilization
processes in which the pore-forming agent becomes irreversibly associated with
membranes after permeabilization, as may be the case with Bax and some pore
forming peptides that are internalized into the inner leaflet of the bilayer
after pore formation ([Pokorny2002]_). 

Non-saturable partitioning of Bax to membranes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To illustrate the application of our modeling approach we begin by examining
the simple unimolecular pore formation models described by Schwarz.  The
simplest possible model, implemented in the function
:py:meth:`tbidbaxlipo.models.core.Builder.build_model_bax_schwarz`, consists of
Bax partitioning to membranes (as described above) followed by pore formation
from monomers::

    Bax(loc='c') ** solution + Vesicles() ** solution >>
    Bax(loc='m') ** ves + Vesicles() ** solution

    Bax(bh3=None, a6=None, loc='m') ** ves >>
    Bax(bh3=None, a6=None, loc='c') ** solution

    Bax(loc='m') >> Bax(loc='p') + Pores()]

First we examine how Bax binds or partitions to liposomes in this model and
others of its type. Here the assumption, as in Schwarz's model, is that Bax
simply partitions between the lipid and aqueous phases. This approach to
modeling the binding of Bax to liposomes has three important features.

First, even though in this simple model the binding capacity of the
vesicles is unlimited, if the amount of vesicles is too small (smaller than the
K_d of Bax for vesicles), then very little Bax will be bound. This can be seen
in this liposome titration plot:

.. plot::

    from tbidbaxlipo.plots.pore_plots import plot_liposome_titration
    plot_liposome_titration()

That is, **the fraction of Bax bound is dependent on the amount of liposomes in
this system and the affinity of Bax and liposomes for each other.** This model
should be identical to the Schwarz model with the exception that it does not
incorporate the possibility of unfavorable interactions between Bax monomers
that would hinder translocation and limit the vesicles' binding capacity.

This brings us to the second point: in this model the fraction of Bax bound is
determined only by the amount of liposomes in the system. In terms of the
familiar binding isotherm, the fraction of bound Bax is given by:

.. math::

    \frac{Bax_{bound}}{Bax_{total}} = \frac{Lipos}{K_D + Lipos}

Now, in a typical protein-protein binding context, we would note that the
variable :math:`Lipos` in the above expression refers to the amount of free
(unbound) liposomes at equilibrium, not to total liposomes. However, since
this model assumes that liposomes have unlimited binding capacity, this is
a moot point--the amount of liposomes free to bind Bax is identical to the
total amount of liposomes. **Thus for any amount of liposomes, the fraction of
Bax bound is determined only by the amount of liposomes, not by the amount
of Bax.**

This has an important consequence, namely that as Bax concentration is
increased, this will result in a proportional increase in the amount of Bax
per liposome. By manipulating the above expression we see that:

.. math::

    Bax_{bound} = Bax_{total} \frac{Lipos}{K_D + Lipos}

And so the amount of bound Bax increases linearly with total Bax, with a slope
determined by the :math:`K_D` and the amount of liposomes.

This raises the final point, **that the amount of Bax per liposome is given by a
Poisson distribution with an average given by the amount of Bax bound divided
by the amount of liposomes (I haven't checked this numerically).**

**Next we discuss dye release.** In this model **even if the amount of protein
and vesicles is equal, and it only takes one protein to form a pore, you won't
get 100% dye release. You will, however, reach 1 pore ver vesicle (on average)
at completion.** This is due to the uneven (Poisson) distribution of the pores
among liposomes. The ``one_cpt`` model actually captures this quite nicely, and
its results are validated by the ``n_cpt`` model, as shown in the figure
below:

.. plot::

    from tbidbaxlipo.plots.bax_heat_stoch_det_comparison import plot
    plot()

Of course, this also means that whenever liposomes are in excess on a molar
basis, it is impossible to get 100% permeabilization. This is because in this
model one can get a maximum of one pore per protein via an irreversible
process.

For example, if we set the concentration of both Bax and Vesicles to 50 nM,
we see that dye release plateaus at around 60%, whereas the average
number of pores per vesicle reaches completion at 1. This is because some
vesicles have more than one pore, whereas others have none:

.. plot::

    from tbidbaxlipo.models.one_cpt import Builder
    from tbidbaxlipo.plots.pore_plots import plot_pores_and_efflux
    from matplotlib import pyplot as plt
    params_dict = {'Bax_0': 50., 'Vesicles_0': 50.}
    b = Builder(params_dict=params_dict)
    b.build_model_bax_schwarz()
    plot_pores_and_efflux(b.model)
    plt.title('Dye release/pores for equimolar Bax and vesicles')

Now we look at the scaling of the **pore formation rate (not dye release rate)
as a function of Bax concentration.** In these plots the concentration of
liposomes is 5 nM, so at the maximum Bax concentration of 100 nM the maximum
achievable number of avg. pores is 20. This model produces a rate-law plot with
a straight line in the log-log plot with slope 1. Put in words, this means that
the **velocity of pore formation increases linearly with the amount of Bax,
never reaching saturation.** Moreover, this means that the total number of
pores that can be produced is equal to the total amount of Bax divided by the
number of Bax molecules required to form a pore. If pores are monomeric, then
there can maximally be as many pores as Bax molecules--steady state in the pore
timecourse will occur at this value. This means that if Bax concentration is
doubled, the steady state number of pores (and the rate) will double as well.

.. plot::

    from tbidbaxlipo.models.one_cpt import Builder
    from tbidbaxlipo.plots.pore_plots import plot_bax_titration
    b = Builder()
    b.build_model_bax_schwarz()
    plot_bax_titration(b.model)

Third, **this reaction scheme can be thought of as simple enzyme-substrate
catalysis where the enzyme, rather than the substrate is consumed.** Bax is the
enzyme, the liposome is the substrate, and the product is the permeabilized
liposome.  That is, it is: ``E + S <-> ES --> EP``. As such, the reaction must,
by necessity, always stop (or rather, asymptotically decelerate); it stops in
the limit when all ``E`` is consumed and all possible pores have been formed.
If the P/L ratio is high (>> 1) then dye release may become experimentally
indistinguishable from 100% well before the reaction is completed in terms of
pore formation. When P/L is high, the kinetic curve for the pores/ves velocity
appears as a straight line for the course of the experiment. When P/L is low,
the protein is rapidly consumed and both dye release and pores/ves plateau
quickly.

If the partitioning of protein to liposomes is fast (as it is expected to be),
then :math:`ES` comes rapidly to steady-state. In this model :math:`S`, the
liposomes, can never be diminished because more pores can always form, hence
this aspect of the Michaelis-Menten assumption applies.

**Fourth, unlike in the reversible model (see below) there can be no linear,
constant phase in the pores/ves plot for this model.** This would require a way
to form pores which did not continue to consume protein.

Reversible pore formation
~~~~~~~~~~~~~~~~~~~~~~~~~

The next case to consider is the same simple model as above but with the
modification that the proteins involved in pore formation can dissociate from a
vesicle and return to solution. If this is the case then a single protein can
permeabilize a (potentially large) number of vesicles.

The reverse rate dramatically effects the shape of the kinetic curves.
In the plot below a series of traces for pores per vesicle and percent dye
release are shown (in each case, as above, both Bax and vesicles are set
to concentrations of 50 nM as shown above for the irreversible case).

.. plot::

    # 50nM Vesicles and Bax, pore formation forward rate of 1e-3
    from tbidbaxlipo.plots.pore_plots import \
         plot_effect_of_pore_reverse_rate
    plot_effect_of_pore_reverse_rate()

As the plot shows, if the reverse rate is slow (1e-6), the pore formation
process is very similar to the irreversible case, in which the pores per
vesicle curve plateaus at 1.

When the pore reverse rate is fast (1e-2), the protein is returned to the
solution essentially immediately after the pore is formed, allowing it to
permeabilize other liposomes. In this case the conversion of liposomes
follows the reaction scheme

.. math::

    E + S \rightleftharpoons ES \rightarrow EP \rightarrow E + P

in which :math:`E` is Bax, :math:`S` is the unpermeabilized liposome, and
:math:`P` is the permeabilized liposome. :math:`EP` is the state in which
Bax remains bound to the liposome after permeabilizing it. However, if the
rates of the pore formation and pore reversal processes are fast (to be defined
formally later) the quantities of :math:`E` and :math:`ES` are relatively
undiminished, and the conversion of :math:`S` to :math:`P` is approximately
a first-order process with a rate proportional to :math:`E`:

.. math::

    S \rightarrow P

However, since in this case the "substrate" :math:`S`, the liposomes, is not
consumed by pore formation, the formation of the product :math:`P` is actually
linear (zero order). This can be seen in the plot as a straight-line velocity
of pore formation for the fast reverse rate.

In the third case, the reverse rate occupies an intermediate value, such that
a significant, and constant, amount of protein :math:`E` is occupied on
permeabilized liposomes.

Saturable Bax Binding
~~~~~~~~~~~~~~~~~~~~~

Next we examine the case where the binding of Bax to liposomes is saturable,
that is, there is a limited number of binding sites on liposomes for Bax.

First we look at the fraction of Bax bound as a function of Bax for simple
partitioning vs. a model in which the finite nature of liposome binding sites
is explicitly accounted for. As discussed above, for the partitioning model,
the fraction of Bax bound is determined only by the amount of liposomes,
whereas in the binding site model, the fraction of Bax bound decreases once the
liposomes become saturated and none of the additional Bax can bind. In the
simulations shown below there is 30 nM of liposomes or liposome "binding
sites".

.. plot::

    from matplotlib import pyplot as plt
    from tbidbaxlipo.models import lipo_sites, one_cpt
    from tbidbaxlipo.plots.pore_plots import plot_fraction_bax_bound

    plt.ion()
    params_dict = {'Vesicles_0': 30}
    b = lipo_sites.Builder(params_dict=params_dict)
    b.translocate_Bax()
    plot_fraction_bax_bound(b.model, figure_id=10)

    b = one_cpt.Builder(params_dict=params_dict)
    b.translocate_Bax()
    plot_fraction_bax_bound(b.model, figure_id=10)
    plt.legend(['Binding site', 'Partitioning'], loc='lower left')

Next we examine the behavior of this model upon incorporating pore formation,
simulating the pore formation timecourse for many Bax concentrations as above.
What these plots show is that not only does the steady-state (maximal) value
for the number of pores saturate with increasing Bax, but the initial velocity
saturates as well. Rather than having a slope of 1 as in the partitioning model,
the log-log plot starts out with a slope of 1 and then saturates.

.. plot::

    from tbidbaxlipo.models import lipo_sites
    from tbidbaxlipo.plots.pore_plots import plot_bax_titration
    params_dict = {'Vesicles_0': 2, 'pore_formation_rate_k':5e-3}
    b = lipo_sites.Builder(params_dict=params_dict)
    b.build_model_bax_schwarz()
    plot_bax_titration(b.model)

The other thing that this plot shows is that at saturation, all curves reach
a final value of 1 pore per "liposome" on average; however in this model the
liposomes really represent liposome binding sites. The reason why the value of
1 is always attained is because once a pore forms, the liposome binding site
remains irreversibly "bound" to the Bax pore. Schematically, this is

.. math::

    E + S \rightleftharpoons ES \rightarrow EP

Because the binding between Bax and liposome binding sites is 1 to 1, there can
ever be as many pores as there are molecules of EP, and hence as many molecules
of S. Thus the average number of pores per site (total pores divided by number
of sites) is 1.

The scheme above also shows that the reaction slows down at late times due not
only to the consumption of S (unbound liposome binding sites) but also due to
the consumption of E (free, non-pore Bax).

A dimeric Bax pore
~~~~~~~~~~~~~~~~~~

I NEED TO REVISIT ALL OF THIS ANALYTICALLY TO MAKE SURE IT IS NOT THE
RESULT OF NUMERICAL ARTIFACTS.

Changing the model to use a dimeric pore has one obvious consequence--the
average number of pores per vesicle, and hence the total number of pores goes
down by half.

But there is another interesting consequence--in the Bax
titration, the slope of the log-log plot starts out at 2 for low concentrations
of Bax, then shifts to 1 at high concentrations of Bax!

Change the rate of dimerization changes this--the rate limiting step is
dimerization only when dimerization is slow relatively to the other processes.
Changing the dimerization rate to be fast makes the log-log slope approach 1.
Notably, when Bax concentrations are low relative to the dimerization rate,
the rate limiting step again becomes dimerization.

Conversely, when the dimerization forward rate is made to be very low, the
slope of the log-log plot is two.

This is true even when the reverse rate for dimerization is 0, so the issue is
not one of the Kd of dimerization, but rather the bimolecularity of the
interaction.

It's not totally clear to me how when the dimerization rate is fast the order
of the rate law is still 1, even though twice the amount of Bax is still
required to permeabilize the same amount of liposomes. I suppose that that
aspect is irrelevant--that is a constant factor of a change (2-fold), but
doesn't speak to the exponent of the rate law. The exponent of the rate law
refers to the order of the rate-limiting reaction, i.e., how the rate scales
with concentration. So if the rate is linear in a dimer of :math:`E`, that is
still a log-log slope of 1 for the rate, even though it scales with the dimer,
not the monomer.

Presumably the bend in the curve comes at around the point where the average
Bax per liposome is around 1?

A tetrameric Bax pore
~~~~~~~~~~~~~~~~~~~~~

Interestingly, using a scheme in which Bax pores consist of tetramers that
assemble by dimerization of dimers, the log-log rate law plot (for the
initial rate) has a slope of 4 four low concentrations of Bax, but this
then bends down to 1 at high concentrations of Bax.

You can get some funny results switching the log-log slope between 2 and 4
depending on the parameters you choose, but to some extent this depends
on the numerical sampling done to get the initial slope. Ideally, it would
be done with the same fitting equation as used for the real data.

This is all very confusing. Ideally I would do this analytically using
perturbation theory.

Other analyses to do
~~~~~~~~~~~~~~~~~~~~

* Do analysis for trimeric vs. dimeric pores, see if they give 3/2 rate laws,
  respectively
* Do with cooperative assembly
* Do with stepwise assembly

Almeida
-------

This model always goes to 100% permeabilization. However, it should be noted
that it was developed specifically to compare all-or-none vs. graded
forms of dye release.

Newmeyer
--------

This model also always goes to 100% permeabilization, even though many of the
authors' own plots show otherwise.

Schlesinger
-----------

Does the assumption about the rate law hold in this case?

