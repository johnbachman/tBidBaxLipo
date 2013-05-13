Unimolecular Dye Release Models
===============================

General considerations: Schwarz
-------------------------------

Schwarz's analysis indicated that when pore formation is independent, then
pores by definition have a Poisson distribution across liposomes, which allows
the efflux data to be used to calculate the average number of pores per vesicle
(and hence also the total number of pores, given that the concentration of
vesicles is known).

In his analysis of the permeabilization kinetics of the peptide melittin,
Schwarz found that the kinetics at each concentration followed a two-phase
process which can be fit with a three-parameter equation.

Moreover, he noted that the log-log plot showed a straight-line fit indicating
a third-order dependence, while that of the final rate had a second-order
dependence (I should check these numbers TODO).

He modeled this using a simple model involving the partitioning of the
proteins to liposomes, followed by the irreversible formation of pores. There
are some notable features of his model and analysis:

* First, his predictions for the amount of protein partitioned to the membrane
  is based on a dimensionless partition coefficient, multiplied by a
  thermodynamic interaction parameter which is meant to take into account
  possible unfavorable interactions between proteins as they accumulate in
  the membrane (e.g., repulsion of charges), multiplied by the amount of
  lipid. This is described by the equation:

.. math::

    r_1 = \frac{\Gamma_1 \cdot c_1}{\alpha_1}

where :math:`r_1` denotes the protein to lipid molar ratio, :math:`\Gamma_1`
is the dimensionless partition coefficient, :math:`\alpha_1` is the
thermodynamic interaction coefficient, and :math:`c_1` is the concentration
of free (aqueous) monomeric peptide.

* Second, he notes that in all cases his assumption is that the amount of
  protein involved in oligomerization or pore formation is very small relative
  to the total amount, implying that the pool of protein is never appreciably
  consumed. This is a bit odd since one of his first observations about the
  melittin data is that the reaction appears to stop at levels that are
  significantly below 100% permeabilization.

* In an analysis which I haven't completely followed, he predicts that the
  oligomerization model which he describes predicts the third-order dependence
  but not the second-order dependence. He hypothesizes an initial deposit of
  dimers that is then consumed as the bilayer properties change rapidly. This
  might be worth trying to reproduce in the model.

The simplest model: monomeric Bax pores
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To illustrate the application of our modeling approach we begin by examining
thesimple unimolecular pore formation models described by Schwarz.  The
simplest possible model, implemented in the function
:py:meth:`tbidbaxlipo.models.core.Builder.build_model_bax_schwarz`, consists of
Bax partitioning to membranes followed by pore formation from monomers::

    Bax(loc='c') ** solution + Vesicles() ** solution >>
    Bax(loc='m') ** ves + Vesicles() ** solution

    Bax(bh3=None, a6=None, loc='m') ** ves >>
    Bax(bh3=None, a6=None, loc='c') ** solution

    Bax(loc='m') >> Bax(loc='p') + Pores()]

This simple model demonstrates a number of important principles. First, even
though in this simple model the binding capacity of the vesicles is unlimited,
if the amount of vesicles is too small (smaller than the K_d of Bax for
vesicles), then very little Bax will be bound. This can be seen in this
liposome titration plot:

.. plot::

    from tbidbaxlipo.pore_comparison.schwarz import plot_liposome_titration
    plot_liposome_titration()

This model should be identical to the Schwarz model with the exception that
it does not incorporate the possibility of unfavorable interactions between
Bax monomers that would hinder translocation and limit the vesicles' binding
capacity.

Second, even if the amount of protein and vesicles is equal, and it only takes
one protein to form a pore, you won't get 100% dye release. You will, however,
reach 1 pore ver vesicle (on average) at completion. This is due to the uneven
(Poisson) distribution of the pores among liposomes. The ``one_cpt`` model
actually captures this quite nicely, and its results are validated by the
``n_cpt`` model. Of course, this also means that whenever liposomes are in
excess, it is impossible to get 100% permeabilization. This is because in this
model one can get a maximum of one pore per protein via an irreversible
process.

For example, if we set the concentration of both Bax and Vesicles to 50 nM,
we see that dye release plateaus at around 60%, whereas the average
number of pores per vesicle reaches completion at 1. This is because some
vesicles have more than one pore, whereas others have none:

.. plot::

    from tbidbaxlipo.models.one_cpt import Builder
    from tbidbaxlipo.pore_comparison.schwarz import plot_pores_and_efflux
    from matplotlib import pyplot as plt
    params_dict = {'Bax_0': 50., 'Vesicles_0': 50.}
    b = Builder(params_dict=params_dict)
    b.build_model_bax_schwarz()
    plot_pores_and_efflux(b.model)
    plt.title('Dye release/pores for equimolar Bax and vesicles')

Second, this model does indeed reproduce a rate-law plot with a straight line
in the log-log plot with slope 1:

.. plot::

    from tbidbaxlipo.models.one_cpt import Builder
    from tbidbaxlipo.pore_comparison.schwarz import plot_bax_titration
    b = Builder()
    b.build_model_bax_schwarz()
    plot_bax_titration(b.model)

Third, this reaction scheme can be thought of as simple enzyme-substrate
catalysis where the enzyme, rather than the substrate is consumed. Bax is the
enzyme, the liposome is the substrate, and the product is the permeabilized
liposome.  That is, it is: ``E + S <-> ES --> EP``. As such, the reaction must,
by necessity always stop (or rather, asymptotically decelerate); it stops in
the limit when all ``E`` is consumed and all possible pores have been formed.
If the P/L ratio is high (>> 1) then dye release may become experimentally
indistinguishable from 100% well before the reaction is completed. When P/L is
high, the kinetic curve for the pores/ves velocity appears as a straight line
for the course of the experiment. When P/L is low, the protein is rapidly
consumed and both dye release and pores/ves plateau quickly.

Ff the partitioning of protein to liposomes is fast (as it is expected to be),
then :math:`ES` comes rapidly to steady-state. *In this model :math:`S`, the
liposomes, can never be diminished because more pores can always form,* hence
this aspect of the Michaelis-Menten assumption applies.

Fourth, unlike in the reversible model (see below) there can be no linear,
constant phase in the pores/ves plot for this model. This would require a way
to form pores which did not continue to consume protein.

Reversible pore formation
~~~~~~~~~~~~~~~~~~~~~~~~~

The next case to consider is where the proteins involved in pore formation can
dissociate from a vesicle and return to solution. If this is the case then
a single protein can permeabilize a (potentially large) number of vesicles.

The reverse rate dramatically effects the shape of the kinetic curves.
In the plot below a series of traces for pores per vesicle and percent dye
release are shown (in each case, as above, both Bax and vesicles are set
to concentrations of 50 nM as shown above for the irreversible case).

.. plot::

    # 50nM Vesicles and Bax, pore formation forward rate of 1e-3
    from tbidbaxlipo.pore_comparison.schwarz import \
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
a first-order decay process with a rate proportional to :math:`E`:

.. math::

    S \rightarrow P

In the third case, the reverse rate occupies an intermediate value, such that
a significant, and constant, amount of protein :math:`E` is occupied on
permeabilized liposomes.

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
** Do with cooperative assembly
** Do with stepwise assembly

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

