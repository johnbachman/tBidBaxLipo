Unimolecular Dye Release Models
===============================

Schwarz
-------

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
  lipid. I am curious to see what a "binding curve" looks like in this
  model.

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

In any case, the simplest model, with partitioning to membranes followed by
pore formation from monomers, already demonstrates a few important principles.

First, even if the amount of protein and vesicles is equal, and it only takes
one protein to form a pore, you won't get 100% dye release. You will, however,
reach 1 pore ver vesicle (on average) at completion. This is due to the uneven
(Poisson) distribution of the pores among liposomes. The ``one_cpt`` model
actually captures this quite nicely, and its results are validated by the
``n_cpt`` model. Of course, this also means that whenever liposomes are in
excess, it is impossible to get 100% permeabilization. This is because in this
model one can get a maximum of one pore per protein via an irreversible
process.

Second, this model does indeed reproduce a rate-law plot with a straight line
in the log-log plot with slope 1:

.. plot::

    from tbidbaxlipo.pore_comparison.schwarz import plot_Bax_titration
    plot_Bax_titration()

Third, this reaction scheme can be thought of as simple enzyme-substrate
catalysis where the enzyme, rather than the substrate is consumed. Bax is the
enzyme, the liposome is the substrate, and the product is the permeabilized
liposome.  That is, it is: ``E + S <-> ES --> ES*``. As such, the reaction
must, by necessity always stop (or rather, asymptotically decelerate); it stops
in the limit when all ``E`` is consumed and all possible pores have been
formed. If the P/L ratio is high (>> 1) then dye release may become
experimentally indistinguishable from 100% well before the reaction is
completed. When P/L is high, the kinetic curve for the pores/ves velocity
appears as a straight line for the course of the experiment. When P/L is low,
the protein is rapidly consumed and both dye release and pores/ves plateau
quickly.

Fourth, there can be no linear, constant phase in the pores/ves plot for this
model. This would require a way to form pores which did not continue to
consume protein.

A dimeric Bax pore
~~~~~~~~~~~~~~~~~~

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

