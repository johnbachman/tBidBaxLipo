Kinetic analysis of Bax pore formation
======================================

Possible title: "Bax recruitment is kinetically limited by Bid, while Bid
recruitment is stoichiometrically limited by membranes"

(Aiming for 5 figs, 3.25 x 4 in each, 50k characters = ~10 pages)
Page length estimator: http://www.biophysics.org/tabid/556/default.aspx

Mitochondrial outer membrane permeabilization, or MOMP, is widely considered to
be the key decision point in apoptosis. MOMP is regulated by the Bcl-2 family
of proteins. Evasion of apoptosis is a hallmark of cancer and can occur by
aberrant Bcl-2 protein levels (ref?).

Extensive work has elucidated many of the structural features of the
interactions of these proteins with membranes (Annis, Westphal?) and with each
other (Czabotar, Walensky). This knowledge of structure has led to the
development of peptides and small-molecules that can trigger or sensitize cells
to commit apoptosis in a targeted fashion (Walensky, Gavathiotis, ABT 199,
263).

However, the kinetic properties of this coordinated system of proteins, and
their origin in structure, is relatively less understood. Particularly notable
is that measured levels of Bcl-2 proteins in cancer cells vary over multiple
orders of magnitude in multiple dimensions. Natural variation in these proteins
in cells define a phase space which determines the apoptotic sensitivity
of the cells. Exploring this phase space using mechanistic in vitro studies helps
us to better understand different mechanistic regimes and rate limiting steps.

Permeabilization kinetics by Bax have been observed both in cells and in vitro
biochemical systems with synthetic liposomes. Observations in single cells have
shown that after delays of up to several hours, the release of mitochondrial
proteins occurs over 1-3 minutes in a highly coordinated fashion (refs); at
faster timescales a accelerating "wave" can be observed (MOMP waves).
Explanations offered for the propagation include Bax auto-activation (Andrews
ref), passive translocation of Bid and Bax between membranes (Edlich,
Shamas-Din, Gilmore-FAK), and the release of reactive oxygen species (ref). 

Reconstituted biochemical systems have been used effectively in the Bcl-2 field
to discover and describe activity (Andrews, Newmeyer, Kuwana, etc.)

Almeida and coworkers have used detailed mathematical models to describe the
mechanism of action of various pore forming peptides (refs). More recently,
Kushnareva et al. Notably, the kinetics of pore formation by Bax have been
described as increasing cooperatively with Bax concentration
(Saito/Schlesinger), linearly (Kushnareva), and in a saturating fashion
(Satsoura), likely due to differences in the experimental systems used and the
analytical frameworks applied.

Bax goes through a series of conformational states. Activation can occur
through interaction with a BH3-only protein such as cBid or with other physical
perturbations such as elevated temperature (Green), incubation with detergent
(Hsu/Youle) or pH (???).

Previous efforts to characterize the thermodynamic and kinetic properties of
membrane permeabilization by pore forming peptides have focused primarily on
pore forming peptides. However, a handful of studies have also focused on the
Bcl-2 proteins Bid and Bax.

This represents the first attempt to quantitatively account for the activity of
multiple interacting proteins in membranes. Moreover, we define our differing
hypotheses explicitly in terms of kinetic mechanisms and rigorously evluate
them against the data.

Results
-------

Bax membrane insertion kinetics saturates at high Bax concentrations.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(Alternatively, start with the 2-D dose response: Bax and liposomes,
Bid constant).

(Or should I start with permeabilization?)

A two-dimensional titration of Ba

To study the kinetics of Bax insertion we used a single-cysteine mutant of Bax
labeled with the environment-sensitive dye NBD as previously described
(Lovell). The NBD at the labeled position on the protein is weakly fluorescent
in an aqueous (polar) environment, but fluoresces more strongly in a membrane
(hydrophobic) environment. We observed that while holding the concentration of
liposomes and cBid constant, titrating the concentration of unlabeled Bax
slowed the fractional recruitment of the NBD-Bax (Figure 1A).

.. figure:: ../doc/_static/130911_c126_bax_titration_9.png
    :width: 6in
    :align: center

    **Figure 1A**

A previous study of Bax-membrane interactions revealed that Bax binding to
membranes unexpectedly saturated at high Bax/lipid ratios (Satsoura). This
study involved measuring Bax localization by a variety of methods including gel
filtration and FCS after a 2 hour incubation period. One of the explanations
proposed included the possibility of limited "binding sites" on the liposome
for Bax. Another possibility is that the saturation of Bax recruitment is
primarily a kinetic phenomenon, and that if allowed to reach equilibrium Bax
would reach levels of recruitment independent of protein concentration.

To estimate the kinetic rate and an extrapolated equilibrium value, we fit to
the two-parameter, single-exponential equation shown in XXX. Fitted curves are
shown as the black lines in Figure 1A. In this equation the parameter
:math:`F_{max}` gives an estimated equilibrium F/F0 fluorescence, while
:math:`k_1` gives an estimate of the recruitment rate, in units of sec\
:sup:`-1`.  The fitted values show that the recruitment rate decrease in a
dose-dependent fashion with Bax concentration (Figure 1B), while the predicted
fraction of Bax inserted is invariant to Bax.

.. figure:: ../doc/_static/130911_c126_bax_titration_10.png
    :width: 6in
    :figwidth: 6in
    :align: center

    **Figure 1B**. Fitted k1 values vs. Bax concentration.

.. figure:: ../doc/_static/130911_c126_bax_titration_11.png
    :width: 6in
    :figwidth: 6in
    :align: center

    **Figure 1C**. Fitted Fmax values vs. Bax concentration.

To explain this data, we considered four possibilities, which we formulated as
mathematical models: 1) Bax recruitment and insertion is mediated by simple
partioning to the membrane phase as previously described for peptides (Schwarz,
Almeida?, Satsoura?); 2) Recruited Bax is dependent on liposome binding sites
or overall liposome binding capacity for its bound state 3) Bax recruitment is
dependent on a limited set of liposome sites for its peripherally bound but not
for its inserted state; 4) Bax recruitment is dependent on the activator Bid in
an enzymatic fashion; 5) Bax recruitment is mediated by dimerization with Bid,
with the two proteins able to bind after the activation of Bax (product
inhibition). We considered this latter possibility because we previously showed
that Bid and Bax remain bound after most Bax is activated, suggesting the
possibility that the accumulation of activated Bax:cBid complexes could retard
Bid's ability to recruit additional Bax. After fitting each of these models to
the underlying data (**Supplemental Figures**), the models were fit with
equation XXX and the fitted values for k1 and Fmax at each Bax concentration
were plotted (**Figures 1B and 1C**). Notably, Model 5 was unable to fit the
underlying data and was only poorly approximated by equation XXX, hence it was
excluded.

**As shown in XXX**, only the models 3 and 4, can fit the underlying kinetic
data. Model 1 predicts that insertion rate will not scale with Bax
concentration, as it clearly does. Model 2 predicts that the kinetics will stay
roughly the same while the equilibrium amount recruited will diminish. Models 3
and 4, both of which depend on the formation of a transient but saturable
complex between Bax and either a liposome or cBid, reproduce the observed data.

Bax recruitment kinetics depends strongly on liposome concentration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No saturation of rate observed in 43C heated Bax. However, FMax saturates
(obviously) near 100% permeabilized. Hence the initial rate would presumably
also saturate since they are composed.

cBid determines the rate, and Bax the extent, of membrane permeabilization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Experiment todo list
--------------------

* Repeat NBD-Bax titration at single Bid concentration (either with or without
  competitor) to show saturation. Do replicates to get error bars on Fmax and
  k1 values.

* Perform ANTS release with Bax titration to see if rate saturates? (can I use
  the data from 7/24 for this?)

* (already done?): WT DKO mitos, incubate with Bid and Bax; then pellet, and
  incubate with IMS-EGFP mitos; expectation is that they don't permeabilize
  much.

* Bax hole? Incubate lipos with Bax + Bid or BH3; then incubate along with
  another set of lipos plus additional Bax. Measure permeabilization of the
  second set of liposomes. Do you get less permeabilization of the second set
  when the first set carries Bax? If so, suggests that the Bax preferentially
  goes to the second set. Ideally this could be done with no activator so that
  all activation was due to Bax auto-activation. Or could treat with heat or
  peptides and do gel filtration.

* Bax hole expt. Incubate unlabeled Bax with unlabeled liposomes. Then incubate
  labeled Bax with mCherry lipos. Does the amount of FRET decrease on the amount
  of Bax put into the unlabeled lipos?

* As an extension, pre-incubate the "hole" lipos with K21E/BH3 mut Bax, or K21E
  activated with mutant Bim BH3 which should be able to autoactivate itself but
  not the complementary Bax.

* Does BclXL prevent the Bax hole phenomenon?

Cover letter
------------

We believe that this work represents a kind of "systems biochemistry" aiming at
characterizing key biochemical processes in quantitative mechanistic detail.


