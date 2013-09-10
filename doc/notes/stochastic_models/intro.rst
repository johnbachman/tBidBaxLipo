.. _stochastic_models_intro:

Introduction
============

Modeling Apoptosis and MOMP
---------------------------

Apoptosis is a form of programmed cell death that plays an important role in
both normal development and disease. Apoptosis can be triggered following the
detection of internal cellular stress by the so-called "intrinsic" pathway, or
by the binding of ligands to death receptors via the "extrinsic" pathway.  In
both pathways the regulation of mitochondrial outer membrane permeabilization
(MOMP) is believed to be the key step in committing the cell to apoptosis
[Tait2010]_.

Mitochondrial outer membrane permeabilization is regulated by the Bcl-2 family
of proteins, which act at the mitochondrial outer membrane to create the lethal
pore that releases pro-apoptotic proteins from the mitochondrial intermembrane
space into the cytosol [Chipuk2010]_, [Youle2008]_. The activity of Bcl-2
proteins is regulated by post-translational modifications, translocation to and
from membranes, conformational changes, and homo- and hetero-oligomerization
among family members.

The regulation of MOMP has been studied in a variety of experimental model
systems, including live cells, permeabilized cells, mitochondria isolated from
cells, and synthetic liposomes. In addition, a number of mathematical models
have been developed to interpret data generated from the study of these
systems, including by our group [Spencer2011]_. Some of these models explicitly
account for the compartmentalization of Bcl-2 interactions at the membrane,
while others do not; in our own work, we have found that this
compartmentalization plays a key role in the regulation of Bcl-2 protein
activitity [Albeck2008b]_, consistent with recent biochemical findings
[Leber2010]_.

Even in models which explicitly treat the compartmentalization of Bcl-2
interactions at membranes, the membrane compartment is treated as an
undifferentiated, bulk compartment with reduced volume relative to the cytosol.
However, diffusion of proteins among mitochondrial membranes is restricted in
cases where the mitochondrial network has undergone significant fission; in
addition, sites of Bcl-2 protein binding may be further limited to specific
contact points between inner and outer membranes. This suggests that the
distribution of Bcl-2 family members across the mitochondrial network may be
non-homogeneous and subject to stochastic effects. However, the role of such
compartmentalization in governing the activity of Bcl-2 proteins, if any, is
unknown.

In addition to its possible relevance to physiological process of apoptosis
within the cell, compartmentalization of Bcl-2 protein interactions to many
discrete membrane compartments may play a significant role in commonly-used
model systems of MOMP, such as isolated mitochondria and synthetic lipid
vesicles. In these systems, a bulk solution of isolated mitochondria or
synthetic lipid vesicles is treated with purified Bcl-2 proteins.  While each
of the components is abundant, the amount of protein per compartment may be
quite small, affecting the quantitative interpretation of the assay results.

In this work we use synthetic lipid vesicles as a model system to study the
effects of reaction compartmentalization on permeabilization kinetics.
Synthetic lipid vesicles have been used for many years as a model system for
studying the activity of pore-forming peptides and proteins, and have resulted
in a number of key findings regarding the activity of Bcl-2 proteins.  In these
studies, lipid vesicles enclosing a permeabilization marker such as a
fluorescent dye are incubated with one or more Bcl-2 family members and the
change in fluorescence is measured over time.

Kinetics of permeabilization processes
--------------------------------------

The kinetics of permeabilization processes differ from canonical enzyme
kinetics in two important ways.

First, from the perspective of a single vesicle, what may be **a continuous
process** of protein interaction and pore accumulation on the vesicle only
**registers as a single event,** the formation of the first pore. This makes
mechanistic interpretation of in vitro dye release assays challenging because
the experimentally observed dye release kinetics must be related to the
kinetics of pore accumulation or protein interaction in some fashion.

Second, because the vesicles serve as discrete reaction compartments, the
process involves **an ensemble of distinct reactions,** which in some cases
will not follow kinetics described in terms of compartment totals or averages.
For example, if a permeabilization process involves two proteins `P1` and `P2`,
and `P1` is roughly equimolar with vesicles, then many vesicles will have no
molecules of `P1` and they may never be permeabilized. A model described in
terms of the total amount of membrane-bound `P1` would predict that permeabilization would proceed on all vesicles, albeit slowly.

Previously published models of permeabilization processes, discussed below,
mainly address the first of these challenges--the interpretation of dye release
kinetics in mechanistic terms. The second challenge, the deviation from
kinetics described in terms of compartment averages, has not previously been
addressed.

Here we develop explicit stochastic models of one or more membrane proteins
interacting on an ensemble of vesicles, and compare the results
from the stochastic model to the expected deterministic equivalent.  We find
that the deterministic approach is an effective representation of the
underlying stochastic process for the kinds of systems that have been studied
to date, namely single pore-forming agents such as cytolytic peptides. However,
we find that the deterministic model deviates from the stochastic model in many
cases of interest, such as multiple interacting membrane proteins that regulate
each other (e.g., Bcl-2 proteins), or even in the case of a single agent that
can self-recruit (as in the case of the Bcl-2 protein Bax). Investigation
of these important systems requires the application of new numerical and
analytical methods.

.. _previously_published_models:

Previously published models of permeabilization kinetics
--------------------------------------------------------

Gerhard Schwarz and colleagues developed an extensive theory of vesicle
permeabilization kinetics in a series of papers in the 1990s ([Schwarz1990]_,
[Schwarz1992a]_, [Schwarz1992b]_, [Schwarz1995]_). The starting point for
Schwarz's analysis is the assumption that the process of pore formation by
peptides is independent, such that throughout the process pores are distributed
across liposomes according to a Poisson distribution. Since the fraction of
unpermeabilized vesicles (that is, vesicles with no pores) is given by the
Poisson equation :math:`e^{-P(t)}`, where :math:`P(t)` is the average number of
pores per liposome at time `t`, one can infer the average number of pores as
:math:`P(t) = -\log E(t)` where :math:`E(t)` is the fraction of unpermeabilized
vesicles, given by the data.  This transformation allows models to be specified
in terms of the rate of pore formation and thus compared to the dye release
data. He goes on to suggest a simple mechanistic model of permeabilization
which can be solved analytically in which pores increase linearly (and hence
dye release follows a single exponential).

Almeida and colleagues ([Pokorny2002]_, [Gregory2008]_) adopt an alternative
approach in which the kinetics of dye release (not pore formation) are
explicitly modeled according to a set of coupled ordinary differential
equations which are solved numerically. In these equations the rate of dye
efflux is proportional to the product of the fraction of unpermeabilized
vesicles :math:`(1 - CF_{out})` and the amount of protein in the pore-competent
state, :math:`P^*` (adopting Almeida's notation):

.. math::

    \frac{dCF_{out}}{dt} = \frac{k_{eflx}}{v_o[L]}(1 - CF_{out})P^*

where :math:`v_o` and :math:`[L]` are constants. Note here that if the
pore-forming species :math:`P^*` is at steady state, then :math:`CF_{out}` will
increase according to an exponential function, corresponding to Schwarz's
simple model.

Schlesinger and colleagues ([Saito2000]_, [Schlesinger2006]_,
[Christenson2006]_), citing the earlier work of Schwarz, observes that dye
release follows a simple exponential function with a fairly small linear
component:

.. math::

    R(t) = A_1(1 - e^{\frac{t}{\tau}}) + mt

Using this equation, Schlesinger analyzes the concentration dependence of the
time constant :math:`\tau` for permeabilization by detergent-activated,
C-terminally truncated Bax, and concludes that the kinetics indicate a tetramer
as the minimal pore size [Saito2000]_.

Most recently, Kushnareva et al. ([Kushnareva2012]_) constructed a model in
which they treat the pore-forming species as an enzyme that converts vesicles
from full to empty. Given a constant amount of active pore former the dye
release kinetics therefore follow a single exponential, as the fraction of
intact vesicles "decays." The focus of their model is the kinetics of
activation of the pore forming species, which they term the "catalyst".

With the exception of [Kushnareva2012]_, these models were developed to model
the action of a single pore-forming agent, such as pore-forming peptides or
truncated Bax. In [Kushnareva2012]_, however, the model ostensibly takes into
account the action of the protein cBid as well as a hypothetical catalyst
protein.

All of these models **share the property that** insofar as they treat mechanism
explicitly, **they treat vesicles as a single undifferentiated compartment,
with deterministic rate equations for the interaction kinetics of
membrane-bound proteins.** These deterministic models are an approximation for
the underlying physical process in which proteins translocate to and from
discrete vesicles, interact with each other, and form pores. However, as
mentioned above, it is possible to imagine scenarios in which the deterministic
approximation might break down.

In the sections that follow we first describe our stochastic modeling approach
which allows us to explicitly describe the behavior of the stochastic system
and thereby explore cases in which the deterministic approximation breaks down.
In addition to testing the validity of the deterministic modeling approach,
this analysis helps us to develop better intuition for how to interpret assays
in cases in which the stochastic regime is dominant.

