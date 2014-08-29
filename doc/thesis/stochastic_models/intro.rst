.. _stochastic_models_intro:

Introduction
============

"Compartmentalization of reactions at membranes plays a critical role in the
rate and extent of apoptotic pore formation by Bax"

"Modeling the complex regulation of pore forming proteins by explicit
enumeration of vesicle reaction compartments."

* Permeabilization of membranes by peptides and proteins is important for
  bacterial toxins, pore forming peptides used as therapeutics, and in the
  processes of apoptosis and necroptosis. Also for prions (Andersson).

* We wish to understand the mechanism by which these PFPs permeabilize
  membranes.

* Often studied by measuring the kinetics of permeabilization of model
  membranes using PFPs.

* The problem: interpreting the leakage data in a way that informs us about the
  pore formation mechanism.

* There is current interest in using these kinetic studies to inform not just
  about single peptide agents, but also regulated multi-protein systems such as
  the Bcl-2 family of proteins.

* The theoretical challenge associated with analyzing leakage data is that the
  leakage measurement is not of a single large vesicle, but rather an ensemble
  of discrete vesicles, each operating as a separate reaction compartment. This
  is generally addressed by treating membrane as a single compartment and
  having reactions proceed according to the relative molar concentrations of
  the proteins there.

* Here we develop a novel simulation approach that explicitly enumerates
  hundreds of compartments and then uses that to simulate bulk observables
  such as fraction of protein bound, fraction of leakage, etc. Using this
  numerical simulation as a reference, we evaluate the previously described
  single-compartment models for their ability to correctly identify the
  underlying pore forming mechanism.

* We find that while the models previously described for peptides work for the
  mechanisms they were meant to address, extending them to new mechanisms of
  interest leads to problems that must be addressed by the multi-compartment
  approach. Such mechanisms include positive or negative regulation by other
  proteins, or positive feedback in pore formation due to auto-activation or
  aggregation. However, the one-compartment approximation can be effective even
  in these cases when certain conditions are met.

* We also find that calculation of Hill coefficients, etc. are not accurate
  determinants of stoichiometries of pore complexes, as previously suggested.

* Given that the Bcl-2 family protein network exhibits all of these features,
  we conclude that the multi-compartment simulation is a necessary tool when
  studying the mechanisms of pore-forming proteins with complex regulation.

Modeling Apoptosis and MOMP
---------------------------

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

