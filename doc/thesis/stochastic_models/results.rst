Results
=======

Comparison with previous simulation approaches
----------------------------------------------

To study the mechanism of multi-protein pore formation processes we developed a
simulation approach that involved the enumeration of many vesicle compartments
(:ref:`Figure 1A <stochastic_models_fig1>`). Previously described mathematical
modeling approaches to studying pore formation treated the model membranes as
an undifferentiated lipid compartment in which the PFPs would operate. In
simple scenarios, the one-compartment (1C) modeling approaches are advantageous
since they are amenable to analytical rather than strictly numerical analysis,
which makes parameter estimation and other things faster and simpler. Thus we
first sought to determine cases in which the 1C approach is sufficient for
describing pore formation mechanisms by comparing results of simulations from
1C and MC models and determining whether the 1C model was an adequate
representation of the underlying multi-compartment system.

To do this, we developed a modeling approach that allows us to produce
corresponding 1C and MC models for any proposed pore formation mechanism.
Maintaining equivalence requires appropriate conversions of the rate parameters
involved (:ref:`Methods <reconciling_rates>`). The 1C model is formulated
in terms of ordinary differential equations, which are integrated numerically,
while the MC model is simulated using the stochastic simulation algorithms
implemented in BioNetGen or Kappa (:ref:`Methods <simulation_methods>`). We
can then analyze the simulation results to determine whether a given 1C model
is equivalent to its corresponding MC model in the limit of many vesicles.

To validate the correctness of our implementation, we formulated a simple
model of vesicle permeabilization, consisting of the following
elementary reactions [Schwarz1992b]_:

.. math::

    P_c + L \rightleftharpoons P_l

.. math::

    m \cdot P_l \rightarrow Pore

Where :math:`P_c` and :math:`P_l` represent the peptide or protein in its
aqueous (cytosolic) or lipid-associated state, :math:`L` represents the
concentration of vesicles, and :math:`m` is the number of monomers involved
in pore formation.

As expected, when rate parameters are set to corresponding values (see
:ref:`Methods <reconciling_rates>`), the MC model duplicates the dye release
kinetic curves produced by two of the three previously described models
(Schwarz, Almeida).

The 1C model shown in (1) above keeps track of the number of pores formed over
time; the fraction of vesicles permeabilized can be calculated by assuming a
Poisson distribution of pores across vesicles [Schwarz]_. In the MC model, we
explicitly keep track of both the number of pores formed on each vesicle as
well as whether it has been permeabilized or not (has > 1 pores). As described
by Schwarz, the number of pores per vesicle follows a Poisson distribution with
the mean obtained from the deterministic simulation of the 1C model
[Schwarz1990]_ (:ref:`Figure 1B <stochastic_models_fig1>`).

In addition to the Poisson/pore approach of Schwarz, several other approaches
that explicitly track the fraction of vesicles permeabilized have also been
described. These include treating the permeabilization of vesicles as an
enzymatic conversion mediated by the pore forming protein [Kushnareva2012]_;
tracking the fraction of permeabilized vesicles in a set of ODEs
[Almeida2009]_; or empirical approaches describing the permeabilization curve
as a two- or three-parameter exponential equation. 

**The problem with Kushnareva** is (1) that it doesn't account for the
saturation of PFPs that occurs at low P/L ratios. In this way it suffers from
the same problems as a pseudo-first order approximation of catalysis. In
addition, the model doesn't account for the rebinding of empty vesicles.

Almeida accounts for rebinding of empty vesicles, **but the problem with
Almeida** is that it doesn't account for stable binding of protein/peptide to
membrane.  In addition, doesn't account for the population of PFPs that are in
the pore state on empty vesicles. If the pore-state is long-lived, this will
affect the free concentration of P and hence the kinetics.

**Both these models have the problem** that they don't explain cases where the
permeabilization curve plateaus at levels less than 100%.

**Figure 1B** For the purposes of this paper, we focus on comparisons with
previous approaches in which the dye release process is modeled as resulting
from the formation of stable pores, rather than transient disruptions in the membrane leading to graded release.

Notably, the pseudo-first order enzymatic approach described by Kushnareva et
al. does not match the results from the multi_cpt simulation. This is due to
the fact that this approximation does not account for depletion of the pool of
P due to binding to either permeabilized or permeabilized vesicles. As a result
it is an effective approximation only when Pbound << P_0.

The multi-compartment approach validates the Poisson assumption/pores approach
described by Schwarz for unimolecular pore formation mechanisms without complex
regulation.

**Figure 1C**. The multi-compartment approach shows that the enzymatic approach
is wrong?

**Figure 1D**. The multi-compartment approach shows that Almeida et al., is
  right?

Table listing models, with references and features

    - Exponential model (1 and 2 and 3 sum exponentials)

    - Schwarz: log transform the data to estimate "number of pores"

    - One and two exponential equations (history of this equation from Almeida,
      Schwarz, Schlesinger

    - Kushnareva/Newmeyer model: enzymatic style

    - European group paper?

Discrepancies in predictions of permeabilization kinetics
---------------------------------------------------------

* Show that reaction topology determines whether the continuum model
  matches the compartment model.

* Adding a second protein breaks other methods when concentrations are
  low, unless...

Discrepancies in predictions of binding
---------------------------------------

* **Adding auto-activation breaks other methods, unless...**

* Bax is believed to auto-activate.


Extraction of rates by fitting exponentials yields incorrect results
--------------------------------------------------------------------

* Non-origin nature of slope of Bax permeabilization.

* **In fitting permeabilization curves with exponentials, it is essential to
  account for Fmax as well as k**

* Coins/buckets argument

    * hinges in part on the fact that the curve is a two-parameter curve, with
      both k and fmax.

    * Both enzyme and pore formation case don't provide explanations for why
      fmax is less than 100%.

**Refute notion that linearity in slope indicates non-saturation and
non-cooperativity!**

    - Show timescale separation analysis??

Inferring stoichiometry
-----------------------

* **Hill coefficient analysis is not a reliable indicator of stoichiometry**

* Perturbation theory explanation?


