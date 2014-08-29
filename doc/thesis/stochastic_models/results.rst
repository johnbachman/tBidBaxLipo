Results
=======

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

Because the MC model tracks the number of pores formed on each individual
compartment, we can also calculate the distribution of pores across the
population of vesicles. As described by Schwarz, the number of pores per
vesicle follows a Poisson distribution with the mean obtained from the
deterministic simulation of the 1C model [Schwarz1990]_ (:ref:`Figure 1B
<stochastic_models_fig1>`).

Using the 
Using the 
Using the evaluation

**Figure 1B** For the purposes of this paper, we focus on comparisons with
previous approaches in which the dye release process is modeled as resulting
from the formation of stable pores, rather than transient disruptions in the membrane leading to graded release.

As a validation that the multi-compartment simulation algorithm has been
correctly implemented, we instantiate a very simple model of pore formation and
compare it to several of the previously described pore formation models. As
expected, when rate parameters are set to corresponding values (see Methods),
the MC model duplicates the dye release kinetic curves produced by two of the
three previously described models (Schwarz, Almeida). Moreover, the
distribution of pores across vesicles in the MC simulation matches a Poisson
distribution with a mean pores per liposome value calculated according the
Schwarz.

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

* **Adding a second protein breaks other methods when concentrations are
  low, unless...**

* Adding an activator protein, such as Bid,

* **Adding auto-activation breaks other methods, unless...**

* Bax is believed to auto-activate.

* **Hill coefficient analysis is not a reliable indicator of stoichiometry**

* Perturbation theory explanation?

* **In fitting permeabilization curves with exponentials, it is essential to
  account for Fmax as well as k**


B + L <> BL >> BL*

What I am trying to explain:

    - non-origin nature of slope of Bax permeabilization?

    - Show that reaction topology determines whether the continuum model
      matches the compartment model.

Need to show experimentally true as well as theoretically true

Coins/buckets argument

    - hinges in part on the fact that the curve is a two-parameter curve, with
      both k and fmax.

    - Both enzyme and pore formation case don't provide explanations for why
      fmax is less than 100%.

**Evaluation of permeabilization models for individual perm. curves**

    - This could potentially go in the liposome perm kinetics chapter.

Figure: Example permeabilization curve.

Table listing models, with references and features

    - Exponential model (1 and 2 and 3 sum exponentials)

    - Schwarz: log transform the data to estimate "number of pores"

    - One and two exponential equations (history of this equation from Almeida,
      Schwarz, Schlesinger

    - Kushnareva/Newmeyer model: enzymatic style

    - European group paper?

Bax specific:

    - Phenomenology: a delay; nearly exponential activity; maximal activity
      below 100% permeabilization; slow rise;

    - At start, you have no pores nucleated, auto-activation helps
      get pores nucleated, hence the acceleration. However, this
      starts to fight against the depletion of Bax due to recruitment
      to existing pores, and eventually depletion wins out.

    - Three velocities: initial, intermediate, final; pore production is
      linear at each one? dp/dt = k

    - Two-phase scaling of the kinetic constant, k

    - Hyperbolic scaling of the Fmax

**Prediction of role of auto-activation**

    - Auto-activation may deplete 

**Refute notion that linearity in slope indicates non-saturation and
non-cooperativity!**

    - Show timescale separation analysis??


