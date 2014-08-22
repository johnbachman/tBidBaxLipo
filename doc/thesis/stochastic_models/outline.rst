Outline
=======

"Compartmentalization of reactions at membranes plays a critical role in the
rate and extent of apoptotic pore formation by Bax"

Intro
-----

* Permeabilization of membranes by peptides and proteins is important for
  bacterial toxins, pore forming peptides used as therapeutics, and in the
  processes of apoptosis and necroptosis.

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

Results
-------

* **The simulation approach, validation**

* Figure 1A: the simulation approach works by enumerating individual reaction
  compartments. For a given pore formation mechanism, we can enumerate a
  multi-compartment version or a single-compartment version. This requires a
  careful accounting for rates (see Methods).

* Figure 1B: the multi-compartment approach validates the Poisson
  assumption/pores approach described by Schwarz for unimolecular pore
  formation mechanisms without complex regulation.

* Figure 1C: the multi-compartment approach shows that the enzymatic approach
  is wrong?

* Figure 1D: the multi-compartment approach shows that Almeida et al., is
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


