Outline
=======

* Methods of modeling liposome permeabilization kinetics
    - Review of previously described methods
    - Comparison of compartmentalized vs. continuous representations of
      liposomes

* Molecular mechanism of Bax membrane insertion during activation and pore
   formation
    - What is the conformational pathway of Bax as it inserts into the
      membrane?  Does it depend on the mechanism of activation?
    - Identifying the number of conformational intermediates by model comparison
    - Bax inserts into membranes before pore formation
    - Dynamics of Bid-Bax rearrangements during activation

* Kinetics of liposome permeabilization by Bid and Bax
    - What determines kinetics of pore formation by Bax.  In particular, what
      are the relative roles of membrane binding, activation by BH3-only
      proteins, and auto-activation by Bax itself
    - Role of liposome binding sites for Bid but not for Bax
    - Saturability of Bid by Bax
    - Role of auto-activation at different activator strengths
    - Measurement of Bax distribution by single-molecule microscopy

* Stability and stoichiometry of the Bax pore
    - What is the nature of the pore? Lipidic vs. proteinaceous, stable vs.
      transient, fixed vs. variable.
    - Determinants of all-or-none vs. graded release
    - Minimal stoichiometry of the pore, inference from bulk measurements
    - Minimal stoichiometry of the pore, single-molecule measurements

Chp 1. Modeling pore formation
------------------------------

"Compartmentalization of reactions at membranes plays a critical role in the
rate and extent of apoptotic pore formation by Bax"

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

Chp 2. Bax insertion pathway
----------------------------

Chp 3. Kinetics of liposome permeabilization by Bid/Bax
-------------------------------------------------------

- Introduction

  - Cells need to be able to reliably execute apoptosis, not dying when they
    don't mean to, but also able to execute apoptosis when needed.

  - The measured concentrations of Bcl-2 family proteins is highly variable.
    How is the apoptosis network able to execute the decision over such a wide
    variation of concentrations?

  - Suggests concept of apoptotic phase space, in which Bid (activator) and
    Bax (effector) represent the protein axes; lipid axis is also important.
    Anti-apoptotic axis is also clearly important, but important to first
    define in the absence of anti-apoptotics.

  - So goal: measure a two-d dose response for Bid and Bax.

  - The Bax dose response has been measured in a variety of fashions, with
    differing conclusions: cooperative, saturating, and linear.

  - In addition, the response surface tells us about the underlying biochemical
    mechanism. Can we use models to explain the underlying response surface
    in terms of the mechanism?

- Results

  - **Pore formation kinetics, and Bax insertion, plateau at high Bax/liposome ratios when Bid is the activator.**

    - Satsoura et al. previously described saturation of liposomes by Bax at
      numbers as low as 20.

    - Bax titration at single liposome concentration, Bid titration, Bim BH3
      titration, ANTS

    - Saturation of initial rate? Saturation of kinetics? difference between
      cBid and Bid BH3?

    - Choose a model and show scaling of parameters (e.g., Kushnareva model)

    - Permeabilization kinetics saturate for cBid but not for Bim BH3

    - H1: due to stoichiometric limitation of Bax insertion sites

    - H2: due to saturation of the Bax:liposome complex (Bax
      encounter sites)

    - H3: due to saturation of the activator cBid.

  - **Liposome concentration affects the rate of Bax insertion, but is not
    stoichiometrically limiting**

    - Liposome titration experiment shows that the amount of liposomes affects
      the efficiency of Bax recruitment, as predicted by simple model; it
      affects the forward rate of Bax insertion, and maybe also the equilibrium
      amount, harder to say, but this could be due to the Bax off-rate.

    - Bax competition experiment shows that H1 and H2 are not true--addition of
      Bax in a Bim BH3 activation scenario does not limit the rate or extent
      of Bax insertion.

  - **Bax saturation is due to the saturation of the activator Bid.**

    - *In Bax insertion assay with Bax titration, adding more Bid increases
      critical concentration of Bax inhibition*

  - **Effect of liposome concentration on Bid and Bax binding**

    - Bid binding, FRET?

    - Bax binding, BH3 peptide

  - **Bid binding to membranes saturates at low stoichiometries**

    - Bid-membrane FRET experiments show competitive binding, however
      the results are not well-fit by single-site competitive binding!
      Linear, rather than hyperbolic competition curve

    - *Bid-membrane gel filtration*

    - *Bid-647 membrane FCS to estimate fraction bound*

    - *Bid488-mito binding by fluorescence*

  - **The extent of liposomes permeabilized is determined exclusively by
    the Bax/liposome ratio, not by the amount of activator.**

    - *Result from Bid titration*

    - *Result from Bax-NBD insertion curves that show that there is no late
      linear phase, as would be expected with product inhibition.*

    - This rules out the possiblity that Bid inhibits Bax at high concentration
      by occupying its BH3 groove.

  - **Incomplete permeabilization is due to ???** 

    - Due to exponential nature of pore formation (i.e., read only first event)

    - Due to progressive, irreversible depletion of the Bax and hence a
      slowdown in rate

    - Due to transient binding of the activated Bax to the soluble Bax that
      effectively competes it away.

  - **Role of auto-activation**

    - From model: relevance of auto-activation is dependent on the baseline
      insertion activity: if baseline probability of insertion is low, then
      auto-activation plays a relatively larger role, but it's role is local.
      On the other hand, if baseline activation is higher, it plays a minimal
      role.

    - and in Bax-limited regimes, auto-activation can actually decrease the
      overall extent of permeabilization(?)

    - Indeed, if there is a Bax-hole, the auto-activation activity of Bax may
      actually limit the overall capability.

    - To test the Bax-hole: do a liposome assay, come to steady-state, then
      add more Bax; if the fractional permeabilization you get is less then
      what you would have gotten at the start, then you have a Bax hole.
      Do this with varying levels of starter Bax.

Chp 4. Stability and stoichiometry of the pore
----------------------------------------------

  - **All-or-none pore formation and stable insertion by Bax is dependent on
    the presence of activator.**

    - *Activation of Bax by heat does not produce substantial amounts of
      inserted Bax.*

    - *Activation of Bax by heat leads to graded release(?)* Tried this
      experiment before, need to repeat due to complete failure.

    - *Activation of Bax by heat leads to*

  - **Stoichiometry of the pore**

    - *Can it be estimated from the amount of permeabilization achieved
      at low Bax/lipo ratios* (when auto-activation would be presumed to be
      minimal?). Idea would be to operate in a regime where the probability of
      being inserted would be high, and relatively independent of the amount of
      Bax already on the liposome (e.g., a high BH3 or Bid situation). In
      this regime the effect of the exponential nature of permeabilization
      would also have a minimal effect, because the probability of getting
      "more than one pore" per liposome would be small.


Miscellaneous
-------------

Motivation: to understand mechanistic basis of Bax pore formation kinetics.

- *** Identifiability of the parameters?***

2. Screening of model ensemble against in vitro kinetic data.

- Fit to primary data, not reduced (feature-based) data.

- All models fit poorly, no mechanistic additions increase fit

- Hypothesis from shape of kinetic curve: saturable Bax-liposome association complex

- Adding this feature allows all models to fit regardless of protein interaction topology

3. Experimental evidence for saturable association complex

- Localization kinetics, measured by NBD-fluorescence or by membrane FRET, slow down at higher Bax (***but approach same steady state level?***)

- Measuring the association complex by FRET: fraction of Bax bound saturates at high Bax/liposome ratios

- Measure by FCS

- ***To be reconciled***:
    - gel filtration shows saturation of equilibrium binding at high Bax/lipo ratios
    - NBD fluorescence appears to show that same steady state level is approached

4. Mechanistic basis for saturability

- ***Repeat above experiments in presence/absence of cardiolipin (negatively charged lipids)***


5. Accounting for non-independent binding across liposomes

- It is generally believed that Bax binding is cooperative

- Cooperatvity implies non-independent binding to liposomes

- Multi-compartment model shows that non-independent binding would change bulk kinetics from what is observed

- Suggests that cooperativity/Bax aggregation may play a less dominant kinetic role

- ***Does K21E mutation inhibit membrane binding or only oligomerization?***

- Need for single vesicle approaches to measure Bax/liposome distribution.
