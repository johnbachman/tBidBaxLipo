Results
=======

**The simulation approach, validation**

For a proposed pore formation mechanism we can produce a two-compartment model
or a corresponding multi-compartment model and determine whether the
two-compartment model is equivalent to the multi-compartment model, in the
limit of many vesicles.

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


MC model duplicates the dye release kinetics
recapitulates previously
validated findings from the the findings previously The multi-compartment
approach validates the Poisson assumption/pores approach described by Schwarz
for unimolecular pore formation mechanisms without complex regulation.

.. plot::

    from tbidbaxlipo.plots.stoch_det_comparison.translocation import \
        plots, jobs, data
    plots.plot_timecourse_comparison(jobs, data, 0)

.. image:: ../../_static/simple_translocation_1.png
    :width: 6in

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


