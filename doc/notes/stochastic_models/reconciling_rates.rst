Reconciling stochastic and deterministic models of dye release
==============================================================


In a continous/deterministic model of dye release (e.g., as implemented in the
module :py:mod:`tbidbaxlipo.models.one_cpt`), all reactions occurring on
vesicles are considered to occur in an undifferentiated compartment, and one
speaks of the bulk concentration of membrane-bound versus aqueous Bax, for
example. In a discrete/stochastic model (e.g., as implemented in
:py:mod:`tbidbaxlipo.models.n_cpt` and :py:mod:`tbidbaxlipo.models.site_cpt`,
the vesicles are explicitly represented as different compartments, each with a
certain (discrete) number of protein molecules associated with it. To explore
conditions under which the deterministic approximation fails to represent the
underlying stochastic reality, it is useful to make parallel models using both
approaches. To determine if the behavior of the stochastic or deterministic
models is truly similar or different, it is necessary to carefully consider the
treatment of rate constants and other quantitative parameters of both types of
models.

Converting bimolecular rate constants
-------------------------------------

The following is a derivation of the relationship between stochastic
and deterministic rate constants that will be familiar to many readers.

We begin by assuming a (solution) reaction volume of :math:`v`. For any two
(Molar) concentrations :math:`C_1` and :math:`C_2` of species :math:`S_1` and
:math:`S_2` in the deterministic system, consider a bimolecular forward
reaction governed by (the expressions in brackets provide the units of the
variables they accompany):

.. math::

    v_f^{det} \left[\frac{mole}{L \cdot s}\right] = k_f^{det} \left[ \frac{L}{mole \cdot sec} \right] \cdot C_1 \left[\frac{mole}{L} \right] \cdot C_2  \left[ \frac{mole}{L} \right]

where :math:`k_f^{det}` is the deterministic forward rate constant. That is,
the velocity of the reaction is determined by the product of the concentrations
of the two species and a rate constant.

The total number of molecules :math:`N_1` of :math:`S_1` in the system is given
by

.. math::

    N_1 = N_A \cdot C_1 \cdot v

where :math:`N_A` is Avogadro's number (which has units of molecules per mole);
similarly for :math:`S_2`. Since the velocity of a stochastic reaction is given
in molecules per second, the deterministic and the stochastic fluxes must be
related by

.. math::

    v_f^{det} \left[\frac{mole}{L \cdot s}\right] = v_f^{stoch} \left[\frac{molec}{s} \right] \cdot \frac{1}{N_A \cdot v} \left[\frac{mole}{molec \cdot L}\right]

Starting from the equation for the deterministic flux :math:`v_f^{det}`, above,
and substituting in for :math:`v_f^{det}, C_1,` and :math:`C_2`, we get 

.. math::

    v_f^{stoch} \left[\frac{molec}{s} \right] \cdot \frac{1}{N_A \cdot v} \left[\frac{mole}{molec \cdot L}\right] = k_f^{det} \cdot \frac{N_1}{N_A\cdot v} \cdot \frac{N_2}{N_A \cdot v}

The equation for the stochastic flux is therefore given by

.. math::

    v_f^{stoch} \left[\frac{molec}{s} \right] = \frac{k_f^{det}}{N_A \cdot v} \cdot  N_1 \cdot N_2

and the stochastic and deterministic rate `constants` are related by

.. math::
    k_f^{stoch} = \frac{k_f^{det}}{N_A \cdot v} \qquad \qquad k_f^{det} = k_f^{stoch} \cdot N_A \cdot v

Rescaling bimolecular stochastic rate constants
-----------------------------------------------

However, we do not wish to simulate all :math:`N_A \cdot C_1 \cdot v \approx
10^{23} \cdot 10^{-9} \cdot 10^{-5} \approx 10^9` molecules in the system in an
agent-based simulation. Hence we `rescale` our system by declaring the number
of molecules/agents that we will use to represent the molar concentrations
:math:`C_1` and :math:`C_2`; let us denote these as :math:`N^R_1` and
:math:`N^R_2`.

The ratio between the actual number of molecules in the system and the rescaled
number is

.. math::

    R = \frac{N^R_1}{N_A \cdot C_1 \cdot v} = \frac{N^R_1}{N_1}

We can therefore convert from the actual number to the rescaled number with
:math:`N_1 \cdot R = N_1^R`. In the rescaled system, our flux is not molecules
per second, but rescaled molecules per second, or :math:`molec^R`.

Writing out the equation for the rescaled bimolecular flux, starting from the
equation for the stochastic flux:

.. math::

    v_f^{R} \cdot \frac{1}{R} = k_f^{stoch} \cdot N_1^R \cdot \frac{1}{R} \cdot N_2^R \cdot \frac{1}{R} 

    v_f^{R} = \frac{k_f^{stoch}}{R} \cdot N_1^R \cdot N_2^R 

So the rescaled rate is equal to the true stochastic rate as

.. math::

    k_f^R = \frac{k_f^{stoch}}{R} \qquad \qquad k_f^{stoch} = k_f^R \cdot R

The deterministic rate can therefore be related directly to the rescaled
stochastic rate as 

.. math::

    k_f^R = \frac{k_f^{stoch}}{R} = \frac{k_f^{det}}{N_A \cdot v \cdot R}

Suppose that we express our value for :math:`N_1^R` (as a number of molecules) as
a multiple :math:`k_{rsf} \cdot C_1` of the numerical value of the value of
:math:`C_1`, regardless of its units (e.g. :math:`nM^{-1}\ s^{-1}`, :math:`\mu
M^{-1}\ s^{-1}`), and we do so for all species in the system (for example, if
:math:`C_1` is 1 nM, we let :math:`N_1^R` be 1000 molecules in the rescaled
system, giving the rate scaling factor :math:`k_{rsf}` a value of 1000). Then
the ratio :math:`R` becomes

.. math::

    R = \frac{k_{rsf} \cdot C_1}{N_A \cdot C_1 \cdot v} = \frac{k_{rsf}}{N_A \cdot v}

and the deterministic and rescaled stochastic rates are related by

.. math::

    k_f^R = \frac{k_f^{det}}{k_{rsf}}

When the rescaled value chosen is equal to the concentration, then
:math:`k_{rsf} = 1` and the rescaled stochastic rates and the deterministic
rates are also equal.

Rescaling translocation rates
-----------------------------

Now consider the basic translocation reaction of a protein :math:`P` to
membranes:

.. math::

    P_c \rightarrow P_m

where the subscripts :math:`c` and :math:`m` denote cytosol and membrane,
respectively. In a deterministic model, the rate of this reaction is dependent
on the bulk concentration of vesicles :math:`[Ves]`, that is

.. math::

    [P_c] + [Ves] \overset{k_f}{\rightarrow} [P_m] + [Ves]

with the ODE for :math:`P_c`:

.. math::

    \frac{dP_c}{dt} = - k_f [P_c][Ves]

In models where vesicles have unlimited capacity for binding protein, the
concentration :math:`[Ves]` is a constant :math:`Ves_0`, and the reaction
is pseudo-first-order.

For discrete, compartmentalized models, the situation is seemingly more
complicated. Here there are many reactions

.. math::

    P_c \rightarrow P_{m1}


    P_c \rightarrow P_{m2}

    \ldots

    P_c \rightarrow P_{mn}

where :math:`n` is the number of vesicles. However, the solution is fairly
simple; we create a list of :math:`n` reactions each with (unimolecular)
forward rate :math:`k_f`---with the same value as the constant used for the
deterministic model. This yields an overall flux of

.. math::

    \frac{dP_c}{dt} = -k_f [P_c] n

If we use a value of :math:`n` (i.e., use a number of discrete compartments)
that is numerically identical to the concentration :math:`[L]`, the overall
flux is equivalent between the discrete and continuous models.

The forward translocation rates used in the stochastic model implementations
are thus standardized as being the stochastic rates for a protein to
translocate to an `individual` compartment (rather than to `any` compartment).
Every compartment gets its own translocation rule ``p@sol -> p@cpt.`` Hence
when these are added together, they generate an aggregated forward flux that is
equal to the sum over all the individual translocations.

Moreover, when the size of the stochastic system is scaled up, the
translocation rate does `not` need to have the rate scaling factor applied,
since the rate is not actually bimolecular (it isn't dependent on the amount of
vesicles except insofar as the amount of vesicles dictates the number of
compartments and reactions). In the scaled-up system, each individual protein
retains the same propensity to find a vesicle; the larger number of proteins is
simply divided across a larger number of vesicles, and the overall relative
flux should be the same **(a little hand wavy)**.


In the deterministic case, the single rule ``p@sol->p@mem`` needs to have a
forward rate multiplied by the amount of vesicles (in Molar).

To elaborate on this, suppose :math:`k_f` is the stochastic rate constant for
translocation of a protein to a `particular` compartment. That is, it is
the rate constant for every rule (in the site-based model) of the type::

    Rule('p_to_cpt_i', p(cpt='solution') >> p(cpt='c_i'), k_f)

with :math:`n` rules of this type for :math:`n` compartments.

Hence the aggregate flux of p from solution to any compartment is

.. math::

    \frac{dp_{sol}}{dt} = \sum_n -k_f \cdot p_{sol} = -k_f \cdot p_{sol} \cdot Ves_0

where :math:`Ves_0` is the (fixed) total number of compartments in the system.

Now suppose we define two new variables, :math:`p^R_{sol}` and :math:`Ves^R_0`,
indicating the rescaled variable where

.. math::

    p^R_{sol} =  k_{rsf} \cdot p_{sol}

    Ves^R_0 = k_{rsf} \cdot Ves_0

The point is not that the two systems should have equal rates of
change--because they should not.  (It would appear that we need to divide the
rescaled forward rate constant by :math:`k_{rsf}`, but I have to follow this
through.

.. todo:: Discrepancy between site and compartment based implementations

    Oddly, the compartment-based approach appears to require the rate scaling
    factor, whereas the site-based approach does not?!!!  TODO something is
    amiss here.

Reactions at Membranes
----------------------

Suppose we assume that there is a fundamental forward rate that defines the
reaction propensity between tBid and Bax on a 100nm diameter liposome,
:math:`k_{s}`, where the :math:`s` stands for "stochastic." In the
discrete/compartment-based model, where we are dealing with numbers of
molecules per compartment, we then have that the forward rate, per compartment,
in the stochastic case is

.. math::

    v_{s} = k_s \cdot tBid \cdot Bax

With :math:`n` compartments, the aggregate flux, in molecules is

.. math::

    v_{s} = \sum_n k_s \cdot tBid_i \cdot Bax_i

with the subscripts for tBid and Bax denoting the number of molecules of tBid
and Bax on compartment :math:`i`. The units for :math:`v_s` resulting from this
expression are :math:`s^{-1}`, as is appropriate. However, the important thing
to note here is that as the number of compartments grows, the number of
molecules per compartment decreases. However, in the stochastic model this does
not need to be kept track of explicitly, since this will automatically be
reflected in the concentrations of tBid and Bax in each individual compartment.
Moreover, the total flux over the :math:`n` compartments also does not need to
be kept track of explicitly, since 

Turning to the continuum case, we wish to write the expression for the
deterministic flux in terms of the original stochastic forward rate,
:math:`k_s`. The resulting units should be in :math:`molar^{-1} s^{-1}`. Here
we don't know how much tBid or Bax is in each compartment, only how much is
associated with membranes in total (denoted with the subscript :math:`m`). We
imagine dividing the pool of membrane-bound tBid and Bax among the :math:`n`
vesicles equally, so the rate in each individual compartment is

.. math::

    v_s = k_s \cdot \frac{[tBid_m]}{[Ves]} \cdot \frac{[Bax_m]}{[Ves]}

There are two problems with this, first, that it only gives us the average
flux, ... Crucially here though, the `total` rate of tBid/Bax binding is 

