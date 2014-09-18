Discussion
==========

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


