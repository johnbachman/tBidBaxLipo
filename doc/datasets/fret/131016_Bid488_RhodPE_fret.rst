Bid-Alexa-488 and Rhod-PE liposome FRET (10/16/13)
==================================================

In this experiment I wanted to see if I would be able to measure FRET between an
Alexa-488 labeled protein (Bid in this case) and Rhodamine-labeled liposomes.
I kept cBid-488 concentrations constant at 50 nM and did a liposome titration
from 0.015 nM to 15.5 nM (1 mg/mL).

I combined the Bid and liposomes by multi-channel pipette and then started
reading the plate. As controls, I had Bid with 0 nM liposomes (donor-only
control) as well as one complete dilution series of the liposomes (acceptor
background control).

Though I don't intend to follow up much on this experiment, since Bid dynamics
is not my main priority, there are some interesting observations that can be
made from this data:

* Bid/liposome FRET has a rapidly equilibrating component as well as a slower
  additional component. This could be due to a conformational change that
  results in increased FRET, either due to a change in the intrinsic FRET
  efficiency of the conformation, or due to a shift in the population to a pool
  of molecules that have a slower off-rate, leading to accumulation at the
  membrane.
* The titration suggests that both the fast pre-equilibrium and slower
  equilibrium FRET may have a conventional Hill-type binding component as well
  as a component that goes up linearly with the liposome concentration. This
  may be "non-specific binding" of the liposomes similarly to what I've
  hypothesized for Bax, that is, some kind of interaction with non-preferred
  lipid components that is nonetheless productive. Alternatively, it could
  simply represent non-productive collisions with liposomes in solution (that
  don't result in further interactions with the membrane).

First we look at the raw, averaged, and background-subtracted data. For the
background subtraction, the fluorescence of the liposomes in the 488 channel
(due to bleedthrough) was subtracted from the experimental wells with the
corresponding liposome concentration.

.. plot::
    :context:

    from tbidbaxlipo.plots.layout_131016 import *
    plot_data()

Next, we fit the Bid-only condition to a line so that we can divide by
smoothed values (so that we don't compound the noise of the experimental
conditions by dividing by something noisy):

.. plot::
    :context:

    close('all')
    bid_only = plot_fit_donor_only()

With the fitted Bid-only values we can now calculate FRET for the
background-subtracted, averaged experimental timecourses. In addition, we
can fit these FRET timecourses to a single exponential function of the form

.. math::

    F_{min} + F_0 e^{-k t}

Interestingly, while all the concentrations of liposomes have already reached
a rapid pre-equilibrium, there is a steady increase in FRET apparent at the
highest liposome concentrations:

.. plot::
    :context:

    close('all')
    fret = plot_fret_exp_fits(bid_only)

We can now plot the value of the FRET efficiency, in this case at the first
timepoint. In addition we can attempting to fit the FRET titration curve with a
variety of functions. The quadratic fit represents a fit to the quadratic
solution of the binding equations, with the concentration of the liposomes
represented in terms of lipids (multiplied by 83775). It turns out to be pretty
much identical to the conventional Hill equation in fit.

.. plot::
    :context:

    close('all')
    plot_fret_titration(fret, 0, 1)

While the exponential curve appears to clearly be a worse fit, the jog in the
curve at the upper two concentrations makes it difficult to argue strongly for
or against the conventional Hill binding isotherm vs. the non-specific binding
curve. The non-specific binding curve is basically a simple binding isotherm
with an additional term that is proportional to the liposome concentration.  It
may simply reflect an increase in FRET that is due to random collisions, since
it goes up linearly in the quencher/acceptor concentration.

.. math::

    F_{max} \frac{x}{K_D + x} + F_{nsb} x

For comparison, we repeat the plot for the FRET efficiency averaged over the
last three timepoints and find that the results are similar:

.. plot::
    :context:

    close('all')
    plot_fret_titration(fret, -3, None)

