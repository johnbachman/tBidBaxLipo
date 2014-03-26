.. _140318_c126_bax_lipo_titration:

NBD-c126-Bax with liposome titration and Bim BH3 (3/18/14)
==========================================================

The purpose of this experiment was to look at how the Bax/liposome ratio
affects the rate of Bax insertion while removing the possible effects of
interaction between Bid and Bax (e.g., saturation of Bid by Bax) and the
effects of interaction between Bid and liposomes. The data from 9/11/13
suggested a saturation phenomenon in which increasing Bax slowed down the rate
of insertion of Bax, but it was not clear from that experiment whether this was
due to the saturation of Bid by Bax, or the saturation of liposomes by Bax.

Therefore in this experiment I used a high concentration (50 uM) of Bim BH3
as the activator, held Bax constant, and did a dilution series of liposomes
from 1 mg/mL lipid down to 0.001 mg/mL.

Ran plate at 37C, pre-incubated at this temp for 30 minutes before starting the
read. The liposomes and Bax were mixed approximately 90 seconds before the moment
of the first read.

As this was a pilot experiment, I only did one replicate.

Conclusions
-----------

* Lower liposome concentrations definitely makes Bax insertion slower.

* This is to be expected since insertion follows from mBax, which even
  according to the usual model of binding, should depend in a saturating
  fashion on the amount of liposomes (L/(L+KD)).

* This duplicates the finding of Justin in his two-D titration, in which
  reduction of liposomes makes Bax insertion slower. However, it indicates that
  failure of Bid to bind liposomes is not the reason for the slowdown--rather,
  it is the failure of Bax to bind the liposomes.

* It is not completely clear if the same steady-state amount of Bax insertion
  would be reached for each concentration of liposomes. This seems to
  be the case for the first few concentrations, however after that it is probably
  impossible to say (quantify "impossible" by MCMC?, i.e., the failure to reduce
  uncertainty about the value of the parameter). If Bax jumping is substantial, one
  might expect that the steady state amount of Bax inserted might be less for the
  lower concentrations due to the predominance of the jumping-off rate.

* One might argue that the slowdown is due to failure of the BH3 peptide to
  "bind" to the liposomes and thus find and activate Bax at the lower liposome
  concentrations. This objection makes sense only if we consider, like Bid,
  that the BH3 must first bind to the membrane before it can interact with Bax.
  How do we imagine that the interaction with BH3s takes place?


Raw data
--------

First, the raw data of all timecourses:

.. ipython::

    In [1]: from tbidbaxlipo.plots.layout_140318 import *

    In [2]: figure()

    In [3]: plot_all(timecourse_wells)

    In [4]: title("Raw timecourses")

    @suppress
    In [5]: plt.savefig('_static/140318_c126_bax_lipo_titration_1.png')

.. image:: ../../_static/140318_c126_bax_lipo_titration_1.png
    :width: 6 in

Background-subtracted timecourses
---------------------------------

The signal from each timecourse consists of fluorescence due to 1) the buffer
and plate, 2) the liposomes, 3) the Bim BH3 and DMSO, and 4) the Bax. We didn't
measure the effects of the Bim BH3 so we cannot specifically control for it;
however we can control for the other effects.

In addition, there is evidently photobleaching during the course of the
experiment, and we expect that the different background components will
photobleach at different rates.  Because the timecourse of the liposome-only
wells incorporates the effects of both the liposome fluorescence and the
expected photobleaching trajectory of the liposomes over time, we can simply
subtract these wells from the Bax-liposome wells to control for these effects.

.. ipython::

    In [2]: figure()

    In [3]: plot_all(bgsub_wells)

    In [4]: title("Background-subtracted Bax/lipo timecourses")

    @suppress
    In [5]: plt.savefig('_static/140318_c126_bax_lipo_titration_2.png')

.. image:: ../../_static/140318_c126_bax_lipo_titration_2.png
    :width: 6 in

