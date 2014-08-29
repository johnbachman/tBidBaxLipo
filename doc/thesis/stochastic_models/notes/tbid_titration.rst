.. _tbid_titration:

Dye release under tBid titration
================================

As shown in the figures for :ref:`tbid_bax_activation`, the deterministic model
approximates the stochastic model at the level of the bulk observables,
but when tBid binds liposomes very tightly, the distribution of tBids among
the liposomes becomes uneven and very non-Poissonian. While this may not
be manifest at the level of bulk observables such as inserted Bax,
it should have an effect on the experimentally observed dye release kinetics.

To examine this I ran simulations with tBid concentrations from 0 to 20 nM,
while keeping liposomes at 5 nM and Bax at 100 nM (and keeping tBid
binding as irreversible):

.. ipython::

    In [1]: from tbidbaxlipo.plots.stoch_det_comparison.tbid_titration \
       ...: import jobs, data

    In [2]: from tbidbaxlipo.plots.stoch_det_comparison.plots import *

    In [3]: [job.params_dict for job in jobs]

Here we plot the dye release timecourses alongside the deterministic equivalent
that is expected via the Poisson assumption of Schwarz [Schwarz1990]_.  Here
the simulation predicts the numbers of pores per liposome(that is, the number
of pores is expected to follow the dye release according to the function
:math:`-log\ E(t) = P(t)`, where `E` is the efflux function (number of
liposomes with no pores), and `P` is the average number of pores per vesicle:

.. ipython::

    In [4]: plot_dye_release_titration(jobs, data);

    @suppress
    In [5]: savefig('_static/tbid_titration_1.png')

.. image:: ../../../_static/tbid_titration_1.png
    :width: 6in

Todos
-----

* Do a liposome experiment using concentrations where you might normally get
  about 100% dye release, but in this case pre-incubate the liposomes with
  tBid. Then add the liposomes 1:2 or so with liposomes that haven't been
  pre-incubated with tBid, and add Bax. Compare to the case where the same
  amount of tBid was incubated with the full amount of liposomes, getting
  more evenly distributed across them. Expectation would be that the Fmax
  for dye release would be somewhere near half for the pre-incubated case.
  On the other hand, one might find that the bulk observables (like NBD-Bax-126
  insertion) would be the same in both cases.
* Do pre-incubation experiment, and compare to Bim BH3 treatment (after
  measuring Bim BH3 off-rate from liposomes using the Octet).
* Do immobilized liposome experiment, looking at distribution of tBids,
  and turnover rate after photobleaching. Follow-up would be to add Bax
  and see if only liposomes containing Bax are immobilized.
