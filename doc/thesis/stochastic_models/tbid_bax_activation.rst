.. _tbid_bax_activation:

Activation of Bax by tBid
=========================

Now we look at a bimolecular system in which tBid and Bax both translocate
to liposomes and Bax is activated by binding with tBid on membranes:

.. ipython:: python
    :suppress:

    from pylab import *
    from tbidbaxlipo.plots.stoch_det_comparison.plots \
        import plot_hist_vs_poisson
    from tbidbaxlipo.plots.stoch_det_comparison.tbid_bax_activation.plots import *

.. ipython:: python

    from tbidbaxlipo.plots.stoch_det_comparison.tbid_bax_activation import *
    job = jobs[2]
    one_cpt_model = job.one_cpt_builder().model
    one_cpt_model.rules

This model has the nominal parameters of 5 nM liposomes, 20 nM tBid, and 100 nM
Bax; it also has an irreversible translocation process for tBid:

.. ipython:: python

    one_cpt_model.parameters

When we simulate it we find good correspondence between the stochastic and
deterministic models in terms of the total levels of cytosolic, mitochondrial,
and inserted (activated) Bax:

.. ipython::

    In [9]: plot_timecourse_comparison(jobs, data, 2)

    @suppress
    In [10]: savefig('_static/tbid_bax_activation_1.png')

.. image:: ../../_static/tbid_bax_activation_1.png
    :width: 6in

However, when we look at the distribution of activated Bax across the
liposomes at the final timepoint, we see that it is not at all Poisson
(that is, insertion events are not independent):

.. ipython::

    In [3]: plot_hist_vs_poisson(jobs, data, 2, 'iBax', -1);

    @suppress
    In [4]: savefig('_static/tbid_bax_activation_2.png')

.. image:: ../../_static/tbid_bax_activation_2.png
    :width: 6in

However, the distribution of the tBid across liposomes is independent:

.. ipython::

    In [5]: plot_hist_vs_poisson(jobs, data, 2, 'mtBid', -1);

    @suppress
    In [6]: savefig('_static/tbid_bax_activation_3.png')

.. image:: ../../_static/tbid_bax_activation_3.png
    :width: 6in

This leads us to two (not necessarily mutually exclusive) hypotheses for the
primary reason why Bax is being differentially recruited to liposomes:

- The irreversibility of tBid binding to liposomes is what sets up the disparity
  between liposomes, since liposomes that never see tBid can never recruit Bax
- The concentration of tBid is the key feature; if tBid is in excess such that
  most liposomes have at least one tBid, then disparities will not result.

It also raises the questions:

- What happens when there is a liposome with large amounts of tBid that recruits
  a large amount of Bax but subsequently runs out of Bax binding sites?
  Does this put an upper bound on the distribution? How does this affect the
  bulk kinetics?

Moderate tBid concentration, reversible binding
-----------------------------------------------

We first look at the first possibility, that it is the artificial setting of
the reverse rate to zero that sets up the non-independence of Bax activation.
We run simulations in which we set the reverse of tBid from liposome binding
to 0.01 per second:

.. ipython::

    In [2]: job = jobs[3]

    In [4]: one_cpt_model = job.one_cpt_builder().model

    In [5]: one_cpt_model.parameters['tBid_transloc_kr']

    In [5]: one_cpt_model.parameters['tBid_0']

Correspondence of bulk observables to the deterministic
model is excellent for this system:

.. ipython::

    In [6]: plot_timecourse_comparison(jobs, data, 3);

    @suppress
    In [7]: savefig('_static/tbid_bax_activation_4.png')

.. image:: ../../_static/tbid_bax_activation_4.png
    :width: 6in

Here it seems that the distribution of activated Bax is much more Poissonian:

.. ipython::

    In [8]: plot_hist_vs_poisson(jobs, data, 3, 'iBax', -1);

    @suppress
    In [9]: savefig('_static/tbid_bax_activation_5.png')

.. image:: ../../_static/tbid_bax_activation_5.png
    :width: 6in

And as before, the distribution of tBids is independent:

.. ipython::

    In [5]: plot_hist_vs_poisson(jobs, data, 3, 'mtBid', -1);

    @suppress
    In [6]: savefig('_static/tbid_bax_activation_6.png')

.. image:: ../../_static/tbid_bax_activation_6.png
    :width: 6in

Low tBid concentration, irreversible binding
--------------------------------------------

To test the effect of tBid concentration on iBax distribution we start with
our initial assumption of irreversible binding and reduce the concentration of
tBid to 1 nM:

.. ipython::

    In [2]: job = jobs[0]

    In [4]: one_cpt_model = job.one_cpt_builder().model

    In [5]: one_cpt_model.parameters['tBid_transloc_kr']

    In [5]: one_cpt_model.parameters['tBid_0']

Interestingly, correspondence of bulk observables to the deterministic
model is excellent, even here:

.. ipython::

    In [6]: plot_timecourse_comparison(jobs, data, 0);

    @suppress
    In [7]: savefig('_static/tbid_bax_activation_7.png')

.. image:: ../../_static/tbid_bax_activation_7.png
    :width: 6in

But the distribution of activated Bax tells a different story--the distribution
is strongly bimodal, with over 80% of liposomes having no activated Bax, while
the remaining fraction are heavily loaded with Bax:

.. ipython::

    In [8]: plot_hist_vs_poisson(jobs, data, 0, 'iBax', -1);

    @suppress
    In [9]: savefig('_static/tbid_bax_activation_8.png')

.. image:: ../../_static/tbid_bax_activation_8.png
    :width: 6in

As expected, the tBid distribution is perfectly Poissonian. Note that the fraction
of liposomes with no tBid appears identical to the fraction of liposomes with
no activated Bax, as expected.

.. ipython::

    In [5]: plot_hist_vs_poisson(jobs, data, 0, 'mtBid', -1);

    @suppress
    In [6]: savefig('_static/tbid_bax_activation_9.png')

.. image:: ../../_static/tbid_bax_activation_9.png
    :width: 6in

Low tBid concentration, reversible binding
------------------------------------------

Now the question is whether moderate off-rates for tBid can cause Bax to
redistribute more evenly across the liposomes even when tBid is very low.  We
set up the model to have 1 nM tBid, with an off rate of 0.01 per second:

.. ipython::

    In [2]: job = jobs[1]

    In [4]: one_cpt_model = job.one_cpt_builder().model

    In [5]: one_cpt_model.parameters['tBid_transloc_kr']

    In [5]: one_cpt_model.parameters['tBid_0']

There appears to be good correspondence between the stochastic and deterministic
models at the level of the bulk observables:

.. ipython::

    In [6]: plot_timecourse_comparison(jobs, data, 1);

    @suppress
    In [10]: savefig('_static/tbid_bax_activation_10.png')

.. image:: ../../_static/tbid_bax_activation_10.png
    :width: 6in

Remarkably, the moderate reverse rate of 0.01 substantially alleviated
the uneven distribution of Bax:

.. ipython::

    In [8]: plot_hist_vs_poisson(jobs, data, 1, 'iBax', -1);

    @suppress
    In [9]: savefig('_static/tbid_bax_activation_11.png')

.. image:: ../../_static/tbid_bax_activation_11.png
    :width: 6in

And as before, the distribution of tBids is independent:

.. ipython::

    In [5]: plot_hist_vs_poisson(jobs, data, 1, 'mtBid', -1);

    @suppress
    In [6]: savefig('_static/tbid_bax_activation_12.png')

.. image:: ../../_static/tbid_bax_activation_12.png
    :width: 6in

Conclusion
----------

Irreversible, or nearly irreversible, binding of tBid to membranes, is a key
factor in determining how evenly activated Bax is distributed across membranes.
Irreversibility appears to be the key factor even when taking into account tBid
concentration; when tBid concentrations are low, however, the larger fraction
of liposomes with no tBid at all cause strong bimodality in the distribution
of activated Bax. Interestingly, the effect of variation in tBid concentration
is manifest even when tBid concentrations are larger, showing that the uneven
distribution of active Bax is not merely a matter of the fraction of liposomes
that have zero tBid.

Todos
-----

.. todo:: Using tBid irreversibility to make predictions

    Can the near-irreversibility of tBid be used to determine the number of
    tBids required to activate Bax? Is there a predicted tBid sensitivity curve
    once tBid drops below less than one per vesicle? Does this predicted curve
    differs between the stochastic and deterministic models and 

.. todo:: Measure tBid turnover and saturation on liposomes by TIRF microscopy

    Bleach the liposomes after tBid binding.

.. todo:: Repeat tBid binding expt on Octet

    Can I get tBid concentration high enough to see saturation?




