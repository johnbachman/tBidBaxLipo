Simple translocation model
==========================

We first examine a model for which the deterministic/continuum model
should work well as an approximation for the stochastic system, one
consisting simply of translocation of Bax to membranes. Code
for the comparison is in
:py:class:`tbidbaxlipo.plots.stoch_det_comparison.translocation`.

.. ipython::

    @suppress
    In [1]: from tbidbaxlipo.models import simulation

    @suppress
    In [2]: from tbidbaxlipo.plots.stoch_det_comparison.translocation import *

    @suppress
    In [3]: from tbidbaxlipo.plots.stoch_det_comparison.plot_funcs import *

    @suppress
    In [4]: from tbidbaxlipo.plots.stoch_det_comparison.translocation \
       ...: import bax_kr_nominal, bax_kr_nominal_t06

    In [5]: b_one = bax_kr_nominal.job.one_cpt_builder()

This model has only two rules, for translocation of Bax to and from the
single membrane compartment:

.. ipython::

    In [6]: b_one.model.rules

We simulate both the deterministic and stochastic versions for 60 seconds
(translocation occurs on a fairly fast timescale):

.. plot::

    from tbidbaxlipo.plots.stoch_det_comparison.translocation import *
    from tbidbaxlipo.plots.stoch_det_comparison.translocation.bax_kr_nominal \
         import job, n_cpt_obs
    plot_timecourse_comparison(job, n_cpt_obs)

In this figure, the error bars represent the standard error across all of the
simulations. As expected, the ``n_cpt`` model matches the single compartment
model. In addition, we sample from the distribution of Bax molecules bound to
liposomes (the ``mBax`` observables) at different timepoints and find that the
distribution is indeed Poisson. We plot the distributions averaged across many
simulations to estimate the true frequency of observing a given number of
Bax molecules per liposome; the error bars in the plots indicate the standard
error of these estimates. The values observed in the stochastic simulations are
compared to a Poisson distribution with the value from the deterministic
simulation as its mean. In addition, we can compare the mean and the variance
of the stochastically sampled distribution, since for a Poisson distribution
the mean and variance should be equal:

For the first timepoint:

.. ipython::

    In [1]: job = bax_kr_nominal.job

    In [2]: n_cpt_obs = bax_kr_nominal.n_cpt_obs

    In [3]: mBax_dist = plot_hist_vs_poisson(job, n_cpt_obs, 'mBax', 1);

    In [4]: print_obs_means_and_vars(job, n_cpt_obs, 'mBax', 1)

    @suppress
    In [5]: savefig('_static/simple_translocation_1.png')

.. image:: ../../_static/simple_translocation_1.png
    :width: 6in

For the 20th timepoint:

.. ipython::

    In [6]: mBax_dist = plot_hist_vs_poisson(job, n_cpt_obs, 'mBax', 20);

    In [7]: print_obs_means_and_vars(job, n_cpt_obs, 'mBax', 20)

    @suppress
    In [8]: savefig('_static/simple_translocation_2.png')

.. image:: ../../_static/simple_translocation_2.png
    :width: 6in

For the final timepoint:

.. ipython::

    In [9]: mBax_dist = plot_hist_vs_poisson(job, n_cpt_obs, 'mBax', \
       ...: job.n_steps);

    In [10]: print_obs_means_and_vars(job, n_cpt_obs, 'mBax', job.n_steps)

    @suppress
    In [11]: savefig('_static/simple_translocation_3.png')

.. image:: ../../_static/simple_translocation_3.png
    :width: 6in

