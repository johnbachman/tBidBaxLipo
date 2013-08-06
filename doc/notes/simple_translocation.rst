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

    In [3]: from tbidbaxlipo.plots.stoch_det_comparison.translocation import *

    In [5]: b_one = job.build(one_cpt)

This model has only two rules, for translocation of Bax to and from the
single membrane compartment:

.. ipython::

    In [6]: b_one.model.rules

We simulate both the deterministic and stochastic versions for 60 seconds
(translocation occurs on a fairly fast timescale):

.. plot::

    from tbidbaxlipo.plots.stoch_det_comparison.translocation import *
    plot_timecourse_comparison()


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

    In [1]: mBax_dist = plot_hist_vs_poisson('mBax', 1);

    In [2]: print_mbax_means_and_vars(1)

    @suppress
    In [3]: savefig('_static/mbax_dist_t1.png')

.. image:: ../_static/mbax_dist_t1.png
    :width: 6in

For the 20th timepoint:

.. ipython::

    In [4]: mBax_dist = plot_hist_vs_poisson('mBax', 20);

    In [5]: print_mbax_means_and_vars(20)

    @suppress
    In [6]: savefig('_static/mbax_dist_t20.png')

.. image:: ../_static/mbax_dist_t20.png
    :width: 6in

For the final (100th) timepoint:

.. ipython::

    In [7]: mBax_dist = plot_hist_vs_poisson('mBax', job.n_steps);

    In [8]: print_mbax_means_and_vars(job.n_steps)

    @suppress
    In [9]: savefig('_static/mbax_dist_t_final.png')

.. image:: ../_static/mbax_dist_t_final.png
    :width: 6in

