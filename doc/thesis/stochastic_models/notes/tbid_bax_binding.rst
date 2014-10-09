.. _tbid_bax_binding:

Binding of Bax and tBid
=========================

Now we look at a bimolecular system in which tBid and Bax both translocate
to liposomes and then binding to each other.

.. ipython:: python

    from pylab import *
    from tbidbaxlipo.plots.stoch_det_comparison.plots \
        import plot_hist_vs_poisson
    from tbidbaxlipo.plots.stoch_det_comparison.tbid_bax_binding import jobs, \
        data
    from tbidbaxlipo.plots.stoch_det_comparison.tbid_bax_binding.plots import \
        plot_timecourse_comparison
    job = jobs[0]
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

    @savefig tbid_bax_binding_1.png
    In [9]: plot_timecourse_comparison(jobs, data, 0)

We now look at the distribution of membrane bounding Bid and Bax:

.. ipython::

    @savefig tbid_bax_binding_2.png
    In [3]: plot_hist_vs_poisson(jobs, data, 0, 'mBax', -1);

    @savefig tbid_bax_binding_3.png
    In [5]: plot_hist_vs_poisson(jobs, data, 0, 'mtBid', -1);


