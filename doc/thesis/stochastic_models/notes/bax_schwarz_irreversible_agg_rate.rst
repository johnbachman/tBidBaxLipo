Simple model with irreversible permeabilization, aggregation
============================================================

.. ipython:: python

    from tbidbaxlipo.plots.stoch_det_comparison.bax_schwarz_irreversible_agg_rate \
             import jobs, data
    m = jobs[0].one_cpt_builder().model
    [r.name for r in m.rules]
    ["%s: %f" % (p.name, p.value) for p in m.parameters]

.. plot::
    :context:

    plt.close('all')
    from tbidbaxlipo.plots.stoch_det_comparison.bax_schwarz_irreversible_agg_rate \
            import jobs, data
    from tbidbaxlipo.plots.stoch_det_comparison.plots import *
    plot_bax_timecourse_comparison(jobs, data, 0)

.. plot::
    :context:

    plt.close('all')
    plot_dye_release_titration(jobs, data)

.. plot::
    :context:

    plt.close('all')
    plot_hist_vs_poisson(jobs, data, 0, 'pores', -1)

.. plot::
    :context:

    plt.close('all')
    plot_hist_vs_poisson(jobs, data, 0, 'pBax', -1)


