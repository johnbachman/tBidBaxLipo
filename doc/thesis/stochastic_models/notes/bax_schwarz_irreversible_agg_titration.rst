Simple model with irreversible permeabilization, aggregation, titration
=======================================================================

.. ipython:: python

    from tbidbaxlipo.plots.stoch_det_comparison.bax_schwarz_irreversible_agg_titration \
             import jobs, data
    m = jobs[0].one_cpt_builder().model
    [r.name for r in m.rules]
    ["%s: %f" % (p.name, p.value) for p in m.parameters]

.. plot::
    :context:

    from tbidbaxlipo.plots.stoch_det_comparison.bax_schwarz_irreversible_agg_titration \
            import jobs, data
    from tbidbaxlipo.plots.stoch_det_comparison.plots import *
    plt.close('all')
    plot_dye_release_titration(jobs, data)
    plot_dr_fmax_vs_bax(jobs, data)

.. plot::
    :context:

    plt.close('all')
    plot_bax_timecourse_comparison(jobs, data, 10)

.. plot::
    :context:

    plt.close('all')
    plot_hist_vs_poisson(jobs, data, 10, 'pores', -1)

.. plot::
    :context:

    plt.close('all')
    plot_hist_vs_poisson(jobs, data, 10, 'pBax', -1)



